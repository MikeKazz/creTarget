"""
executor.py

Runs a batch of independent shell commands ("chunks of work") either:
  - locally, in parallel, via a process pool (default, --executor local)
  - as a SLURM array job (--executor slurm)

This generalizes the pattern already proven for launch-dCas9 scoring
(chunk -> submit array -> wait -> verify) to any chunked work, starting
with FlashFry discover/score.

Each "chunk" is just a shell command string. The executor doesn't know or
care what the command does -- it only handles running it, waiting, and
reporting success/failure per chunk.
"""

import concurrent.futures
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ChunkResult:
    index: int
    label: str
    returncode: int
    log_path: str


class LocalExecutor:
    """Runs chunks in parallel local subprocesses via a thread pool (each
    chunk is itself an OS process via subprocess.run, so this parallelizes
    real work despite Python's GIL -- the threads are just supervising
    subprocesses, not doing CPU work themselves)."""

    def __init__(self, n_workers: int | None = None):
        self.n_workers = n_workers or max(1, (os.cpu_count() or 4) - 1)

    def run_batch(self, commands: list[str], logdir: str, labels: list[str] | None = None) -> list[ChunkResult]:
        Path(logdir).mkdir(parents=True, exist_ok=True)
        labels = labels or [f"chunk_{i:04d}" for i in range(len(commands))]
        results: list[ChunkResult] = [None] * len(commands)  # type: ignore

        def _run_one(i, cmd, label):
            log_path = os.path.join(logdir, f"{label}.log")
            with open(log_path, "w") as logf:
                proc = subprocess.run(cmd, shell=True, stdout=logf, stderr=subprocess.STDOUT)
            return ChunkResult(index=i, label=label, returncode=proc.returncode, log_path=log_path)

        print(f"  [LocalExecutor] running {len(commands)} chunk(s) with {self.n_workers} parallel workers", file=sys.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.n_workers) as ex:
            futures = {
                ex.submit(_run_one, i, cmd, labels[i]): i
                for i, cmd in enumerate(commands)
            }
            done = 0
            for fut in concurrent.futures.as_completed(futures):
                r = fut.result()
                results[r.index] = r
                done += 1
                status = "OK" if r.returncode == 0 else f"FAILED (exit {r.returncode})"
                print(f"  [LocalExecutor] {done}/{len(commands)} {r.label}: {status}", file=sys.stderr)

        return results


class SlurmExecutor:
    """Submits an array job where each task runs one chunk's command (each
    command is written to its own tiny shell script, dispatched by array
    index). Waits via squeue polling, then verifies REAL per-task exit
    codes via sacct -- squeue only tells you a task left the queue, not
    whether it succeeded, and a log file existing is not evidence of
    success either (SLURM creates it the moment the task starts, so a
    segfault/OOM/walltime-kill still leaves a log file behind)."""

    def __init__(self, partition="", gres="", mem="16G", cpus="8", time_="2:00:00",
                 poll_interval=30, max_consecutive_poll_failures=10,
                 sacct_max_retries=5, sacct_retry_delay=5):
        self.partition = partition
        self.gres = gres
        self.mem = mem
        self.cpus = cpus
        self.time = time_
        self.poll_interval = poll_interval
        self.max_consecutive_poll_failures = max_consecutive_poll_failures
        self.sacct_max_retries = sacct_max_retries
        self.sacct_retry_delay = sacct_retry_delay

    def run_batch(self, commands: list[str], logdir: str, labels: list[str] | None = None) -> list[ChunkResult]:
        logdir_path = Path(logdir)
        logdir_path.mkdir(parents=True, exist_ok=True)
        labels = labels or [f"chunk_{i:04d}" for i in range(len(commands))]

        task_dir = logdir_path / "tasks"
        task_dir.mkdir(exist_ok=True)
        for i, cmd in enumerate(commands):
            task_script = task_dir / f"task_{i:04d}.sh"
            task_script.write_text(f"#!/bin/bash\nset -euo pipefail\n{cmd}\n")
            task_script.chmod(0o755)

        array_range = f"0-{len(commands) - 1}"
        sbatch_cmd = ["sbatch", "--parsable", "--array", array_range,
                      "--mem", self.mem, "--cpus-per-task", self.cpus, "--time", self.time,
                      "-o", str(logdir_path / "%a.log"), "-e", str(logdir_path / "%a.err")]
        if self.partition:
            sbatch_cmd += ["-p", self.partition]
        if self.gres:
            sbatch_cmd += ["--gres", self.gres]

        wrapper = task_dir / "dispatch.sh"
        wrapper.write_text(
            "#!/bin/bash\n"
            "set -euo pipefail\n"
            f'TASK_ID=$(printf "%04d" "${{SLURM_ARRAY_TASK_ID}}")\n'
            f'bash "{task_dir}/task_${{TASK_ID}}.sh"\n'
        )
        wrapper.chmod(0o755)
        sbatch_cmd.append(str(wrapper))

        print(f"  [SlurmExecutor] submitting array job ({array_range}): {' '.join(sbatch_cmd)}", file=sys.stderr)
        submit = subprocess.run(sbatch_cmd, capture_output=True, text=True, check=True)
        job_id = submit.stdout.strip()
        print(f"  [SlurmExecutor] submitted job {job_id}, waiting...", file=sys.stderr)

        self._wait_for_completion(job_id)

        print(f"  [SlurmExecutor] job {job_id} no longer in queue, verifying real exit codes via sacct...", file=sys.stderr)
        states = self._sacct_states(job_id, len(commands))

        results = []
        for i, label in enumerate(labels):
            log_path = str(logdir_path / f"{i}.log")
            if i not in states:
                print(f"  [SlurmExecutor] WARNING: no sacct record found for task {i} "
                      f"(job {job_id}_{i}) after {self.sacct_max_retries} attempts -- "
                      f"treating as FAILED rather than assuming success", file=sys.stderr)
                results.append(ChunkResult(index=i, label=label, returncode=1, log_path=log_path))
                continue

            state, exitcode = states[i]
            success = (state == "COMPLETED")
            rc = 0 if success else 1
            status_str = "OK" if success else f"FAILED (state={state}, exitcode={exitcode})"
            print(f"  [SlurmExecutor] task {i} ({label}): {status_str}", file=sys.stderr)
            results.append(ChunkResult(index=i, label=label, returncode=rc, log_path=log_path))

        return results

    def _wait_for_completion(self, job_id: str):
        """Polls squeue until the job leaves the queue. Treats a squeue
        FAILURE (nonzero returncode -- e.g. a transient controller hiccup,
        common on busy shared clusters) as "still waiting, try again", NOT
        as "job is done" -- empty/failed squeue output must not be
        confused with an empty *successful* result."""
        consecutive_failures = 0
        while True:
            pending = subprocess.run(["squeue", "-j", job_id, "-h"], capture_output=True, text=True)

            if pending.returncode != 0:
                consecutive_failures += 1
                print(f"  [SlurmExecutor] squeue failed (exit {pending.returncode}, attempt "
                      f"{consecutive_failures}/{self.max_consecutive_poll_failures}): "
                      f"{pending.stderr.strip()}", file=sys.stderr)
                if consecutive_failures >= self.max_consecutive_poll_failures:
                    raise RuntimeError(
                        f"squeue failed {self.max_consecutive_poll_failures} times in a row while "
                        f"waiting for job {job_id}. This does NOT necessarily mean the job failed -- "
                        f"it means we couldn't check. Verify job {job_id}'s status manually (squeue/sacct) "
                        f"before assuming anything about its outcome."
                    )
                time.sleep(self.poll_interval)
                continue

            consecutive_failures = 0  # reset only on a genuinely successful squeue call

            if not pending.stdout.strip():
                return  # job has left the queue

            n = len(pending.stdout.strip().splitlines())
            print(f"  [SlurmExecutor] {n} task(s) still running, checking again in {self.poll_interval}s", file=sys.stderr)
            time.sleep(self.poll_interval)

    def _sacct_states(self, job_id: str, n_tasks: int) -> dict:
        """Returns {task_index: (state, exitcode_str)} via sacct, retrying
        on transient sacct failures or incomplete accounting records (which
        can lag briefly after squeue shows the job gone -- the accounting
        DB update is not synchronous with the job leaving the queue)."""
        states: dict = {}
        for attempt in range(1, self.sacct_max_retries + 1):
            proc = subprocess.run(
                ["sacct", "-j", job_id, "--parsable2", "-n", "-o", "JobID,State,ExitCode"],
                capture_output=True, text=True,
            )
            if proc.returncode != 0:
                print(f"  [SlurmExecutor] sacct failed (attempt {attempt}/{self.sacct_max_retries}): "
                      f"{proc.stderr.strip()}", file=sys.stderr)
                time.sleep(self.sacct_retry_delay)
                continue

            states = {}
            for line in proc.stdout.strip().splitlines():
                parts = line.split("|")
                if len(parts) != 3:
                    continue
                task_jobid, state, exitcode = parts
                # Array task records look like "<job_id>_<i>"; sub-steps
                # ("<job_id>_<i>.batch", ".extern") are not what we want --
                # only the task-level record reflects the actual command's
                # real exit status.
                if "." in task_jobid or "_" not in task_jobid:
                    continue
                prefix, _, idx_str = task_jobid.partition("_")
                if prefix != job_id:
                    continue
                try:
                    idx = int(idx_str)
                except ValueError:
                    continue
                states[idx] = (state, exitcode)

            if len(states) >= n_tasks:
                return states

            print(f"  [SlurmExecutor] sacct returned {len(states)}/{n_tasks} task record(s) so far "
                  f"(attempt {attempt}/{self.sacct_max_retries}), retrying...", file=sys.stderr)
            time.sleep(self.sacct_retry_delay)

        return states  # best effort -- missing entries are handled by the caller as failures


def get_executor(kind: str, **kwargs):
    if kind == "local":
        return LocalExecutor(n_workers=kwargs.get("n_workers"))
    elif kind == "slurm":
        return SlurmExecutor(
            partition=kwargs.get("partition", ""),
            gres=kwargs.get("gres", ""),
            mem=kwargs.get("mem", "16G"),
            cpus=kwargs.get("cpus", "8"),
            time_=kwargs.get("time", "2:00:00"),
        )
    else:
        raise ValueError(f"Unknown executor kind: {kind}")
