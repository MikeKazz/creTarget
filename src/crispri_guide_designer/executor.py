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
    index). Waits via squeue polling, same pattern already proven for
    launch-dCas9."""

    def __init__(self, partition="", gres="", mem="16G", cpus="8", time_="2:00:00", poll_interval=30):
        self.partition = partition
        self.gres = gres
        self.mem = mem
        self.cpus = cpus
        self.time = time_
        self.poll_interval = poll_interval

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
        job_id = subprocess.run(sbatch_cmd, capture_output=True, text=True, check=True).stdout.strip()
        print(f"  [SlurmExecutor] submitted job {job_id}, waiting...", file=sys.stderr)

        while True:
            pending = subprocess.run(["squeue", "-j", job_id, "-h"], capture_output=True, text=True)
            if not pending.stdout.strip():
                break
            n = len(pending.stdout.strip().splitlines())
            print(f"  [SlurmExecutor] {n} task(s) still running, checking again in {self.poll_interval}s", file=sys.stderr)
            time.sleep(self.poll_interval)

        results = []
        for i, label in enumerate(labels):
            log_path = str(logdir_path / f"{i}.log")
            # SLURM %a expands to array task id, not zero-padded -- our logs are named by raw index
            raw_log = logdir_path / f"{i}.log"
            returncode = 0 if raw_log.exists() else 1
            results.append(ChunkResult(index=i, label=label, returncode=returncode, log_path=str(raw_log)))
        return results


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
