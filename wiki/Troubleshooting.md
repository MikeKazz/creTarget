# Troubleshooting

Start with the diagnostic command:

```bash
cretarget doctor
```

It checks external tools, configured reference paths, the `.fai` and FlashFry
header, cached models, and the seed reference. It fails on missing required
pieces and reports optional ones without failing.

For a failed pipeline run, the most useful files are:

| File | What it tells you |
|---|---|
| `.nextflow.log` | Launch-time errors, parameter validation, profile resolution |
| `results/pipeline_info/trace.txt` | Per-task status, exit code, peak memory, runtime |
| `results/run_manifest.json` | Version, commit, parameters, command line |
| `work/<hash>/.command.err` | The failing task's stderr |
| `work/<hash>/.command.log` | The full task log |

Add `-with-report report.html` for an interactive resource report.

## Installation and launch

### Nextflow reports an unsupported version

creTarget requires **Nextflow 24.04 or newer**. The pipeline uses
`process.resourceLimits`, which earlier versions ignore, so resource ceilings
would silently not apply.

### Java version is rejected

Use Java 17 or newer. Nextflow 23.x caps at Java 19, another reason to use a
current Nextflow.

### `No FASTA index at ...fai`

```bash
samtools faidx /path/to/hg38.fa
```

The pipeline checks before launching, rather than failing inside bedtools.

### `FlashFry index header not found`

`hg38_standard_cas9ngg_db.header` must sit next to `hg38_standard_cas9ngg_db`.
Point `--flashfry_db` at the database file, not the header.

### `Profile ... enables the ... container engine but no image was given`

Pass `--container_image ghcr.io/mikekazz/cretarget:0.5.1`, or drop the container
profile. This is a deliberate hard error: running with an engine enabled and no
image would use host tools while reporting a container.

### `--run_seed_finder requires --seed_finder_min_length`

Those parameters have no defaults on purpose. Choose them deliberately; see
[Seed screening](Seed-Screening).

## Input problems

### `no usable regions parsed`

Every line was blank, a comment, or malformed. Check that the parser matches the
content. A `chr:pos` line in `--enhancer_bed` is skipped, and vice versa.

### `no input regions remain after contig validation`

The input contig names do not match the FASTA `.fai`. The common case is
Ensembl-style `1` against UCSC-style `chr1`, or the reverse. creTarget prints a
hint. Rename the contigs or use a matching reference.

### A SNP BED row has zero length

BED is 0-based, half-open. A one-base SNP at 1-based position 120,000 is
`chr1  119999  120000`, not `chr1  120000  120000`.

### `input file not found or not a regular file`

A process substitution such as `<(zcat peaks.bed.gz)` or a pipe was passed.
Nextflow stages inputs by file name, so the substitution does not survive. Write
a real file first.

### Coordinates look shifted by one

Check the source convention. `chr:pos` is 1-based. `chr:start-end` and BED are
0-based, half-open. See [Inputs](Inputs).

## Runtime problems

### A task is OOM-killed or times out

The pipeline retries with escalated memory and time, up to `maxRetries = 2` and
the `max_*` ceilings. If it still fails:

- Lower `--chunk_size` so each task holds less.
- Raise `--max_memory` / `--max_time` if the partition grants more.
- Check the trace's `peak_rss`. `PREP_REGIONS` grows with the **total** number of
  targets, not chunk size.
- For seed screening, set `--seed_finder_max_guides` to shard the screen.

### Tasks are rejected rather than queued

A request escalated past what the partition allows. Lower `max_cpus`,
`max_memory`, or `max_time` in a site profile to the largest single allocation
the partition can grant.

### `conda activate` fails inside a batch job

A batch job gets a non-interactive shell, so the conda block in `~/.bashrc`
never runs. Source `.../etc/profile.d/conda.sh` first.
`assets/run_cretarget.sbatch` does this.

### Singularity uses a stale cached image

A stale `NXF_SINGULARITY_CACHEDIR` in the environment. Override it to a path
under your data directory.

### `flock` or work-directory errors on NFS

Nextflow's cache needs file locking. NFS without locking fails in confusing
ways. Put the work directory on a filesystem that supports locking, or fix the
mount options. `assets/inspect_cluster.sh` tests this.

### `-resume` does not reuse anything

The cache key includes parameters and input paths. Changing a flag, moving the
work directory, or pointing at a different absolute path invalidates the
affected tasks. Keep the same `-w` and input paths.

### `no guides to screen`

Every guide failed the hard filters. That is a real result. If you intentionally
want to screen them anyway, pass `--seed_finder_all_guides`.

### Out of memory in seed screening

The seed index is ~1 GB of RSS per million guides. Set
`--seed_finder_max_guides` to split the table into smaller tasks. Each shard
repeats the reference pass, so choose the largest shard that fits.

### `checkpoint(s) not found`

Run `cretarget fetch-models`, or pass the correct `--ckpt_dir`. The pipeline
does not download mid-run; this is deliberate for compute nodes without
network.

### `distance_to_atac_summit` is empty / NaN

There were no peaks on that guide's contig in the narrowPeak file. `NaN` is not
zero, and the step reports how many guides were affected.

### Output ordering differs between runs

It should not. Chunks are collated by numeric chunk index, not by work-directory
path, precisely so identical runs produce byte-identical output. If you see a
difference, open an issue with both tables and both `run_manifest.json` files.

## Getting help

- [FAQ](FAQ)
- [Open a bug report](https://github.com/MikeKazz/creTarget/issues/new?template=bug_report.yml)
  with the version, command, `doctor` output, trace, and smallest reproducing
  input.
- For security problems, follow
  [SECURITY.md](https://github.com/MikeKazz/creTarget/blob/main/SECURITY.md),
  not a public issue.

Full guide: [Troubleshooting](https://mikekazz.github.io/creTarget/guides/troubleshooting/)
