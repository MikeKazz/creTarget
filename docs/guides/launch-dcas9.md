# Optional: launch-dCas9 on-target scoring

`--run_launch_dcas9` adds two CNN-based efficacy predictions to every guide:

- `launch_dCas9_enhancer` — enhancer fitness
- `launch_dCas9_promoter` — promoter fitness

The enhancer score contributes **30%** of the composite ranking score when
present. The promoter score is reported for reference and is intentionally
**not** part of the ranking. See
[Output files](outputs.md#composite_score).

The model and the inference script are from
[launch-dCas9](https://github.com/Wancen/launch-dCas9). creTarget does not
retrain or modify them.

## What you need

1. The creTarget package and pipeline (see [Installation](../getting-started/installation.md)).
2. A checkout of launch-dCas9 with its `launch-dCas9.py` inference script.
3. A Python environment with `torch` that can run that script.
4. The two pretrained checkpoints, fetched by creTarget.

### Fetch the checkpoints

```bash
cretarget fetch-models
```

This downloads `ckpt-seq-enhancerFitness.pth` and
`ckpt-seq-promoterFitness.pth` from Zenodo record
[10.5281/zenodo.22729308](https://doi.org/10.5281/zenodo.22729308), verifies
their SHA-256 hashes, and prints the attribution. They land in
`~/.cretarget/models` by default. They are **not** in the Python package or the
container image and remain the work of the launch-dCas9 authors.

### Prepare the inference environment

Follow the launch-dCas9 repository's instructions for its Python environment.
creTarget needs two paths from you:

- `--launch_script`: the path to `launch-dCas9.py`
- `--launch_python`: the Python interpreter from the torch environment

creTarget itself runs under its own interpreter; only `launch-dCas9.py` is
invoked with `--launch_python`. That means you do **not** need to install
creTarget into the torch environment.

```bash
git clone https://github.com/Wancen/launch-dCas9.git
# Create the environment as described in that repository, then:
/path/to/torch/env/bin/python -c "import torch; print(torch.__version__)"
```

## Run

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile slurm,conda \
    --enhancer_bed enhancers.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --run_launch_dcas9 \
    --launch_script /path/to/launch-dCas9.py \
    --launch_python /path/to/torch/env/bin/python \
    --ckpt_dir ~/.cretarget/models \
    --outdir results
```

If the dtype or CUDA configuration differs between your environment and the
model, check the launch-dCas9 repository; creTarget passes the same
`--model CNN --variant seq` interface their script documents.

### GPU requests

- **SLURM:** `--gpu_gres gpu:1` adds `--gres=gpu:1` to the launch-dCas9 task.
  The scheduler sets `CUDA_VISIBLE_DEVICES`; creTarget does not override it.
- **Local executor:** `--gpu_count N` pins concurrent tasks round-robin across
  `N` devices. Pair it with `maxForks` so two tasks do not share a GPU.
- **CPU-only:** works, just more slowly. No `--gpu_*` flag is required.

## How it runs

launch-dCas9 scoring is per chunk. Within each chunk, creTarget deduplicates
spacers, scores each unique sequence once, and maps predictions back to every
row. A spacer that appears in two chunks is scored twice; the prediction is a
property of the sequence, so the result is identical.

Rows with non-ATGC characters in the spacer cannot be one-hot encoded and get
`NaN` scores rather than crashing the step.

## Output columns

| Column | Meaning |
|---|---|
| `launch_dCas9_enhancer` | Predicted enhancer fitness |
| `launch_dCas9_promoter` | Predicted promoter fitness |

Both are added to `guides.ranked.tsv` and carried into `guides.final.tsv` when
seed screening is enabled.

## Citation

If you use this step, cite launch-dCas9 in addition to creTarget. The
checkpoints are redistributed with attribution; they remain under the original
authors' terms. See [Citation](../project/citation.md).

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `checkpoint(s) not found` | Run `cretarget fetch-models`, or pass the correct `--ckpt_dir` |
| `launch-dCas9.py failed` | Check `--launch_python`, the torch environment, and the model path |
| `expected prediction output not found` | The launch-dCas9 version or output path differs; check its CLI |
| Slow CPU-only run | Expected. Request a GPU or accept the runtime |
| CUDA out of memory | Lower concurrency (`maxForks`), use `--gpu_count`, or use a larger GPU |
