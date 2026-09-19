# VCC 2026 — Submission Builder

A minimal repository for building a Virtual Cell Challenge **2026** submission
from the 2026 control-cell data. It contains marimo notebook that
reads the data, explores it, builds a baseline prediction, validates it against
the submission spec, and writes the submission `.h5ad` file.


## Setup

1. Copy the environment template and edit the paths:

   ```bash
   cp .env.example .env
   ```

   - `DATA_DIR` — directory containing the `vcc_2026/` data (see below).
   - `PROJECT` — path to this repository root (predictions are written under
     `$PROJECT/results/`).

2. Install dependencies with [Pixi](https://pixi.sh):

   ```bash
   pixi install
   ```

## How the data is read

The notebook reads paths from `.env` (via `python-dotenv`) and loads everything
from `$DATA_DIR/vcc_2026/`:

| File | What it is |
|------|------------|
| `gene_names.csv` | 18,533 rows, **has** a header `gene_name`. Defines the submission gene order. Read with polars. |
| `pert_counts.csv` | 300 rows, header `target_gene`. The 300 target genes to predict. Read with polars. |
| `context_{A,B,C}.h5ad` | 18,400 control (non-targeting) cells × 18,533 genes each. `adata.X` holds raw counts; `adata.var_names` is the gene panel. Read with `anndata.read_h5ad`. |

## Normalization of scRNA-seq data

Raw single-cell RNA-seq counts are dominated by per-cell **sequencing depth**
(total UMIs), not biology: a more deeply-sequenced cell shows larger counts for
essentially every gene. Comparing cells directly therefore compares depth, not
expression.

**Library-size normalization** removes that: divide each cell by its total count
and multiply by a fixed target sum, so every cell sums to the same value and its
profile becomes a per-gene proportion. A common form is CP10K (counts per
10,000); many pipelines then apply `log1p`.

What this notebook does (in the profile cell): scales each control cell to a
fixed target sum of `5e4` (`_counts / _depth * 5 * 1e4`), then takes the median
of the 1,000 lowest-expressing control cells per target gene. There is **no
`log1p`** — the 2026 scorer works in counts space, so normalization stays on the
linear count scale and the final profile is rounded to integer counts
(`np.rint`) before being written to the submission.

## What data we predict on

Predictions must cover every combination of:

- **300 target genes** (from `pert_counts.csv`), ×
- **3 contexts** (A, B, C), ×
- **400 cells** each

= **360,000 cells × 18,533 genes**.

The model input is drawn from that context's `context_*.h5ad` control cells. The
baseline in this repo is a placeholder: for each target gene it takes the median
profile of the 1,000 lowest-expressing control cells (scaled to the context's
median depth). Replace that cell with a real model to improve the score.

## Submission format

One AnnData `.h5ad` file:

- `X` — raw **integer** counts, shape `(360000, 18533)`.
- `obs` — two columns: `target_gene` and `context`.
- `var` — indexed by gene name, in `gene_names.csv` order.

Enforced constraints (all asserted in the notebook's validation cell, which
fails loudly if any break):

- counts are whole, non-negative, finite;
- per-cell total ≤ 1,000,000;
- stored non-zero entries ≤ 4.75 × 10⁹;
- no explicitly-stored zeros (they count against the cap);
- exactly 400 cells per (context, target_gene).

## About marimo notebooks

The notebook is a [marimo](https://marimo.io) notebook, stored as a plain `.py`
file (not `.ipynb`), so it diffs and version-controls like normal Python. A
marimo notebook is **reactive**: cells form a dependency graph, and changing one
cell automatically re-runs the cells that depend on it — there is no stale
hidden state and no need to remember run order.

Two ways to use it:

- `marimo edit <file>` — open the interactive editor in the browser to read,
  change, and step through cells.
- `marimo run <file>` — execute it top to bottom (headless / as an app), no
  editing.

You can also run it as an ordinary script (`python <file>`) because of the
`if __name__ == "__main__": app.run()` block at the bottom.

## Run

```bash
pixi run marimo edit notebook/preliminary-data-analysis-2026.py
# or, to run headless:
pixi run marimo run notebook/preliminary-data-analysis-2026.py
```

The validation cell prints `OK  360,000 cells x 18,533 genes ...` on success and
writes `results/2026/prediction_2026_<ddmmyy>.h5ad`.

## Prep and submit

Use the `vcc` CLI (from `vcc-cli`), **not** `cell-eval prep` — the latter is
2025-era and mangles 2026 input (wrong gene dimension, header handling, and it
log-normalises the counts the scorer needs).

```bash
pixi run vcc prep -i results/2026/<file>.h5ad \
    -g $DATA_DIR/vcc_2026/gene_names.csv \
    --perts $DATA_DIR/vcc_2026/pert_counts.csv --context-col context
pixi run vcc submit results/2026/<file>.prep.vcc
```
