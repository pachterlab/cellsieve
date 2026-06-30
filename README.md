# Cell Sieve
Combinatorially compressed gene panels to distinguish cell types in scRNAseq data.

## Install

```bash
conda env create -f environment.yml
conda activate cellsieve
python -m pip install -e .
```

## Usage

Cell Sieve selects compact gene panels for cell type classification from
single-cell RNA-seq data using an Elastic Net prior and a Gurobi-backed MILP
solver.

```bash
cellsieve --adata-path DATA.h5ad [options]
```

| Option | Default | What it does |
| --- | --- | --- |
| `--adata-path` | Required | Path to the input AnnData `.h5ad` file. |
| `--celltype-column` | `celltype` | Column in `adata.obs` containing cell type labels. |
| `--out-dir` | `./elastic_sieve_output` | Directory where selected genes, reports, metrics, and optional corrected data are written. |
| `--random-seed` | `42` | Random seed used for stratified subsampling, model fitting, and CV splits. |
| `--alpha` | `0.7` | Elastic Net L1 ratio. Larger values emphasize sparsity; smaller values behave more like ridge regularization. |
| `--panel-size` | `20` | Number of genes to select when not using `--titrate-panels`. |
| `--titrate-panels` | Not set | Run multiple panel sizes in sequence, for example `--titrate-panels 5 10 20 40`. |
| `--max-cells` | `15000` | Maximum number of cells used to fit the global Elastic Net prior. Larger datasets are stratified downsampled. |
| `--correct-misannotations` | Off | Uses out-of-fold predictions from the selected panel to flag likely annotation errors, writes a report, and saves a corrected `.h5ad`. |
| `--min-correction-margin` | `0.30` | Minimum probability margin required before changing an annotation during correction. |
| `--find-redundant-swaps` | Off | Finds highly correlated substitute genes for selected panel genes. |
| `--corr-threshold` | `0.75` | Minimum Pearson correlation for reporting redundant gene swaps. |
| `--gurobi-license-file` | Not set | Optional path to `gurobi.lic`; sets `GRB_LICENSE_FILE` for the current run. |
| `--mip-gap` | `0.10` | Relative MIP optimality gap tolerance for Gurobi. |
| `--time-limit` | `300` | Gurobi time limit in seconds for each panel-size solve. |
| `--n-workers` | `24` | Number of worker threads/jobs used by scikit-learn, joblib, and Gurobi. |

## Examples

Run the default workflow:

```bash
cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --n-workers 8
```

Select a custom single panel size:

```bash
cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --panel-size 40
```

Run a panel-size titration:

```bash
cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --titrate-panels 5 15 30
```

Flag likely misannotations and save a corrected AnnData file:

```bash
cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --correct-misannotations --panel-size 20
```

Find correlated substitute genes for the selected panel:

```bash
cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --find-redundant-swaps --corr-threshold 0.80
```

## Gurobi license

Keep your `gurobi.lic` file outside the repository. Pass the path at runtime:

```bash
cellsieve --gurobi-license-file /path/to/gurobi.lic --adata-path data.h5ad
```
You can also set `GRB_LICENSE_FILE` in your shell instead of using the flag:

```bash
export GRB_LICENSE_FILE=/path/to/gurobi.lic
cellsieve --adata-path data.h5ad
```

Academic users can request a free Gurobi license through the official
[Gurobi Academic Program](https://www.gurobi.com/academics). Gurobi also
maintains a step-by-step Help Center article:
[How do I obtain a free academic license?](https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license)

For local machine use, look for the Academic Named-User License. For cloud,
container, or multi-machine workflows, Gurobi's Academic WLS License may be a
better fit.
