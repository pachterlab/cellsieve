# Cell Sieve
Combinatorially compressed gene panels to distinguish cell types in scRNAseq data.

## Install

```bash
conda env create -f environment.yml
conda activate cellsieve
python -m pip install -e .
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
