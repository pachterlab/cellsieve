# cellsieve
Cell Sieve: Selecting minimal cell type marker panels

## Install

```bash
conda env create -f environment.yml
conda activate cellsieve
python -m pip install -e . --no-deps --no-build-isolation
```

## Gurobi license

Keep your `gurobi.lic` file outside the repository. Pass the path at runtime:

```bash
cellsieve --gurobi-license-file /path/to/gurobi.lic --adata-path data.h5ad
```

On Nikki's machine, for example:

```bash
cellsieve --gurobi-license-file /home/nikki/gurobi.lic --adata-path data.h5ad
```

You can also set `GRB_LICENSE_FILE` in your shell instead of using the flag:

```bash
export GRB_LICENSE_FILE=/path/to/gurobi.lic
cellsieve --adata-path data.h5ad
```
