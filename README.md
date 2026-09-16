# LymphNode

Vertex model simulations of 2D lymph node tissue, comparing two mechanical models
of the fibroblastic reticular cell (FRC) network:

- **Model 1** — cell mechanics governed by area elasticity and length elasticity only.
- **Model 2** — adds line tension on top of area and length elasticity.

Simulations cover parameter sweeps, laser ablation with and without collagenase treatment,
tissue homeostasis (edge collapse + vertex division), tissue expansion
with apoptosis and division and tissue expansion with mechanical perturbations.

## Setup

**Note: this environment is Linux-only.** It has not been tested on macOS or Windows.

A pre-built, tested Linux environment is available for download:
https://doi.org/10.5281/zenodo.22705658 (`tyssue_env.tar.gz`, ~741.7 MB)

```bash
mkdir tyssue_env
tar -xzf tyssue_env.tar.gz -C tyssue_env
source tyssue_env/bin/activate
conda-unpack
```
