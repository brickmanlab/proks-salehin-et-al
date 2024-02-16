---
title: proks-salehin-et-al-2023
date: "20230623"
author: Martin Proks
---

Pre-implanation AI models

## Checklist before submitting

- [ ] Were the data uploaded to the GEO?
- [ ] Were the processed data uploaded to Zenodo? Paste DOI here
- [ ] Did you create `.gitignore` file?
- [ ] Did you create `gh-pages` for deployment?
- [ ] Push the code and run `quarto publish gh-pages`

## Installation

```bash
mamba env create -p /home/fdb589/projects/data/Brickman/conda/envs/scvi-1.0.0 -f environment.yml
```

Package `f2` has to be installed manually, follow [these steps](https://github.com/bhargavchippada/forceatlas2/issues/34#issuecomment-1102409914).

**Note**: If you have issues with cuda make sure the loaded modules are on the same version
as the pip `jaxlib` and `nvidia-cudnn-cu11`.

## Running

```bash
module load miniconda/latest cuda/11.8 cudnn/8.6.0
source activate brickman

jupyter-lab --no-browser
```

## Validation datasets

- Mouse
  - [GSE134240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134240)