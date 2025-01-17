---
title: proks-salehin-et-al
date: "16-02-2024"
author: Martin Proks, Nazmus Salehin
---

Deep Learning Based Models for Preimplantation Mouse and Human Development

- Tutorial: [![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/brickmanlab/proks-salehin-et-al/blob/master/docs/00_tutorial.ipynb)
- Datasets: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10669600.svg)](https://doi.org/10.5281/zenodo.10669600)
- Trained models: [https://huggingface.co/brickmanlab](https://huggingface.co/brickmanlab)
- Portal: [https://brickman-preimplantation.streamlit.app](https://brickman-preimplantation.streamlit.app) 
- Portal codebase: [https://github.com/brickmanlab/preimplantation-portal](https://github.com/brickmanlab/preimplantation-portal)
- bioRxiv paper: [https://www.biorxiv.org/content/10.1101/2024.02.16.580649v1](https://www.biorxiv.org/content/10.1101/2024.02.16.580649v1)

## Installation

```bash
mamba env create -p /home/fdb589/projects/data/Brickman/conda/envs/scvi-1.0.0 -f environment.yml
```

Package `f2` has to be installed manually, follow [these steps](https://github.com/bhargavchippada/forceatlas2/issues/34#issuecomment-1102409914).

**Note**: If you have issues with cuda make sure the loaded modules are on the same version
as the pip `jaxlib` and `nvidia-cudnn-cu11`.

## Running

```bash
module load cuda/11.8-dangpu cudnn/8.6.0-dangpu miniconda/latest && source activate brickman
jupyter-lab --no-browser
```
