# 📄 paper_LORA

[![DOI](https://zenodo.org/badge/1067413189.svg)](https://doi.org/10.5281/zenodo.17250176)

This repository contains the Jupyter notebooks supporting the **LORA (LOng Read Assembly)** paper.  
LORA is part of the Sequana project:

- LORA repository: https://github.com/sequana/lora  
- Sequana project: https://github.com/sequana/sequana

The notebooks reproduce all figures shown in the paper and automatically download the data required for each analysis.

---

## 🚀 Running the notebooks locally (recommended)

Because several analyses require downloading large datasets, we **strongly recommend running notebooks locally** rather than through Binder.

### 1. Create the conda environment

Use the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate paper_lora
```

Install JupyterLab::

```bash
pip install jupyterlab
```

You can open the notebooks interactively with

```bash
jupyter-lab
```


```bash
pip install papermill
papermill --log-output veillonella.ipynb output.ipynb
```

This allows you to monitor progress and identify potential slow steps (e.g., data downloads). This will generate the
figures to be found in the paper.

## ⚠️ A note about Binder

We provide Binder links for convenience, but some notebooks require downloading large datasets, which may exceed Binder’s time or bandwidth limits.

- resources.ipynb and veillonella.ipynb generally run successfully.
- Other notebooks may fail or run very slowly, depending on network conditions.

Binder is therefore best suited for quick inspection rather than full execution.

## 📓 Available notebooks

Each notebook contains the complete workflow used to produce the figures in the paper.

- resources.ipynb – CPU and memory usage - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fresources.ipynb)
- veillonella.ipynb - Assembly, BUSCO, checkm, statitstics - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fveillonella.ipynb)
- leishmania.ipynb - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fleishmania.ipynb) 
- cyanobacteria.ipynb - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fcyanobacteria.ipynb)
- streptococcus.ipynb - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fstreptococcus.ipynb) 
- bacteroides_fragilis.ipynb - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cokelaer/paper_LORA/HEAD?urlpath=%2Fdoc%2Ftree%2Fbacteroides.ipynb) 

