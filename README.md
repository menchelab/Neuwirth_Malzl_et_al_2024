# Neuwirth_Malzl_et_al_2024
[![DOI](https://zenodo.org/badge/698192309.svg)](https://zenodo.org/doi/10.5281/zenodo.10849496)
This repository contains the code for all bioinformatic analyses shown in Neuwirth and Malzl et al. 2024. 
Most code is presented in the `notebooks/sctools` package and should easily be executable after installing the accompanying conda environment
given in `environment.yml`.

## Installing the environment
To install the environment you will first need to install a variant of Anaconda. We recommend [`miniconda`](https://docs.anaconda.com/free/miniconda/) for this.
After this you can simply type
```
conda env create -f environment.yml
```
This will prompt conda to install all necessary packages as they were used for the presented analyses. Please note that, to properly use the environment with [`JupyterLab`](https://jupyter.org/)
you will have to install the ipy-kernel which is done by typing
```
conda activate scpython
python -m ipykernel install --user --name scpython
```
