# Neuwirth_Malzl_et_al_2024
[![DOI](https://zenodo.org/badge/698192309.svg)](https://zenodo.org/doi/10.5281/zenodo.10849496)

This repository contains the code for all bioinformatic analyses shown in Neuwirth and Malzl et al. 2024. 
Most code is presented in the `notebooks/sctools` package and should easily be executable after installing the accompanying conda environment
given in `envs/scpython.yml`.

## Installing the environments
Installation of the prvided environments is illustrated by the scpython environment. To install an environment you will first need to install a variant of Anaconda. We recommend [`miniconda`](https://docs.anaconda.com/free/miniconda/) for this.
After this you can simply type
```
conda env create -f envs/scpython.yml
```
This will prompt conda to install all necessary packages as they were used for the presented analyses. Please note that, to properly use the environment with [`JupyterLab`](https://jupyter.org/)
you will have to install the ipy-kernel which is done by typing
```
conda activate scpython
python -m ipykernel install --user --name scpython
```
The table below details which environments in the enclosed `envs` folder was used for which parts of the analyses
|envrionment|part of the analyses|
|-----------|--------------------|
|scpython|general scRNA-seq analysis like integration and cell type annotation|
|scenic  |gene regulatory network inference of skin diseases Tregs            |
|recombat|batch effect correction prior to SCENIC GRN inference               |

Note on scpython: To run Milo analyses with Python. milopy has to be installed manually from github source as described [here](https://github.com/emdann/milopy). We used v0.1.0

Please note that the presented cell trajectory analysis was conducted in a separate R environment with R 4.3.0 and Bioconductor 3.17 in which Monocle3 1.3.7 was installed as described [here](https://cole-trapnell-lab.github.io/monocle3/docs/installation/)
## Note on running pySCENIC
pySCENIC seems to be ill maintained (at least when it comes to the PyPI packages). At the time of this writing both [pySCENIC](https://github.com/aertslab/pySCENIC) and the [arboreto](https://github.com/aertslab/arboreto) packages used by pySCENIC did not run when installing them from their PyPI distribution. This was due to the version not being bumped after fixing the encountered bugs which prevents the update of the PyPI packages. Therefore, I recommend to install those packages from the clones git repositories for the code to run properly.

## Note on the bundled `sctools` library
The bundled `sctools` is there for reproducibility reasons and corresponds to the snapshot of the code used throughout this study. Furthermore, it is only usable when located in the current working directory (which means notebooks have to reside in `notebooks` in order to work properly. In case you are looking for a more maintained and easier to work with library consider the [`sctools` repository](https://github.com/menchelab/sctools).
