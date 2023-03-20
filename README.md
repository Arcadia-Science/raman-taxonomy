# raman-taxonomy

This repository contains analysis code to accompany the publication 'Identifying the phylogenetic utility of Raman spectra'.<br>

Notebooks containing expanded methods and the code + analyses for generating all figures in the publication can be found in [this Jupyter notebook](01_notebooks/notebook-1.ipynb).

This repository uses conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the run environment.

```
mamba env create --file environment.yml
mamba activate raman-taxonomy
```

Following this, you can run the Jupyter notebook as usual.

**Note:** The conda environment has been tested with a Mac OSX running the Intel x64 architecture. While most packages should be available on Linux distributions, not all packages are available on the new Apple M1/M2 ARM64 architecture.

## Directory structure

[00_data/](00_data/) Raman spectra and associated metadata.<br>
[01_notebooks/](01_notebooks/) Jupyter notebooks outlining the full suite of analyses presented in the pub.<br>
[02_outputs/](02_outputs/) Outputs associated with the pub.<br>

## Data

TODO: Add a description of how you acquired the NCBI data, Ho et al data and the timetree data.
TODO: For the Ho et al, data maybe talk about why this dataset and its relevance/importance.

## R packages and versions used in this repo:

`TreeDist v2.6.0` <br/>
`ArcadiaColorBrewer@c6cea41` <br/>
`scales v1.2.1` <br/>
`phytools v1.2-0` <br/>
`ape v3.4.1` <br/>
`lsa v0.73.3` <br/>
`plotrix v3.8-2` <br/>
`reticulate v1.27` <br/>
`here v1.0.1` <br/>

![github_front_page_figure](fig_1.png)
