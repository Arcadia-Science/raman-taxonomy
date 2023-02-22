# raman-taxonomy

This repository contains analysis code to accompany the publication 'Identifying the phylogenetic utility of Raman spectra'.<br>

Notebooks containing expanded methods and the code + analyses for generating all figures in the publication can be run in Binder: <br>
[![Binder](https://aws-uswest2-binder.pangeo.io/badge_logo.svg)](https://aws-uswest2-binder.pangeo.io/v2/gh/Arcadia-Science/raman-taxonomy/main)

Binder is a service that turns a Git repo into a collection of interactive notebooks that can be executed on a temporary cloud machine.
Binder is currently free, so there is no cost to using it!
We built the Binder for this repository on [Pangeo Binder](https://pangeo-binder.readthedocs.io/en/prod/).
This BinderHub provides slightly more powerful compute, but you have to login via GitHub and can only have one instance running at a time.
GitHub accounts are also free, so you can create one if you don't already have one.

<details>
    <summary>More information on binder and what happens when you click the launch binder button.</summary>
When a repository is configured to run as a binder, passing the GitHub repository URL to binder starts the binder-building process.
Binder first builds a docker image that contains all of the software installations specified by a special set of files in the GitHub repository.
A docker image is a set of instructions that are used to create a docker container.
A docker container is a runnable instance of a docker image -- it's an encapsulated computing environment that can be used to reproducibly install sets of software on diverse computers.
Armed with the docker container, binder launches an "instance" in the cloud (either on Google Cloud or AWS typically) on which it runs the docker container.
Binder does some additional work in the background -- if no software configuration files are provided in the GitHub repo, or if those contain a minimal set of software, binder will by default include JupyterHub in the docker.
When the cloud instance is launched, this is the screen you interact with.
You interact with the cloud instance in your browser.
Binders are ephemeral instances -- after a period of inactivity, the instance is automatically shut down, and any work you have done will be lost.
You're able to download files from your work before the instance is shut down if you do want to save anything.
</details>

## Directory structure

[00_data/](00_data/) Raman spectra and associated metadata.<br>
[01_notebooks/](01_notebooks/) Jupyter notebooks outlining the full suite of analyses presented in the pub.<br>
[02_outputs/](02_outputs/) Outputs associated with the pub.<br>

## R packages and versions used in this repo:

`pracma v2.4.2`<br>
`hyperSpec v0.100.0`<br>
`here v1.0.1`<br>
`ggplot v2_3.4.0`<br>
`lattice v0.20-45`<br>
`TreeDist v2.5.0`<br>
`ArcadiaColorBrewer v0.0.0.9000`<br>
`scales v0.0.0.9000`<br>
`phylosignal v1.3`<br>
`phytools v1.2-0`<br>
`maps v0.0.0.9000`<br>
`ape v3.4.1`<br>
`lsa v0.73.3`<br>
`SnowballC v0.7.0`<br>
`plotrix v3.8-2`<br>
`gplots v3.1.3`<br>
`reticulate v1.27`<br>

![github_front_page_figure](fig_1.png)
