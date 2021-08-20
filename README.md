# metagenomics_pooling_notebook

A Jupyter notebook to assist wet lab shotgun pipeline.

## Installation

To install this notebook, first clone the repository from GitHub:

```bash
git clone https://github.com/tanaes/metagenomics_pooling_notebook.git
```

Create a Python3 Conda environment in which to run the notebook:

```bash
conda create -n pooling_nb 'python>=3.6'
```

Activate the Conda environment:

```bash
source activate pooling_nb
```

Change directory to the downloaded repository folder and install:

```bash
cd metagenomics_pooling_notebook
pip install -e .
```

Finally, to enable the notebook Table of Contents, activate the
required nbextensions:

```bash
jupyter contrib nbextension install --sys-prefix --skip-running-check
jupyter nbextension enable toc2/main
```


## Use

The notebook itself contains detailed usage instructions. 

Ideally, you will copy this original notebook and use it for each plate you
create with the shotgun pipeline.

Just make sure to activate the right Conda environment before starting the
notebook.

Here's a quick example:

```bash
cp metagenomics_pooling.ipynb ~/New_project/new_project_pooling.ipynb
cd ~/New_project
source activate pooling_nb
jupyter notebook
```
