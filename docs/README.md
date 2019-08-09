# Kamodo Project Page

Kamodo is a new CCMC tool for access, interpolation, and visualization of space weather models and data in python. Kamodo allows model developers to represent simulation results as mathematical functions which may be manipulated directly by end users. Kamodo handles unit conversion transparently and supports interactive science discovery through jupyter notebooks with minimal coding and is accessible through python.

## Usage
Suppose we have a vector field defined by a function of positions in the x-y plane:

```python
from kamodo import kamodofy
import numpy as np

x = np.linspace(-np.pi, np.pi, 25)
y = np.linspace(-np.pi, np.pi, 30)
xx, yy = np.meshgrid(x,y)
points = np.array(zip(xx.ravel(), yy.ravel()))

@kamodofy(units = 'km/s')
def fvec(rvec = points):
    ux = np.sin(rvec[:,0])
    uy = np.cos(rvec[:,1])
    return np.vstack((ux,uy)).T
```

The @kamodofy decorator lets us register this field with units to enable unit-conversion downstream:
```python
from kamodo import Kamodo

kamodo = Kamodo(fvec = fvec)
kamodo
```
When run in a jupyter notebook, the above kamodo object will render as a set of equations:

$$\vec{f}{\left (\vec{r} \right )} [km/s] = \lambda{\left (\vec{r} \right )}$$

We can now evaluate our function using dot notation:

```python
kamodo.fvec(np.array([[-1,1]]))
```
```console
array([[-0.84147098,  0.54030231]])
```
We can perform unit conversion by function composition:
```python
kamodo['gvec[m/s]'] = 'fvec'
```
kamodo automatically generates the appropriate multiplicative factors:
$$\vec{g}{\left (\vec{r} \right )} [m/s] = 1000 \vec{f}{\left (\vec{r} \right )}$$
we can verify these results through evaluation

```python
kamodo.gvec(np.array([[-1,1]]))
```
```console
array([[-841.47098481,  540.30230587]])
```
Kamodo also generates quick-look graphics via function inspection.
```python
import plotly.io as pio

fig = kamodo.plot('fvec')
pio.write_image(fig, 'images/fig2d-usage.svg')
```
![usage](notebooks/images/fig2d-usage.svg)

Head over to the [Introduction](notebooks/Kamodo.ipynb) page for more details.

## Getting started

Kamodo is now available through pip
```console
pip install kamodo
```

### Where to download

If you have a nasa.developer.gov account, you may access the kamodo repository with [git](https://git-scm.com/):

    git pull https://developer.nasa.gov/CCMC/Kameleon2.0

!!! note
    Kamodo is currently only available to users with NASA credentials. We are in the process of making Kamodo open to the public as an open-source project. 

#### Download (mini)conda

We strongly recommend using the conda environment system to avoid library conflicts with your host machine's python.

Download and install miniconda from [here](https://conda.io/miniconda.html).

#### Create Kamodo environment

Create a new environment for kamodo

    conda create -n kamodo-user jupyter

#### Activate new environment

From a bash shell:

    source activate kamodo

Requirements

* numpy
* scipy
* sympy
* pandas
* plotly==3.3 
* pytest
* psutil
* conda install antlr-python-runtime (rendering latex)
* conda install -c plotly plotly-orca (for writing images)

!!! note
    plotly version in flux

#### Loading a notebook

Start the notebook server in the ```Prototypes``` subdirectory

    jupyter notebook Prototypes

#### Documentation

To generate this documentation site, you will need to install a few extra packages:

    pip install mkdocs mknotebooks
    pip install python-markdown-math

To start the documentation server, navigate to the root of the repository, then:

    mkdocs serve
