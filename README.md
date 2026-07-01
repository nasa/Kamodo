# Kamodo

![Python 3.6](https://img.shields.io/badge/Python-3.6-blue.svg)
![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg)
![Python 3.8](https://img.shields.io/badge/Python-3.8-blue.svg)
![Python 3.9](https://img.shields.io/badge/Python-3.9-blue.svg)
![Python 3.10](https://img.shields.io/badge/Python-3.10-blue.svg)
![Python 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)



[![codecov](https://codecov.io/gh/asherp/Kamodo/branch/master/graph/badge.svg?token=W1B3L19REF)](https://codecov.io/gh/asherp/Kamodo)

Kamodo is an open source CCMC tool for access, interpolation, and visualization of space weather models and data in python.
Kamodo allows model developers to represent simulation results as mathematical functions which may be manipulated directly by end users.
This general approach allows observational data to be represented functionally, through the use of interpolators.
Kamodo handles unit conversion transparently and supports interactive science discovery in a low coding environment through jupyter notebooks.
These features allow Kamodo to be used in other fields of study and as a teaching tool for working with real world physical data.


This repository hosts the core Kamodo libraries under a permissive [NASA open source license](https://github.com/nasa/Kamodo-core//blob/master/LICENSE).
The core library supports function registration, composition, unit conversion, automated plotting, LaTeX I/O, and remote procedure call (RPC) interfaces.

Space weather simulation readers are implemented as subclasses of the Kamodo base class and are developed and maintained by the Community Coordinated Modeling Center, located at NASA Goddard Space Flight Center. CCMC's Kamodo readers may be found here [https://github.com/nasa/Kamodo/](https://github.com/nasa/Kamodo/)


## Usage
Suppose we have a vector field defined by a function of positions in the x-y plane:

```python
from kamodo import kamodofy
import numpy as np

x = np.linspace(-np.pi, np.pi, 25)
y = np.linspace(-np.pi, np.pi, 30)
xx, yy = np.meshgrid(x,y)
points = np.array(list(zip(xx.ravel(), yy.ravel())))

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
![usage](https://raw.github.com/nasa/Kamodo-core/main/docs/notebooks/images/fig2d-usage.svg)

Head over to the [Introduction](notebooks/Kamodo.ipynb) page for more details.


## Getting started

Kamodo may be installed from pip

```console
pip install kamodo
```

To get the latest version of Kamodo Core, install from the NASA git repo:

```console
pip install git+https://github.com/nasa/Kamodo-core.git
```

### Kamodo Environment

We strongly recommend using the conda environment system to avoid library conflicts with your host machine's python.

Download and install miniconda from [here](https://conda.io/miniconda.html). The advantage to using miniconda is that each new environment includes the bare-minimum for a project. This allows you to keep many different projects on a single work station.

#### Create Kamodo environment

Create a new environment for kamodo

```console
conda create -n kamodo python=3.10
conda activate kamodo
(kamodo) pip install kamodo
```
!!! note
    The leading (kamodo) in your prompt indicates that you have activated the `kamodo` environment.
    From here on, anything you install will be isolated to the `kamodo` environment.

#### Loading example notebooks

If you want to run any of the notebooks in docs, you will need to install `jupyter`:

```console
(kamodo) conda install jupyter
```

Navigate to the top-level of the kamodo repo, then point jupyter to `docs/notebooks`:

    (kamodo) jupyter notebook docs/notebooks

This should open a browser window that will allow you to load any of the example notebooks.

#### Requirements

The following (minimum) requirements are obtained by running `pip install kamodo`

* decorator>=4.4.2
* numpy
* scipy
* sympy==1.5.1
* pandas
* plotly
* pytest
* hydra-core==0.11.3
* Flask==1.1.2
* flask-cors
* flask-restful==0.3.8
* antlr4-python3-runtime==4.7
* python-forge
* requests
* incremental
* pycapnp
* pyOpenSSL


The antlr package may be necessary for rendering latex in a notebook

```sh
conda install antlr-python-runtime
```

Plotly-orca may be needed for proper image export

```sh
conda install -c plotly plotly-orca (for writing images)
```


## Test Suite

Kamodo's unit tests are run with [pytest](https://docs.pytest.org/en/7.0.x/). To install pytest with code coverage

```sh
python -m pip install flake8 pytest
pip install pytest-cov
```

Then, from the base of the git repo:

```sh
pytest --cov kamodo.kamodo --cov kamodo.util --cov plotting kamodo/test_plotting.py kamodo/test_kamodo.py kamodo/test_utils.py
```

This will generate a test report and coverage of the `kamodo` module.

To run RPC tests, you must first generate a self-signed certificate.

```sh
python kamodo/rpc/gen_self_signed_cert.py certfile
# certfile.key and certfile.cert will be placed in your local directory
pytest kamodo/rpc/test_rpc_threaded.py
```

