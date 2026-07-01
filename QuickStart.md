## Kamodo Installation Instructions

Kamodo is built to run with at least 16 GB of RAM. Attempting to run Kamodo with less memory may result in errors.

### Quick Install (Recommended)

**Step 1:** Create a Python environment (optional but recommended)

Note that Anaconda, MiniConda, etc. are not free for all to use. We recommend using [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) instead. Just replace 'conda' with 'micromamba' in the commands below.

```bash
conda create -n Kamodo_env python=3.10
conda activate Kamodo_env
```

**Step 2:** Install Kamodo

Choose one of the following options:

**Option A: Install from PyPI** (when available - coming soon!)
```bash
pip install kamodo-ccmc
```

**Option B: Install from Git Repository** (current recommended method)
```bash
git clone https://github.com/nasa/Kamodo.git
pip install ./Kamodo
```

**Option C: Developer Installation** (for contributing to Kamodo)
```bash
git clone https://github.com/nasa/Kamodo.git
pip install -e ./Kamodo
```

**That's it!** C and Fortran extensions will be automatically compiled during installation if you have the appropriate compilers installed (gcc for C, gfortran for Fortran).

**Step 3:** (Optional) Install Jupyter for interactive work
```bash
pip install ipython jupyter
```

### Compiler Requirements

**For full functionality**, you need:
- **gcc** (C compiler) - Required for SWMF-GM and GAMER-AM readers
- **gfortran** (Fortran compiler) - Required for OpenGGCM reader

**If you don't have compilers:** Don't worry! Kamodo will still install successfully. You'll see warnings about which readers are unavailable, but all other model readers will work normally.

**To install compilers:**
- **Linux:** `sudo apt-get install gcc gfortran`
- **macOS:** `xcode-select --install` (for gcc), `brew install gcc` (for gfortran)
- **Windows:** `conda install -c conda-forge m2w64-gcc-fortran libpython`  

### Testing commands:
```
from kamodo import Kamodo
k = Kamodo()  
import kamodo_ccmc.flythrough.model_wrapper as MW  
MW.Model_Variables('OpenGGCM_GM')
```

Output should be similar to this:
```
The model accepts the standardized variable names listed below.
-----------------------------------------------------------------------------------
B_x : '['x component of magnetic field', 0, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
B_y : '['y component of magnetic field', 1, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
B_z : '['z component of magnetic field', 2, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
B1_x : '['x component of magnetic field (on grid cell faces)', 3, 'GSE', 'car', ['time', 'x', 'x', 'x'], 'nT']'
B1_y : '['y component of magnetic field (on grid cell faces)', 4, 'GSE', 'car', ['time', 'y', 'y', 'y'], 'nT']'
B1_z : '['z component of magnetic field (on grid cell faces)', 5, 'GSE', 'car', ['time', 'z', 'z', 'z'], 'nT']'
E_x : '['x component of electric field (on grid cell edges)', 6, 'GSE', 'car', ['time', 'x', 'x', 'x'], 'mV/m']'
E_y : '['y component of electric field (on grid cell edges)', 7, 'GSE', 'car', ['time', 'y', 'y', 'y'], 'mV/m']'
E_z : '['z component of electric field (on grid cell edges)', 8, 'GSE', 'car', ['time', 'z', 'z', 'z'], 'mV/m']'
V_x : '['x component of plasma velocity', 9, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
V_y : '['y component of plasma velocity', 10, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
V_z : '['z component of plasma velocity', 11, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
N_plasma : '['plasma number denstity (hydrogen equivalent)', 12, 'GSE', 'car', ['time', 'x', 'y', 'z'], '1/cm**3']'
eta : '['resistivity', 13, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'm**2/s']'
P_plasma : '['plasma pressure', 14, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'pPa']'
J_x : '['x component of current density', 15, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
J_y : '['y component of current density', 16, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
J_z : '['z component of current density', 17, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
```

### Caveats
Kamodo does not generate model outputs. Users need to acquire the desired model outputs before they can be functionalized by Kamodo.

### Video Tutorial Channel  
https://www.youtube.com/playlist?list=PLBWJQ5-pik_yBBcrpDRPM2hLluh-jreFa
