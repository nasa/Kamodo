## Kamodo Installation Instructions
Kamodo is built to run with at least 16 GB of RAM. Attempting to run Kamodo with less memory may result in errors.  

In your Python environment: 
1. If you wish to create a new Python environment, use this command (replace Kamodo_env as desired):  
Note that Anaconda, MiniConda, etc. are not free for all to use, unfortunately. We recommend using [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) instead. Just replace 'conda' with 'micromamba' in the commands below.
> conda create -n Kamodo_env python=3.10  
> conda activate Kamodo_env  
2. Install Kamodo from pip (without SWMF-GM): Note this is currently out of date, use step 3.
> python -m pip install kamodo-ccmc  
3. Or you can download the latest Kamodo to the current directory and build: 
> git clone https://github.com/nasa/Kamodo.git  
> python -m pip install ./Kamodo
4. To build the SWMF-GM reader from the git clone (currently requires an editable pip install):  
> cd ./Kamodo/kamodo_ccmc/readers/OCTREE_BLOCK_GRID  
> python interpolate_amrdata_extension_build.py  
> cd ../../../..  
> python -m pip install -e ./Kamodo
5. To work with Kamodo you may also need iPython and/or Jupyter notebooks.  
> python -m pip install ipython jupyter  
-OR-  
> python -m pip install -r ./Kamodo/requirements.txt  

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
