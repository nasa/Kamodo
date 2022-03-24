## Kamodo CCMC Installation Instructions

In conda command prompt:
1. Move to the directory where you want the Kamodo package to be stored or if you wish to create a new environment, use this command:
> conda create -n Kamodo_env python=3.7  
2. Add the packages needed by the CCMC readers to the desired environment (replace 'Kamodo_env' with your environment name):
> conda install -n Kamodo_env -c conda-forge netCDF4 xarray dask astropy ipython jupyter
3. Activate the desired environment. 
> conda activate Kamodo_env
4. Install remaining dependencies:
> python -m pip install --upgrade spacepy  
> python -m pip install hapiclient
5. Download Kamodo to the current directory:
> git clone https://github.com/nasa/Kamodo.git
6. Install the Kamodo package. (Check the directory structure before using this command. The ./Kamodo directory should contain the kamodo_ccmc directory.)
> python -m pip install ./Kamodo     

### Testing commands from ipython or notebook session:
```
from kamodo import Kamodo
k = Kamodo()  
import kamodo_ccmc.flythrough.model_wrapper as MW  
MW.Model_Variables('OpenGGCM_GM')
```

Correct output:
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

### Video Tutorial Links  
Kamodo Flythrough Tutorial: https://www.youtube.com/watch?v=1I2BZBl-wl4  
Kamodo Onboarding Tutorial: https://www.youtube.com/watch?v=nvl61pklEuU  
