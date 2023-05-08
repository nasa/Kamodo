## Kamodo Installation Instructions
Kamodo is built to run with at least 16 GB of RAM. Attempting to run Kamodo with less memory may result in errors.  

In your Python environment: 
1. If you wish to create a new conda environment, use this command (replace Kamodo_env with your own name): 
> conda create -n Kamodo_env python=3.10  
> conda activate Kamodo_env  
2. Install Kamodo from pip:
> python -m pip install kamodo-ccmc  
3. Or you can download Kamodo to the current directory: 
> git clone https://github.com/nasa/Kamodo.git  
> python -m pip install ./Kamodo  
4. To work with Kamodo you may also need iPython and/or Jupyter notebooks.  
> python -m pip install ipython jupyter  

NOTE: Sometimes an error will occur installing the spacepy dependency if numpy is not yet installed. 
Running 'python -m pip install numpy' then repeating the Kamodo pip install usually resolves it. 
If that does not resolve the issue, check out the spacepy troubleshooting page: 
https://spacepy.github.io/install.html#troubleshooting

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
