![CCMC](docs/notebooks/Files/CCMC.png) ![Kamodo](docs/notebooks/Files/Kamodo.png)
# The CCMC Kamodo Analysis Suite
## Vision Statement
Kamodo is an official NASA open-source python package built upon the functionalization of datasets. Once a dataset is functionalized in Kamodo, several important capabilities are then available to the user, including data analysis via function composition, automatic unit conversions, and publication quality graphics all using intuitive and simplistic syntax. By applying these capabilities to heliophysics model outputs, we aim to:
-	Drastically simplify the currently complex data utilization process for model outputs,
-	Provide interactive access to functionalized model outputs for users ranging in programming skill from beginners – via code-free interfaces and video tutorials – to advanced users – via thorough documentation, Jupyter notebook examples and sample workflows,
-	Layer multiple functionalities on top of the functionalized model outputs, all with model-agnostic and uniform syntax, including but not limited to:
    - Flythrough tools,
    - Vector field tracing (including magnetic field mapping),
    - Coordinate conversions,
    - Domain-specific interactive plots of publication quality,
    - Modular driver swapping,
    - Satellite constellation mission planning tools,
    - Simulated imagery, and
    - A line of sight calculation tool,
-	Greatly reduce the programming skill currently required outside of Kamodo to perform model validation studies and model-data comparisons,
-	Enable model output utilization both on the cloud and on personal laptops in a variety of methods (e.g. through HAPI and interactive calls from the command line),
-	Streamline the CCMC user workflow by becoming interoperable with other CCMC services (e.g. CAMEL and the various scoreboards),
-	And become the next generation interface for CCMC users to interact with and analyze model outputs (e.g. through ROR and IR),

...all while keeping the developed software open-source and freely available. The Kamodo team also supports the heliophysics community by pursuing interoperability with commonly-used python packages, collaborating with community members to add model outputs and new functionalities, and remaining involved with community events (e.g. conferences, challenges, and research support). As the library of supported model outputs types expands and new model-agnostic tools are added, Kamodo will become a staple software package in the heliophysics community to transform current workflows into a more efficient and productive process. We are building the next generation of capability with Kamodo. Join us!

## Kamodo Installation Instructions   

### Conda prompt commands: 
- Move to the directory where you want the Kamodo package to be stored or if you wish to create a new environment, use this command:

> conda create -n Kamodo_env python=3.7  

- Add the packages needed by the CCMC readers to the desired environment (replace 'Kamodo_env' with your environment name):

> conda install -n Kamodo_env -c conda-forge netCDF4 cdflib astropy ipython jupyter h5py sgp4

- Activate the desired environment. 

> conda activate Kamodo_env

- Install remaining dependencies:

> python -m pip install --upgrade spacepy  
> python -m pip install hapiclient    

- Download CCMC Kamodo to the current directory:

> git clone https://github.com/nasa/Kamodo.git

- Install the CCMC Kamodo package. (Check the directory structure before using this command. The ./Kamodo directory should contain the kamodo_ccmc directory.)

> python -m pip install ./Kamodo 

Note: Developers should install CCMC Kamodo with the -e option

### Testing commands from an ipython or notebook session


```python
from kamodo import Kamodo
k = Kamodo()  
import kamodo_ccmc.flythrough.model_wrapper as MW  
MW.Model_Variables('OpenGGCM_GM')
```

    
    The OpenGGCM_GM model accepts the standardized variable names listed below.
    -----------------------------------------------------------------------------------
    B_x : '['x component of magnetic field', 0, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
    B_y : '['y component of magnetic field', 1, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
    B_z : '['z component of magnetic field', 2, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'nT']'
    E_x : '['x component of electric field', 6, 'GSE', 'car', ['time', 'x', 'x', 'x'], 'mV/m']'
    E_y : '['y component of electric field', 7, 'GSE', 'car', ['time', 'y', 'y', 'y'], 'mV/m']'
    E_z : '['z component of electric field', 8, 'GSE', 'car', ['time', 'z', 'z', 'z'], 'mV/m']'
    N_plasma : '['number density of plasma (hydrogen equivalent)', 12, 'GSE', 'car', ['time', 'x', 'y', 'z'], '1/cm**3']'
    P_plasma : '['plasma pressure', 14, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'pPa']'
    eta : '['resistivity', 13, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'm**2/s']'
    j_x : '['current density, x component', 15, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
    j_y : '['current density, y component', 16, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
    j_z : '['current density, z component', 17, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'muA/m**2']'
    v_plasmax : '['x component of plasma velocity', 9, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
    v_plasmay : '['y component of plasma velocity', 10, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
    v_plasmaz : '['z component of plasma velocity', 11, 'GSE', 'car', ['time', 'x', 'y', 'z'], 'km/s']'
    
    

## Citing Kamodo

When publishing research which used Kamodo, please provide appropriate credit to the CCMC and the Kamodo team via citation or acknowledgment. Please also let the team know of publications or presentations that use Kamodo. Below is list of publications for Kamodo.

- Pembroke, A., D. De Zeeuw, L. Rastaetter, R. Ringuette, O. Gerland, D. Patel and M. Contreras (2022). Kamodo: A functional API for space weather models and data. JOSS 7, 75, 4053, https://doi.org/10.21105/joss.04053.

- Ringuette, R., D. De Zeeuw, L. Rastaetter, A. Pembroke, O. Gerland, K. Garcia-Sage (2022). Kamodo’s model-agnostic satellite flythrough: Lowering the utilization barrier for heliophysics model outputs, Frontiers in Astronomy and Space Sciences, vol 9. http://dx.doi.org/10.3389/fspas.2022.1005977.

- Ringuette, R., D. De Zeeuw, L. Rastaetter, A. Pembroke, O. Gerland, K. Garcia-Sage (2022). Kamodo’s model-agnostic satellite flythrough: Lowering the utilization barrier for heliophysics model outputs, Frontiers in Astronomy and Space Sciences, vol 9. http://dx.doi.org/10.3389/fspas.2022.1005977.

- Ringuette, R., L. Rastaetter, D. De Zeeuw, K. Garcia-Sage, R. Robinson, and O. Gerland (2022). Kamodo's Satellite Constellation Mission Planning Tool, poster presentation presented by L. Rastaetter at the 2022 Fall meeting of AGU, Dec 12-16, Chicago, IL, USA. https://doi.org/10.22541/essoar.167214257.73153757/v1.

- Ringuette, R., L. Rastaetter, D. De Zeeuw, A. Pembroke, and O. Gerland (2023). Simplifying model data access and utilization. Adv. Space. Res. under review.Use of Kamodo in published work or in other development projects requires citing one or more of the publications below as indicated.


## Resources
- CCMC's Kamodo Official website - https://ccmc.gsfc.nasa.gov/tools/kamodo/  
- CCMC's Kamodo Documentation page - https://nasa.github.io/Kamodo/  
- Sample model outputs - https://ccmc.gsfc.nasa.gov/RoR_WWW/output_files/KAMODO_DEMO/  
- Youtube tutorial channel - https://www.youtube.com/playlist?list=PLBWJQ5-pik_yBBcrpDRPM2hLluh-jreFa  

## The Kamodo team
**Dr. Rebecca Ringuette**  
- ORCiD: https://orcid.org/0000-0003-0875-2023  
- NASA Staff Page: https://ccmc.gsfc.nasa.gov/staff/rebecca-ringuette/

**Dr. Lutz Rastaetter**  
- ORCiD: https://orcid.org/0000-0002-7343-4147  
- NASA Staff Page: https://ccmc.gsfc.nasa.gov/staff/lutz-rastaetter/

**Dr. Darren De Zeeuw**  
- ORCiD: https://orcid.org/0000-0002-4313-5998  
- NASA Staff Page: https://ccmc.gsfc.nasa.gov/staff/darren-de-zeeuw/

**Dr. Katherine Garcia-Sage**  
- ORCiD: https://orcid.org/0000-0001-6398-8755  
- NASA Staff Page: https://ccmc.gsfc.nasa.gov/staff/katherine-garcia-sage/


## Open-Source License
Kamodo is an official NASA open source software package. Kamodo's official source code is hosted on github under a permissive NASA open source license: For more details, go here: https://github.com/nasa/Kamodo/blob/master/LICENSE
