# How to Write a Model Reader for Diffusion Models

This guide demonstrates how the `verb3d` model reader is structured. 
Using this example, you can build and maintain your own diffusion model reader efficiently. 
Diffusion models, such as the VERB code, simulate the dynamics of particles (e.g., electrons) in multi-dimensional space, involving coordinate systems like L-shell, adiabatic invariants, and particle energies en pitch-angles.
Note, those coordinate systems are not natively supported by Kamodo, limiting the ability to interact with the diffusion model within Kamodo environment.

Unlike simpler models, diffusion models require modular handling of code, including careful management of multi-dimensional interpolation and integration with Kamodo’s function-based interface.
This documentation provides a practical approach for structuring your code to write your diffusion model reader.

In Kamodo, the model reader converts raw model output into files organized by model timesteps. 
This allows for efficient data access and processing by splitting the data into manageable chunks, improving performance and memory usage. 
The data or variable is then functionalized and exposed as an **interpolator** — a function that allows continuous data access across the model’s grid points for analysis and visualization.

## Key Components of a Diffusion Model Reader
The following are the essential components of a diffusion model reader, as illustrated by the `verb3d` model reader example:

### 1. **Variable Definitions and Metadata**
Like in other Kamodo model readers, in diffusion model readers, all variables must be defined in the `model_varnames` dictionary. This dictionary contains the metadata for each variable, including its LaTeX representation, description, coordinate system, and units. The `model_varnames` dictionary is placed at the top of the model reader module, following a consistent format across all model readers.

The structure of the `model_varnames` dictionary is as follows:

 ```python
 # model_varnames = {
 #     'Name in file': ['LaTeX representation', 'Description', integer, 
 #                      'Coordinate system', 'Coordinate grid', 
 #                      [Coordinate list], 'Units']
 # }
 ```
 
   - **Name in file**: The variable's name as it appears in the converted NetCDF files.
   - **LaTeX representation**: A LaTeX-formatted representation of the variable.
   - **Description**: A brief description of the variable.
   - **Integer**: An identifier for the variable (index number).
   - **Coordinate system**: The coordinate system associated with the variable (e.g., `LEA`).
   - **Coordinate grid**: Specifies the type of grid used. 'rb' represents the coordinate grid used in radiation belt models.
   - **Coordinate list**: The dimensions of the variable (e.g., `['time', 'L', 'E_e', 'alpha_e']`).
   - **Units**: The units of measurement.


**Example**:

```python
model_varnames = {
    'PSD': ['PSD_lea', 'Phase Space Density in (L, E_e, alpha_e)', 0, 'LEA', 'rb', ['time', 'L', 'E_e', 'alpha_e'], '(c/MeV/cm)**3'],
    'Flux': ['flux_lea', 'Electron flux (L, E_e, alpha_e)', 1, 'LEA', 'rb', ['time', 'L', 'E_e', 'alpha_e'], '1/(s*cm**2*keV*sr)'],
}
```
### 2. **Data Conversion**
Model reader needs to convert raw model output into a Kamodo-compatible format, such as NetCDF.

**Modular Processing**: To ensure modularity, data conversion is handled by a separate module. For instance, the `verb3d_tocdf.py` provides functionality to process the raw output, convert it to NetCDF, and store it for use by the main model reader.

### 3. **Handling Coordinate Systems**
**Coordinate System Handling**: Diffusion models use coordinate systems that are different from the space-based systems (e.g., GEO, GSM). Instead of spatial coordinates, these systems are often related to physical properties like adiabatic invariants or particle motion. For instance, the `verb3d` reader operates in two specific coordinate systems: 
 - **LEA (L-shell, Electron Energy, Equatorial Pitch Angle)**: This system is based on the L-shell `L`, electron energy (`E_e`), and equatorial pitch angle (`alpha_e`). These are common parameters used in radiation belt models.
 - **LMK (L-shell, Mu, K)**: In this system, `L` represents the L-shell, `Mu` and `K` are the first and second adiabatic invariants. This system is commonly used to analyze electron phase space density.

**Grid Definitions**: Kamodo does not automatically convert between these complex, physics-based coordinate systems and spatial systems like GEO or GSM. Therefore, to prevent unwanted coordinate conversion, the type `'rb'` (radiation belts) is used to indicate that these grids should remain in their specialized coordinate systems.

### 4. **Modularity**
Diffusion models often require complex data processing, which can quickly clutter the code based on the other model reader examples. Therefore, writing the diffusion model reader requires a modular approach, where key tasks are separated into individual functions within the same model python module.

## Refactored Functions in the VERB Code Model Reader
This section shows how the `verb3d_4D.py` code has been restructured from simpler models like `dtm` by explaining the **refactored functions** in the diffusion model reader.

### 1. **Initialization**
- **`__init__(...)`**: Constructor is rewritten to be modular, but it still initializes the model reader by setting up the file directory, loading the requested variables, and preparing the model for Kamodo funtionalization.

### 2. **Handling Time and File Lists**
Check for the presence of the list and time files in the of the model output data. If they are missing, the model reader creates them and converts model's raw output to NetCDF files. 

- **`time_list_files(file_dir)`**:  Check if '_list.txt' and '_time.txt' files exist. Creates dataset files if needed.
- **Parameters**:
  - **`file_dir`**: Directory of the model output data.

### 3. **Initial Variable Check and Validation**
Performs an initial check on the requested variables to ensure they are valid and processes requests for all variables.
- **`_variable_requested_check(variables_requested)`**: Checks if the requested variables are valid by comparing them to the available variables in the model.
- **Parameters**:
  - **`variables_requested`**: A list of requested variables.

### 4. **Storing and Organizing Variables**
The model organizes the variables by mapping them to their respective files and stores this information for later use.

**User refactored functionality**:
- **`_model_vars.file_variables(...)`**: Retrieves the list of file variables that match the requested variables from the dataset.
- **`varfiles[p] = [_model_vars.vars[v].var for v in _gvar_list]`**: Maps the LaTeX names of the requested variables to their corresponding file pattern (`p`). 
- **`self.gvarfiles_full.update(...)`**: Updates the dictionary with all available file variables for the current file pattern, storing them for future use.

### 5. **Error Handling and loging**
Handle preliminary errors and logs information about the variables and files processed, helping to identify issues.
- **`print_err_list(variables_requested, varfiles)`**: Prints a list of errors related to the requested variables, helping identify which variables failed to load or register properly.
- **`_print_files(printfiles)`**: Prints a list of files that have been processed, aiding in debugging by showing which files were successfully read and loaded into the model.
- **Parameters**:
  - **`variables_requested`**: The list of requested variables.
  - **`varfiles`**: The list of files associated with the requested variables.
  - **`printfiles`**: The list of files to be printed for debugging.

### **6. Variable Definitions and Metadata**
Create and update the internal variable dictionary to ensure all relevant data is prepared for variables to be functionalized.

- **`_update_var_dict(varfiles)`**: Updates the internal variable dictionary according to the model output files. Refactored to be called when `variables_requested == 'all'`.
- **`_create_variables()`**: Initialized storage structure and preparing them for use in Kamodo.
- **Parameters**:
  - **`varfiles`**: The list of files associated with the requested variables.

### 7. **Loading Grids and Static Variables**
Load the necessary grid data and static variables from files and store them as class properties.
- **`load_grid(file_dir)`**: Loads grid data and making it available as class properties.
- **`load_static(file_dir)`**: Loads static data (e.g., `pc`) and making it available as class properties.
- **Parameters**:
  - **`file_dir`**: Directory of the model output data.

### 8. **Registering Variables**
Register the requested variables in Kamodo. 
This includes handling different types of variables, such as PSD, 1D, and grid variables, each requiring different registration methods.
- **`register_ncfile_variable(varname, gridded_int, file_dir)`**: Registers a variable from an NetCDF by extracting its metadata and grid information and prepares it for Kamodo.
- **`register_grid_variable(varname, gridded_int, file_dir)`**: Similar to `register_ncfile_variable`, but specifically handles variables that are part of a grid.
- **`register_1d_variable(varname, gridded_int, file_dir)`**: Registers 1D variables, which are simpler and do not involve complex grids, and prepares them for use in Kamodo.
- **Parameters**:
  - **`file_dir`**: Directory of the model output data.
  - **`varname`**: The name of the registered variable.
  - **`gridded_int`**: Kamodo flag, indicates whether the variable should be gridded.

### 9. **Final Timing and Verbose Output**
Log the time taken to register the variables, providing performance metrics if verbosity is enabled.

### **Handling interpolations**
This section explains interpolation functions for the data within different coordinate systems (LEA and LMK).
The registration of variables depends on the type of variable being functionalized (e.g., PSD variables, 1D variables, or grid variables). 
Based on the variable's characteristics, different `register_*` functions are called, which handle the interpolation process accordingly.

For **PSD and flux variables** such as `PSD_lmk`, `PSD_lea`, and `flux_lea`, **multi-dimensional interpolation** is used.
This method interpolates smoothly across complex grids in the LEA (`L`, `E_e`, `alpha_e`) or LMK (`L`, `mu`, `K`) coordinate systems.

For **1D variables** like `bf`, `kp`, and `Lpp`, as well as **grid variables**, **nearest-neighbor interpolation** is used.
This method identifies the closest data point, which is usually sufficient for these types of variables.

- **`_nearest_xvec_data_LEA(xvec, data)`**: Finds the nearest data point to the given input vector `xvec` in the LEA coordinate system, using nearest-neighbor interpolation to retrieve the data.
- **`_nearest_xvec_data_LMK(xvec, data)`**: Similar to `_nearest_xvec_data_LEA`, but operates in the LMK coordinate system.
- **`_interp_xvec_data_LEA(xvec, data)`**: Interpolates data based on the input vector `xvec` in the LEA coordinate system, preparing the data for use in Kamodo functions.
- **`_interp_xvec_data_LMK(xvec, data)`**: Interpolates data based on `xvec` in the LMK coordinate system, ensuring smooth data handling for Kamodo.
- **`_nearest_index(arr, val)`**: Finds the nearest index in the array `arr` for the given value `val`, used in various interpolation and nearest-neighbor functions to efficiently locate data points.
- **`_interpolate(X, Y, Z, x, y, z, Psi_tt)`**: Handles multi-dimensional interpolation of data points in 3D space, ensuring accurate data retrieval across complex grids and coordinate systems.
- **`_intrep_stack_hstack(self, Xit, Xi, Yit, Yi, Zit, Zi, PSDit, PSDi)`**: This function handles stacking grid points and interpolation values using `np.hstack`.
- **`_intrep_stack_list(self, Xit_list, Xi, Yit_list, Yi, Zit_list, Zi, PSDit_list, PSDi)`**: Similar to `_intrep_stack_hstack`, but stores grid points and interpolation values in lists for further processing.
- **`_finalize_stack(Xit_list, Yit_list, Zit_list, PSDit_list)`**: Converts lists of stacked grid points and values into arrays after all data has been processed.

- **Parameters**:
  - **`xvec`**: The input vector used for interpolation calculations.
  - **`data`**: The dataset to be interpolated.
  - **`arr`**: The array in which the nearest index is to be found.
  - **`val`**: The value for which the nearest index is sought.
  - **`X`, `Y`, `Z`**: The multi-dimensional coordinates used for interpolation.
  - **`x`, `y`, `z`**: The specific input coordinates for interpolation.
  - **`Psi_tt`**: The data array being interpolated.
  - **`Xit`, `Yit`, `Zit`**: Arrays or lists of grid coordinates, typically representing different dimensions in the interpolation (e.g., `L`, `E_e`, `alpha_e` or `L`, `mu`, `K`).
  - **`Xi`, `Yi`, `Zi`**: Individual sets of grid coordinates, corresponding to specific points in the data grid.
  - **`PSDit`**: The array or list of phase space density (or other variable) values that correspond to the grid points.
  - **`PSDi`**: Individual phase space density (or other variable) values at specific grid points.
  - **`Xit_list`, `Yit_list`, `Zit_list`, `PSDit_list`**: Lists that accumulate grid coordinates and corresponding values for further stacking and processing.

### **Model Variable Management**
In the VERB code model reader, the `ModelVariable` and `ModelVariables` classes are introduced to simplify handling and managing `model_varnames` associated variables. 

#### `ModelVariable` Class
A data class representing a single variable within the diffusion model, encapsulating its metadata and providing access to essential attributes. This class simplifies the management of Kamodo variables by replacing the previous list-based structure with an object-oriented approach. 
**Methods**:
- **`to_list()`**: Converts the `ModelVariable` instance into a list format, representing `model_varnames`.

#### `ModelVariables`  Class
A container class that manages multiple `ModelVariable` instances. It provides methods for accessing, filtering, and retrieving variable metadata by name or other criteria.

**Attributes**:
- **`vars`** (`dict`): A dictionary storing all model variables, where the key is the variable's LaTeX representation name and the value is an instance of `ModelVariable`.
- **`keys`** (`dict`): A dictionary that maps the LaTeX representation to the variable's internal name.

**Methods**:
- **`__init__(model_varnames)`**: Initializes the `ModelVariables` class with a dictionary of variable metadata (`model_varnames`).
  - **Parameters**: 
    - `model_varnames` (`dict`): A `model_varnames` dictionary corresponding to the model. 
- **`latex_variables(variables_requested=None)`**: Returns a list of LaTeX representations for all variables or filtered by a requested variable name.
  - **Parameters**:
    - `variables_requested` (`str`, optional): The name of a specific variable to filter by.
  - **Returns**: (`list`): A list of LaTeX-formatted variable names.

- **`file_variables(variables_requested=None, cdf_keys=None)`**: Returns a list "Name in file" representations that match the requested variable name and/or "Name in file" keys.
  - **Parameters**:
    - `variables_requested` (`str`, optional): The name of the variable to filter by.
    - `cdf_keys` (`list`, optional): A list of "Name in file" keys to match with the variables.
  - **Returns**: (`list`): A list of "Name in file" variables.

- **`model_var_by_name(name)`**: Retrieves a specific `ModelVariable` instance by its internal name.
  - **Parameters**:
    - `name` (`str`): The internal name of the model variable.
  - **Returns**: (`ModelVariable`): The corresponding `ModelVariable` instance.

#### Usage
In the `verb3d_4D.py` model reader, `_model_vars` is initialized and used to efficiently manage the model's variables by providing structured access to the variables' metadata. It enables functions to retrieve specific variables, access their coordinate systems, units, and descriptions, and ensures proper alignment of data arrays. For example, `_model_vars` is utilized to access coordinate grids, units for variables like `L` and `E_e`, and to map variables from the file to their respective LaTeX representations.
```python
_model_vars = ModelVariables(model_varnames)
```

### **Static Control Properties in the `MODEL(Kamodo)` Class**

The `MODEL(Kamodo)` class in the VERB reader includes several static properties that control default behaviors and configurations for processing and functionalizing data.

**Note**: Many of these properties are primarily intended for development and debugging purposes, and should not be modified. Some properties like `_interpolation_method` will not affect the model reader's behavior unless explicitly modified within the class's code.
 
- **`modelname`**: The name of the model being processed. In this case, it is set to `'VERB-3D'`, identifying the specific diffusion model.
- **`_grid_prefix`**: This prefix, set to `'_grid'`, is used to identify and store grid variables within the class. It helps differentiate grid data from other types of variables during processing.
- **`_force_convert_all`**: A boolean flag that, when set to `True`, forces the reader to recompute all NetCDF files, even if they already exist. By default, it is set to `False`, indicating that existing files will be used unless otherwise specified.
- **`_start_date`**: This property allows the user to specify a forced start date for the dataset. By default, it is set to `None`, meaning the start date is determined from the data itself.
- **`_interpolation_method`**: This property defines the default interpolation method used for functionalizing PSD and flux variables. It is set to `'interp'` for multi-dimensional interpolation. The commented-out option `'nearest'` indicates that nearest-neighbor interpolation could be used instead if required.


## Conclusion

Writing a model reader for diffusion models requires careful consideration of data conversion, coordinate systems, and modularity. By following the example of the `verb3d` model reader, you can create a robust and maintainable reader that integrates with Kamodo.
