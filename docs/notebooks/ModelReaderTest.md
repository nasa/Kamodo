SWMF model data reader for magnetosphere outputs.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
                Note: This reader 'walks' the entire dataset in the directory.
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                    (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
            filetime = boolean (default = False)
                - If False, the script fully executes.
                - If True, the script only executes far enough to determine the
                    time values associated with the chosen data.
            printfiles = boolean (default = False)
                - If False, the filenames associated with the data retrieved
                    ARE NOT printed.
                - If True, the filenames associated with the data retrieved ARE
                    printed.
            gridded_int = boolean (default = True)
                - If True, the variables chosen are functionalized in both the
                    standard method and a gridded method.
                - If False, the variables chosen are functionalized in only the
                    standard method.
            verbose = boolean (False)
                - If False, script execution and the underlying Kamodo
                    execution is quiet except for specified messages.
                - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.
            The run constants c, g, and R are saved in the kamodo.constants
            dictionary. R is the radius of the minimum inner boundary in R_E,
            c is the reduced speed of light, and g is the ratio of specific
            heats. Consult model documentation for more information.

        Notes and instructions:
        - The SWMF global magnetosphere outputs are given in one or more files
          per timestep in *.out files. The structure of these files are
          advantageous for efficient interpolation, so no file conversion is
          attempted.
        - The files are small and contain one time step per file, so
          interpolation method 1 is chosen. A custom interpolator is required.
        - The reader requires code in readers/OCTREE_BLOCK_GRID/:
          Compile the library using the command
            python interpolate_amrdata_extension_build.py
          inside the anaconda/miniconda3/miniforge3 environment for Kamodo.
          The shared library file name contains the python version and the name
          of the operation system, e.g.,
            _interpolate_amrdata.cpython-39-darwin.so
        - There is a pull request pending to implement sort_unstructured_data
          in the SpacePy software package. Until that PR is merged into
          SpacePy, you will need to manually modify the package. The
          lib/site-packages/spacepy/pybats/__init__.py script needs to be
          modified in your installation to not sort unstructured data in
          _read_Idl_bin() by default. See detailed instructions on the 2022 AGU
          poster https://doi.org/10.22541/essoar.167214301.16153548/v1 , at the
          bottom of the second to last column.