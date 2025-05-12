# Kamodo_FileDriver


The following scripts take in input from specific magnetospheric models and generates input files for specific ionosphere/thermosphere models. Kamodo is required for the code to run (https://github.com/nasa/Kamodo).

The current supported magnetosphere models are: SWMF_IE.

The current supported ionosphere/thermosphere models are: GITM, WACCMX, TIE-GCM

To run the file use the following command:

fileforcing_oneway(Model_A, Model_B, 'Directory where the input files for Model_A are', 'Directory to output files')

Ex:

fileforcing_oneway('SWMF_IE', 'GITM', "C:/Users/Directory_where_SWMF_IE_files_are", "C:/Users/Directory_where_to_output_files")

Notes: The masterscript.py is required for the individual filedriver scripts to work.
