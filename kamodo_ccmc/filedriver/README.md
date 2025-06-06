# Kamodo_FileDriver

The following scripts read supported magnetospheric model output and generates input files for supported  ionosphere/thermosphere models. Many features of Kamodo are used to support the reading, interpolation, and processing required. 

The current supported magnetosphere models are: 
- SWMF_IE. 

The current supported ionosphere/thermosphere models are: 
- GITM
- WACCMX
- TIE-GCM 

To run this code, use the following command: 

> import kamodo_ccmc.filedriver.master_script_v2 as fd  
> fd.fileforcing_oneway(Model_A, Model_B, 'Directory for Model_A files', 'Directory for new output files') 

For example:

> import kamodo_ccmc.filedriver.master_script_v2 as fd  
> fd.fileforcing_oneway('SWMF_IE', 'WACCMX', "DATA/SWMF-IE_Data", "OUTPUTfiles") 

NOTE: The file **KP_Daily.npz** needs to be in the working directory as it is needed to compute some terms. A sample file along with a script to update the data as needed are included here. Just copy the included file to your working directory.  
