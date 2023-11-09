"""
@author: xandrd
"""

# model_varnames = {'Name in file': ['LaTeX representation', 'Description', integer, 'Coordinate system',
#                            'Coordinate grid', [Coordinate list], 'Units'],

model_varnames = {'PSD': ['PSD_{lea}', 'Phase Space Density in (L, E, A)', 0, 'LEA',
                            'car', ['time', 'L', 'Energy', 'Alpha'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'PSD_2': ['PSD_{lmk}', 'Phase Space Density in (L, \\mu, K)', 0, 'LEA',
                            'car', ['time', 'L', 'Mu', 'K'], '1/(s*cm**2*keV*sr*MeV**2)'],
                  'L': ['L', 'L-shell', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], ''],
                  'L_2': ['L', 'L-shell', 1, 'LEA',
                        'car', ['L', 'Mu', 'K'], ''],
                  'Energy': ['E', 'Energy', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], 'MeV'],
                  'Alpha': ['Alpha', 'Pitch angle', 1, 'LEA',
                        'car', ['L', 'Energy', 'Alpha'], 'deg'],
                  'Mu': ['Mu', '1st adiabatic invariant \\mu', 1, 'LMK',
                         'car', ['L', 'Mu', 'K'], 'MeV/G'],
                  'K': ['K', '2dn adiabatic invariant K', 1, 'LMK',
                        'car', ['L', 'Mu', 'K'], 'GR_E**1/2']
                  }

def MODEL():
    from kamodo import Kamodo

    # main class
    class MODEL(Kamodo):
        '''VERB-3D model data reader.

        Inputs:
            file_dir: a string representing the file directory of the
                model output data.
            variables_requested = a list of variable name strings chosen from
                the model_varnames dictionary in this script, specifically the
                first item in the list associated with a given key.
                - If empty, the reader functionalizes all possible variables
                    (default)
                - If 'all', the reader returns the model_varnames dictionary
                    above for only the variables present in the given files.
            # filetime = boolean (default = False)
            #     - If False, the script fully executes.
            #     - If True, the script only executes far enough to determine the
            #         time values associated with the chosen data.
            # printfiles = boolean (default = False)
            #     - If False, the filenames associated with the data retrieved
            #         ARE NOT printed.
            #     - If True, the filenames associated with the data retrieved ARE
            #         printed.
            # gridded_int = boolean (default = True)
            #     - If True, the variables chosen are functionalized in both the
            #         standard method and a gridded method.
            #     - If False, the variables chosen are functionalized in only the
            #         standard method.
            # verbose = boolean (False)
            #     - If False, script execution and the underlying Kamodo
            #         execution is quiet except for specified messages.
            #     - If True, be prepared for a plethora of messages.
        All inputs are described in further detail in
            KamodoOnboardingInstructions.pdf.

        Returns: a kamodo object (see Kamodo core documentation) containing all
            requested variables in functionalized form.

        Notes:
            -
        '''
        def __init__(self, file_dir, variables_requested=[],
                     filetime=False, verbose=False, gridded_int=True,
                     printfiles=False, **kwargs):
            super(MODEL, self).__init__()
            self.modelname = 'VERB-3D'

        #


    # TODO: plt to common format conversion
    # if you need more than a few lines of code to access the data stored in the file,
    # then a file conversion routine is likely needed.
    # If it takes longer to read in the data from the current file format than from a cdf or h5 file,
    # then the file conversion step should be developed.

    return MODEL