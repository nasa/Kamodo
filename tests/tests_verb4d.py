import unittest
from unittest import TestCase
from unittest.mock import patch

from io import StringIO
import os
import numpy as np

import kamodo_ccmc.flythrough.model_wrapper as MW

class TestVerb01BasicSetup(unittest.TestCase):
    """ This class checks basic setup of the VERB Kamodo model """

    # Static model name
    model = 'VERB-3D'



    @patch('sys.stdout', new_callable=StringIO)
    def test00_Choose_Model_Global(self, mock_stdout):
        """ Test that VERB is present in the global list of models"""
        # List all models
        MW.Choose_Model()
        # Get the entire output
        output = mock_stdout.getvalue()
        self.assertIn("VERB-3D", output, "Model name VERB-3D is missing")

    def test00_Choose_Model(self):
        """ Confirming the model existence """
        model = MW.Choose_Model(self.model)
        self.assertEqual(model.MODEL().modelname, self.model, "Model class does not have correct modelname property")

    def test00_Model_Reader(self):
        """ Confirming the model existence """
        model = MW.Model_Reader(self.model)
        self.assertEqual(model.modelname, self.model, "Model class does not have correct modelname property")


    @patch('sys.stdout', new_callable=StringIO)
    def test01_Variable_Search_Global(self, mock_stdout):
        """ Test that VERB is present in the global list of variables"""
        # List all variables
        MW.Variable_Search('')

        # Get the entire output
        output = mock_stdout.getvalue()

        self.assertIn("VERB-3D", output, "Model name VERB-3D is missing")
        self.assertIn("PSD", output, "Variable PSD is missing")
        self.assertIn("Phase Space Density", "Variable description (Phase Space Density) is missing")


    @patch('sys.stdout', new_callable=StringIO)
    def test03_Variable_Search(self, mock_stdout):
        """ Test for specific variables when model is specified"""
        # List all variables
        MW.Variable_Search('', self.model)

        # Get the entire output
        output = mock_stdout.getvalue()

        self.assertIn("PSD", output, "Variable PSD is missing")
        self.assertIn("Phase Space Density", "Variable description (Phase Space Density) is missing")

    @patch('sys.stdout', new_callable=StringIO)
    def test04_Variable_Search_Local(self, mock_stdout):
        MW.Variable_Search('Phase Space Density', self.model)

        # Get the entire output
        output = mock_stdout.getvalue()

        self.assertIn("PSD", output, "Variable PSD is missing")
        self.assertIn("Phase Space Density", "Variable description (Phase Space Density) is missing")

    def test05_Ver_3D(self):
        """ Test if the model has expected 3D variables description """
        expected_output = [
            'L-shell',
            'Electron energy',
            'Equatorial pitch angle',
            'Momentum times speed of light',
            '1st adiabatic invariant mu',
            '2dn adiabatic invariant K'
        ]

        result = MW.Var_3D(self.model)

        for item in expected_output:
            self.assertIn(item, result, f"{item} is missing in the result")

    def test06_Ver_units(self):
        def test_Var_units(self):
            variables_requested = [
                'PSD_lea', 'flux_lea', 'L', 'E_e', 'alpha_e', 'pc', 'PSD_lmk', 'mu', 'K'
            ]
            expected_output = {
                'PSD_lea': '(c/MeV/cm)**3',
                'flux_lea': '1/(s*cm**2*keV*sr)',
                'L': '',
                'E_e': 'MeV',
                'alpha_e': 'deg',
                'pc': 'MeV',
                'PSD_lmk': '(c/MeV/cm)**3',
                'mu': 'MeV/G',
                'K': '10**-4*T*(m/m)**1/2'
            }

            result = MW.Var_units(self.model, ['Trash', variables_requested])

            for var in expected_output:
                self.assertIn(var, result, f'{var} is missing')
                # For now, we just test that units are present. We do not test the unites themselves yet.
                self.assertIsInstance(result[var], str, f'Error with {var} units')

            # Just test a single unit
            self.assertIsInstance(result['E_e'], expected_output['E_e'], f'Error with E_e units values')

            # Ensure 'Trash' is not in the result
            self.assertNotIn('Trash', result, 'Unexpected variable is returned')

class FakeDataGenerator:
    """ This class generates fake VERB output """
    output_dir = "output"
    output_path = output_dir + '/'  # because Kamodo cannot work without slash

    @staticmethod
    def generate_fake_out1d(file_path):
        header = 'Variables = "time", "Kp", "Boundary_fluxes", "Lpp", "tau", "DCT_Daa WT_CHORUS_DAY_2016", "DCT_Daa WT_CHORUS_NIGHT_2016", "DCT_Dpp WT_CHORUS_DAY_2016", "DCT_Dpp WT_CHORUS_NIGHT_2016", "DCT_Dpa WT_CHORUS_DAY_2016", "DCT_Dpa WT_CHORUS_NIGHT_2016", "DCT_DaaLpp WT_HISS_2016", "DCT_DppLpp WT_HISS_2016", "DCT_DpaLpp WT_HISS_2016"\n'
        zone = 'ZONE T="1d-output"\n'
        time_steps = np.linspace(0, 4, 5)
        data = np.random.random((5, 13))

        with open(file_path, 'w') as file:
            file.write(header)
            file.write(zone)
            for i, time in enumerate(time_steps):
                file.write(f'{time}\t' + '\t'.join(map(str, data[i])) + '\n')

    @staticmethod
    def generate_fake_outpsd(file_path):
        header = 'VARIABLES = "Phase_Space_Density"\n'
        zones = 5
        data_shape = (4, 5, 6)

        with open(file_path, 'w') as file:
            file.write(header)
            for t in range(zones):
                file.write(f'ZONE T="{t}" I={data_shape[2]}, J={data_shape[1]}, K={data_shape[0]}\n')
                data = np.random.random(np.prod(data_shape))
                for value in data:
                    file.write(f'{value:e}\n')

    @staticmethod
    def generate_fake_perp_grid(file_path):
        header = 'VARIABLES = "L, Earth radius", "Energy, MeV", "Pitch-angle, deg", "pc"\nZONE T="Grid"  I=6, J=5, K=4\n'
        L_values = np.linspace(1, 4, 4)
        E_values = 10.**np.linspace(-2, 2, 5)
        A_values = np.linspace(0.1, np.pi/2-0.1, 6)
        pc_values = np.random.random((4, 5, 6))

        with open(file_path, 'w') as file:
            file.write(header)
            for i, L in enumerate(L_values):
                for j, E in enumerate(E_values):
                    for k, A in enumerate(A_values):
                        pc = pc_values[i, j, k]
                        file.write(f'{L:e}\t{E:e}\t{A:e}\t{pc:e}\n')

    @staticmethod
    def setup_fake_data():
        os.makedirs(FakeDataGenerator.output_dir, exist_ok=True)
        FakeDataGenerator.generate_fake_out1d(os.path.join(FakeDataGenerator.output_dir, 'out1d.dat'))
        FakeDataGenerator.generate_fake_outpsd(os.path.join(FakeDataGenerator.output_dir, 'OutPSD.dat'))
        FakeDataGenerator.generate_fake_perp_grid(os.path.join(FakeDataGenerator.output_dir, 'perp_grid.plt'))

    @staticmethod
    def teardown_fake_data():
        for file_name in os.listdir(FakeDataGenerator.output_dir):
            file_path = os.path.join(FakeDataGenerator.output_dir, file_name)
            os.remove(file_path)
        os.rmdir(FakeDataGenerator.output_dir)

# TODO: Create class that will test the nc files creating with diffrent methods:
# def Model_Variables(model, file_dir=None, return_dict=False):
# def Variable_Search(search_string='', model='', file_dir='', return_dict=False):
# def File_Times(model, file_dir, print_output=True):
# def File_List(model, file_dir, print_output=False):
class TestVerb03DatasetGeneration(TestCase):
    def setUp(self):
        pass

# TODO: Make this class all about testing the actual data?
class TestVerb03DatasetCheck(TestCase):
    """ This class is advanced test with fake dataset """

    model = 'VERB-3D'
    teardown = False

    @classmethod
    def setUpClass(cls):
        if os.path.exists(FakeDataGenerator.output_dir) and cls.teardown:
            raise unittest.SkipTest("Test cannot be done because the folder 'output' already exists.")
        FakeDataGenerator.setup_fake_data()

    @classmethod
    def tearDownClass(cls):
        if cls.teardown:
            FakeDataGenerator.teardown_fake_data()

    def test01_verb_dataset(self):
        pass

    @patch('sys.stdout', new_callable=StringIO)
    def test02_File_Times(self, mock_stdout):
        MW.File_Times(self.model, FakeDataGenerator.output_path)

    def test_generate_fake_files(self):
        out1d_path = os.path.join(FakeDataGenerator.output_dir, 'out1d.dat')
        outpsd_path = os.path.join(FakeDataGenerator.output_dir, 'OutPSD.dat')
        perp_grid_path = os.path.join(FakeDataGenerator.output_dir, 'perp_grid.plt')

        # Verify files are created
        self.assertTrue(os.path.isfile(out1d_path))
        self.assertTrue(os.path.isfile(outpsd_path))
        self.assertTrue(os.path.isfile(perp_grid_path))

# TODO: Make this class checking the output. use plotpy generator to control the output somehow.
# class TestVerb04plot(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()
