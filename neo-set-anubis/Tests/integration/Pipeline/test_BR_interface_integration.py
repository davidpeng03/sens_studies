import unittest
import sys, os
sys.path.append(os.getcwd())
from Pipeline.BR_calculator import BRInterface
import numpy as np

class TestBRInterfaceIntegration(unittest.TestCase):

    def test_integration_with_python_strategy(self):
        br_interface = BRInterface()

        br_interface.set_calculation_method('Python', True, False)
        br_interface.set_model('HNL')
        br_interface.set_params({"Ve": 1, "Vmu": 1, "Vta": 1})
        br_interface.set_masses({"N1": 1})

        result_decay_tot = br_interface.calculate("DecayTot", "N1")
        self.assertIsNotNone(result_decay_tot)
        print("result decay tot", result_decay_tot)
        self.assertIsInstance(result_decay_tot, np.float64)

        result_br = br_interface.calculate("BR", "N1", channel=(211, 13))
        self.assertIsNotNone(result_br)
        self.assertIsInstance(result_br, float)

        result_decay_tot_with_channel = br_interface.calculate("DecayTot", "N1", channel=(211, 13))
        print("result with problem", result_decay_tot_with_channel, type(result_decay_tot_with_channel))
        self.assertIsNotNone(result_decay_tot_with_channel)
        self.assertIsInstance(result_decay_tot_with_channel, np.float64)

        result_prod_br = br_interface.calculate("ProdBR", "N1", mother_particle=24)
        self.assertIsNotNone(result_prod_br)
        self.assertIsInstance(result_prod_br, float)

    def test_integration_with_file_strategy(self):
        br_interface = BRInterface()

        br_interface.set_calculation_method('File')
        br_interface.set_model('HNL')
        br_interface.set_params({"Ve": 1, "Vmu": 1, "Vta": 1})
        br_interface.set_masses({"N1": 1})

        result_decay_tot = br_interface.calculate("DecayTot", "N1")
        self.assertIsNotNone(result_decay_tot)
        print("result decay tot with file", result_decay_tot, type(result_decay_tot))
        self.assertIsInstance(result_decay_tot, float)

        result_br = br_interface.calculate("BR", "N1", channel=(11,-11,12))
        self.assertIsNotNone(result_br)
        self.assertIsInstance(result_br, float)

        result_decay_tot_with_channel = br_interface.calculate("DecayTot", "N1", channel=(211, 13))
        self.assertIsNotNone(result_decay_tot_with_channel)
        self.assertIsInstance(result_decay_tot_with_channel, float)

        result_prod_br = br_interface.calculate("ProdBR", "N1", mother_particle=24)
        self.assertIsNotNone(result_prod_br)
        self.assertIsInstance(result_prod_br, float)

if __name__ == '__main__':
    unittest.main()
