import unittest
import os
import json
import sys
sys.path.append("")  # Adjust the path as necessary
from Core.Paramaters import SimulationParameters

class TestSimulationParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_db_file = 'Tests/test_db/test_db.json'
        if os.path.exists(cls.test_db_file):
            os.remove(cls.test_db_file)
        cls.sim_params = SimulationParameters(db_file=cls.test_db_file)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.test_db_file):
            os.remove(cls.test_db_file)

    def test_load_from_particle_package(self):
        self.sim_params.load_from_particle_package()
        self.assertIn(211, self.sim_params.params['mass'])
        self.assertIn(211, self.sim_params.params['pdg'])

    def test_save_to_file(self):
        self.sim_params.save_to_file()
        self.assertTrue(os.path.exists(self.test_db_file))

    def test_load_from_file(self):
        self.sim_params.save_to_file()
        self.sim_params.params = {}
        self.sim_params.load_from_file()
        self.assertIn(211, self.sim_params.params['mass'])
        self.assertIn(211, self.sim_params.params['pdg'])

    def test_get_parameter(self):
        self.sim_params.set_parameter('mass', 211, 0.13957)
        param = self.sim_params.get_parameter('mass', 211)
        self.assertIsNotNone(param)

    def test_set_parameter(self):
        self.sim_params.set_parameter('mass', 999, 1.234)
        param = self.sim_params.get_parameter('mass', 999)
        self.assertEqual(param, 1.234)

    def test_add_parameter(self):
        with self.assertRaises(KeyError):
            self.sim_params.add_parameter('mass', 211, 0.13957)  # Assuming it already exists

        self.sim_params.add_parameter('mass', 888, 0.888)
        param = self.sim_params.get_parameter('mass', 888)
        self.assertEqual(param, 0.888)

    def test_update_parameter(self):
        with self.assertRaises(KeyError):
            self.sim_params.update_parameter('mass', 777, 0.777)  # Assuming it does not exist

        self.sim_params.add_parameter('mass', 666, 0.666)
        self.sim_params.update_parameter('mass', 666, 0.999)
        param = self.sim_params.get_parameter('mass', 666)
        self.assertEqual(param, 0.999)

if __name__ == '__main__':
    unittest.main()


