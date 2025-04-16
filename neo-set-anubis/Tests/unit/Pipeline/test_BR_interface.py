import unittest
import sys, os
sys.path.append(os.getcwd())
from unittest.mock import MagicMock, patch
from Pipeline.BR_calculator import BRInterface
from Pipeline.Template.BR_strategies import PythonCalculationStrategy, FileCalculationStrategy

class TestBRInterface(unittest.TestCase):

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_calculation_method_python(self, MockPythonStrategy):
        br_interface = BRInterface()

        # Test setting Python calculation method
        br_interface.set_calculation_method('Python', True, False)
        self.assertIsInstance(br_interface.calculation_strategy, MockPythonStrategy)
        MockPythonStrategy.assert_called_once_with(True, False)

    @patch('Pipeline.Template.BR_strategies.FileCalculationStrategy')
    def test_set_calculation_method_file(self, MockFileStrategy):
        br_interface = BRInterface()

        # Test setting File calculation method
        br_interface.set_calculation_method('File')
        self.assertIsInstance(br_interface.calculation_strategy, MockFileStrategy)

    def test_set_calculation_method_invalid(self):
        br_interface = BRInterface()

        # Test setting an invalid calculation method
        with self.assertRaises(ValueError):
            br_interface.set_calculation_method('InvalidMethod')

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_model(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        # Mock the set_model method of the strategy
        br_interface.calculation_strategy.set_model = MagicMock()

        # Call set_model on the interface
        br_interface.set_model('HNL')

        # Verify that set_model was called on the strategy with the correct argument
        br_interface.calculation_strategy.set_model.assert_called_once_with('HNL')

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_model_without_strategy(self, MockPythonStrategy):
        br_interface = BRInterface()

        # Test setting model without setting calculation strategy first
        with self.assertRaises(Exception):
            br_interface.set_model('HNL')

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_params(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        mock_params = {'Ve': 1, 'Vmu': 1, 'Vta': 1}
        br_interface.calculation_strategy.set_params = MagicMock()

        br_interface.set_params(mock_params)

        br_interface.calculation_strategy.set_params.assert_called_once_with(mock_params)

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_one_param(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        br_interface.calculation_strategy.set_one_param = MagicMock()

        br_interface.set_one_param('Ve', 1)

        br_interface.calculation_strategy.set_one_param.assert_called_once_with('Ve', 1)

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_set_masses(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        mock_masses = {'N1': 1}
        br_interface.calculation_strategy.set_masses = MagicMock()

        br_interface.set_masses(mock_masses)

        br_interface.calculation_strategy.set_masses.assert_called_once_with(mock_masses)

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_get_params(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        mock_params = {'Ve': 1, 'Vmu': 1, 'Vta': 1}
        br_interface.calculation_strategy.get_params = MagicMock(return_value=mock_params)

        result = br_interface.get_params()

        br_interface.calculation_strategy.get_params.assert_called_once()
        self.assertEqual(result, mock_params)

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_calculate(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        br_interface.calculation_strategy.calculate = MagicMock(return_value=42)

        result = br_interface.calculate('DecayTot', 'N1')

        br_interface.calculation_strategy.calculate.assert_called_once_with('DecayTot', 'N1', channels=None, mother_particle=None)
        self.assertEqual(result, 42)

    @patch('Pipeline.Template.BR_strategies.PythonCalculationStrategy')
    def test_calculate_with_channel_and_mother_particle(self, MockPythonStrategy):
        br_interface = BRInterface()
        br_interface.set_calculation_method('Python', True, False)

        br_interface.calculation_strategy.calculate = MagicMock(return_value=42)

        result = br_interface.calculate('ProdBR', 'N1', channel=(211, 13), mother_particle=24)

        br_interface.calculation_strategy.calculate.assert_called_once_with('ProdBR', 'N1', channels=(211, 13), mother_particle=24)
        self.assertEqual(result, 42)

if __name__ == '__main__':
    unittest.main()
