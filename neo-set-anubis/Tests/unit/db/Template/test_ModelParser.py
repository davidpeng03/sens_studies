import unittest
from unittest.mock import patch, mock_open
from db.Template.ModelParser import ModelParser

class TestModelParser(unittest.TestCase):
    def setUp(self):
        self.parser = ModelParser()

    def test_parse_file(self):
        mock_data = "::Particle1\n(1,2): yes # comment\n:ProductionMode\nprod1\n"
        with patch("builtins.open", mock_open(read_data=mock_data)):
            self.parser.parse("dummy_path.txt")
            self.assertIn("Particle1", self.parser.get_new_particles())
            self.assertEqual(self.parser.get_channels("Particle1"), [(1, 2)])
            self.assertEqual(self.parser.get_prod("Particle1"), ["prod1"])

    def test_modify_channel(self):
        mock_data = "::Particle1\n(1,2): yes # comment\n:ProductionMode\nprod1\n"
        with patch("builtins.open", mock_open(read_data=mock_data)) as mock_file:
            self.parser.parse("dummy_path.txt")
            self.parser.modify_channel("Particle1", (1, 2), "no")
            mock_file().write.assert_called()

    def test_get_new_particles(self):
        mock_data = "::Particle1\n(1,2): yes # comment\n"
        with patch("builtins.open", mock_open(read_data=mock_data)):
            self.parser.parse("dummy_path.txt")
            self.assertListEqual(self.parser.get_new_particles(), ["Particle1"])