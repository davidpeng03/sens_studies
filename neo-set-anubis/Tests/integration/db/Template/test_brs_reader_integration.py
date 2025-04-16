from db.Template.brs_reader import BRParser
import unittest

class TestIntegrationBRParser(unittest.TestCase):
    def test_full_integration(self):
        parser = BRParser()
        with open("test_eq.txt", "w") as f:
            f.write("x + y + z")

        func = parser.parse_equation("test_eq.txt")
        self.assertEqual(func(1, 2, 3), 6)