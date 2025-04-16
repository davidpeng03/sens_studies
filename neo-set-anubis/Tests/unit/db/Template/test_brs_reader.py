import unittest
from unittest.mock import MagicMock, patch
from db.Template.brs_reader import BRParser, ExpressionParser, FileReader
import sympy as sp
from abc import ABC, abstractmethod
from typing import List, Tuple, Callable

class TestBRParser(unittest.TestCase):
    def setUp(self):
        self.parser = BRParser()

    def test_parse_equation_valid_file(self):
        mock_file_content = "x + y + z"
        with patch("db.Template.brs_reader.FileReader.read_file", return_value=mock_file_content):
            func = self.parser.parse_equation("dummy_path.txt")
            self.assertEqual(func(1, 2, 3), 6)

    def test_parse_equation_file_not_found(self):
        with patch("db.Template.brs_reader.FileReader.read_file", side_effect=FileNotFoundError):
            with self.assertRaises(FileNotFoundError):
                self.parser.parse_equation("nonexistent.txt")

    def test_extract_variables(self):
        expression = "x**2 + y"
        variables, expr = ExpressionParser.extract_variables(expression)
        self.assertEqual(variables, [sp.Symbol("x"), sp.Symbol("y")])
        self.assertEqual(str(expr), "x**2 + y")

    def test_generate_function(self):
        variables = [sp.Symbol("x"), sp.Symbol("y")]
        expression = "x**2 + y"
        func = ExpressionParser.generate_function(expression, variables)
        self.assertEqual(func(2, 3), 7)

    def test_parse_empty(self):
        self.assertEqual(self.parser.parse(), [])