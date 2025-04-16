import sympy as sp
import os,sys
from typing import List, Tuple, Callable, Dict
sys.path.append(os.getcwd())
from db.Template.parser import IParser
class FileReader:
    """Utility class for reading file content."""

    @staticmethod
    def read_file(filename: str) -> str:
        """
        Read and return the content of a file.

        Args:
            filename (str): Path to the file.

        Returns:
            str: Content of the file.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        try:
            with open(filename, 'r') as file:
                return file.read().strip()
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filename}")


class ExpressionParser:
    """Utility class for parsing and evaluating mathematical expressions."""

    @staticmethod
    def extract_variables(expression: str) -> Tuple[List[sp.Symbol], sp.Basic]:
        """
        Extract variables and sympy expression from a string.

        Args:
            expression (str): Mathematical expression as a string.

        Returns:
            Tuple[List[sp.Symbol], sp.Basic]: List of variables and sympy expression.
        """
        expr = sp.sympify(expression)
        variables = sorted(expr.free_symbols, key=lambda x: str(x))
        return variables, expr

    @staticmethod
    def generate_function(expression: str, variables: List[sp.Symbol]) -> Callable:
        """
        Generate a function from an expression and variables that accepts a dictionary as input.

        Args:
            expression (str): Mathematical expression as a string.
            variables (List[sp.Symbol]): List of variables.

        Returns:
            Callable: A function that evaluates the expression using a dictionary.
        """
        expr = sp.sympify(expression)
        var_list = [sp.symbols(str(var)) for var in variables]
        func = sp.lambdify(var_list, expr, modules="numpy")

        def wrapped_func(values: Dict[str, float]) -> float:
            """
            Evaluates the generated function using a dictionary.

            Args:
                values (Dict[str, float]): Dictionary mapping variable names to values.

            Returns:
                float: Evaluated result of the expression.
            """
            args = [values[str(var)] for var in var_list]
            return func(*args)

        return wrapped_func


class BRParser(IParser):
    """
    Parser for branching ratio (BR) equations stored in files.

    Inherits:
        IParser: Abstract base class for parsers.
    """

    def __init__(self, file_reader: FileReader = FileReader, parser: ExpressionParser = ExpressionParser):
        """
        Initialize the BRParser.

        Args:
            file_reader (FileReader): Utility for reading files. Defaults to FileReader.
            parser (ExpressionParser): Utility for parsing expressions. Defaults to ExpressionParser.
        """
        self.file_reader = file_reader
        self.parser = parser
        self.filename = None
        self.variables, self.expression = None, None

    def parse_equation(self, filename: str) -> Callable:
        """
        Parse a file containing a mathematical equation and return a callable function.

        Args:
            filename (str): Path to the file.

        Returns:
            Callable: A function that evaluates the equation.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        self.filename = filename
        content = self.file_reader.read_file(filename)
        self.variables, self.expression = self.parser.extract_variables(content)
        return self.parser.generate_function(self.expression, self.variables)

    def parse(self) -> list:
        """
        Abstract method implementation. Returns an empty list as a placeholder.

        Returns:
            list: Parsed data.
        """
        # Implement logic to parse multiple equations if needed.
        return []
    
    def _extract_variables(self):
        vars = [str(x) for x in self.variables]
        return vars

    def _extract_eq(self):
        return self.expression
    
if __name__ == "__main__":
    BRparsers = BRParser()
    file = BRparsers._parse_file("db/HNL/Prod_brs/Prod_br_H.txt")
    print(file)
    var, exp = BRparsers._extract_variables(file)
    print(var, exp)
    func = BRparsers._generate_function(exp, var)
    print(func)
    print(func(1,1,2,3))
    