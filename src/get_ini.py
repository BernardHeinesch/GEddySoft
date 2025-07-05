"""Module for intelligent INI file parsing with automatic type conversion.

This module extends Python's configparser to provide automatic type inference
for configuration values. It supports:

- Standard INI file parsing
- Automatic conversion of strings to appropriate Python types
- Safe evaluation of complex data structures (lists, dicts, etc.)
- Fallback to string values when type conversion fails

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import configparser
import ast


class MyConfigParser(configparser.ConfigParser):
    """Enhanced ConfigParser with automatic type inference.

    This class extends the standard ConfigParser to automatically convert
    string values to appropriate Python types using ast.literal_eval.
    It safely handles both simple types (int, float, bool) and complex
    types (lists, dicts, tuples).

    The conversion is attempted using ast.literal_eval, which safely
    evaluates strings containing Python literals. If conversion fails,
    the original string value is returned unchanged.
    """
    def get(self, section, option, *, raw=False, vars=None, fallback=configparser._UNSET):
        """Get an option value with automatic type conversion.

        This method extends the standard ConfigParser.get() by attempting
        to convert string values to appropriate Python types.

        Parameters
        ----------
        section : str
            Section name in the configuration
        option : str
            Option name in the specified section
        raw : bool, optional
            If True, no interpolation is performed
        vars : dict, optional
            Dictionary of substitution variables
        fallback : any, optional
            Value to return if the option is not found

        Returns
        -------
        any
            The option value converted to its appropriate Python type,
            or the original string if conversion fails
        """
        value = super().get(section, option, raw=raw, vars=vars, fallback=fallback)
        try:
            # Try to safely evaluate the value using ast.literal_eval
            return ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # If evaluation fails, return the original string
            return value


def get_ini(filename):
    """Read and parse an INI file with automatic type conversion.

    This function reads an INI configuration file and returns a parser
    that automatically converts values to appropriate Python types.
    It uses MyConfigParser to handle type inference, making it ideal
    for scientific applications with complex configuration needs.

    Parameters
    ----------
    filename : str
        Path to the INI configuration file

    Returns
    -------
    MyConfigParser
        Configured parser instance with loaded and type-converted values

    Notes
    -----
    The function supports all standard INI file features plus:
    - Section-based organization
    - Key-value pairs with automatic type conversion
    - Complex data structures (lists, dicts, etc.)
    - Fallback to string values when conversion fails
    """


    # Create an instance of the custom ConfigParser
    ini = MyConfigParser()

    # Read the INI file
    ini.read(filename)

    # Return the MyConfigParser instance
    return ini
