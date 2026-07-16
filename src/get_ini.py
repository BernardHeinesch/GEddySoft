import configparser
import ast


class MyConfigParser(configparser.ConfigParser):
    def get(self, section, option, *, raw=False, vars=None, fallback=configparser._UNSET):
        value = super().get(section, option, raw=raw, vars=vars, fallback=fallback)
        try:
            # Try to safely evaluate the value using ast.literal_eval
            return ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # If evaluation fails, return the original string
            return value


def get_ini(filename):
    """
    read the ini file using configparser, but end with a dictionary
    with the right data types, not all strings

    parameters
    ----------
    filename: string, the finename of the ini file

    returns
    -------
    ini: dictionnary with all initialisation information

    comments
    --------

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    # Create an instance of the custom ConfigParser
    ini = MyConfigParser()

    # Read the INI file
    ini.read(filename)

    # Return the MyConfigParser instance
    return ini
