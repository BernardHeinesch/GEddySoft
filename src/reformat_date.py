from datetime import datetime
import re


def reformat_date(date_str, input_format_str, output_format_str):
    # Define a dictionary to map the custom format to strftime/strptime formats
    format_mapping = {
        'yyyy': '%Y',
        'yy': '%y',
        'mm': '%m',
        'dd': '%d',
        'HH': '%H',
        'MM': '%M',
        'SS': '%S'
    }
    
    # Helper function to convert custom format to strftime/strptime format
    def convert_format(custom_format):
        # Match all format parts (e.g., yyyy, mm, dd, HH, MM, SS)
        format_parts = re.findall(r'(yyyy|yy|mm|dd|HH|MM|SS)', custom_format)
        date_format = custom_format
        for part in format_parts:
            date_format = date_format.replace(part, format_mapping[part])
        return date_format
    
    # Convert input and output formats to strftime/strptime formats
    input_date_format = convert_format(input_format_str)
    output_date_format = convert_format(output_format_str)
    
    # Convert the input date string to a datetime object
    date_obj = datetime.strptime(date_str, input_date_format)
    
    # Reformat the date according to the output format provided by the user
    reformatted_date = date_obj.strftime(output_date_format)
    
    return reformatted_date

