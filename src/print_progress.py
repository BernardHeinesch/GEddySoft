# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:00:54 2022

@author: Ariane Faures
"""

def print_progress(processed_, total_to_process, text, flush_=True, end_='\r',
                   width=30):
    """
    Print progress of process.

    Parameters
    ----------
    processed_ : Int
        Number of files/data processed.
    total_to_process : Int
        Total number of files/data to process.
    text : String
        Text to display.
    width : int, optional
        The default is 30. Width of the progressing bar

    Returns
    -------
    percent_processed : float
        Percentage processed
    """

    percent_processed = processed_ / total_to_process * 100
    left = width * int(percent_processed) // 100
    right = width - left
    print('\r' + text + ': [',
          '#' * left, ' ' * right, ']',
          f' {percent_processed:.0f}%',
          sep='', end=end_, flush=flush_)
    return percent_processed
