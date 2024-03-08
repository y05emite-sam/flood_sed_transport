# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:41:10 2023

@author: angel
"""

from os import getcwd, path, mkdir
from shutil import rmtree

def clean_old_output():
    """ Clean old files and the output folder before starting the simulation """
    
    output_folder = 'output'
    cwd = getcwd()
    
    output_path = path.join(cwd, 'output')
    if path.exists(output_path):
        print(f'The folder {output_folder} exists in the working folder and will be removed.\n')
        rmtree(output_path)
    mkdir(output_path)
    
    return output_folder, output_path, cwd