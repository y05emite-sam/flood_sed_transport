# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:53:51 2023

@author: angel
"""

# Manages data storage

import numpy as np
from os import path, chdir

def save_data_to_file(filename, timestamp, data, node_or_link_list):
    formatted_data = np.hstack([timestamp, data[node_or_link_list]])
    reshaped_data = np.reshape(formatted_data, [1, formatted_data.shape[0]])
    with open(filename, "ab") as f:
        np.savetxt(f, reshaped_data, '%.3f')
        
def save_rainfall_intensity_value(filename, timestamp, data):
    formatted_data = np.hstack([timestamp, data])
    reshaped_data = np.reshape(formatted_data, [1, formatted_data.shape[0]])
    with open(filename, "ab") as f:
        np.savetxt(f, reshaped_data, '%.3f')        

def save_data (save_data_time_interval, 
               of, 
               dt_precision, 
               save_data_now, 
               t, 
               output_path, 
               grid, 
               link_list, 
               node_list, 
               save_data_time_interval_original,
               cwd,
               rainfall_intensity_value):
    
    save_data_time_interval = round(save_data_time_interval-of.dt, dt_precision)
    if (save_data_time_interval <=0) or save_data_now:

        #print('Writing results at time :',np.round(t,1), ' \n')

        # Discharge in m3/s
        data_output_path = path.join(output_path, "output0_link_surface_water__discharge.txt")
        save_data_to_file(data_output_path, t, np.abs(of._q * grid.dx), link_list)
        
        # water depth
        data_output_path = path.join(output_path, "output1_node_surface_water__depth.txt")
        save_data_to_file(data_output_path, t, of._h, node_list)
        
        # Topographic elevation
        data_output_path = path.join(output_path, "output2_node_topographic__elevation.txt")
        save_data_to_file(data_output_path, t, grid.at_node["topographic__elevation"], node_list)
        
        # rainfall intensity
        data_output_path = path.join(output_path, "output3_rainfall_intensity.txt")
        save_rainfall_intensity_value(data_output_path, t, rainfall_intensity_value)
 
        save_data_time_interval = round(save_data_time_interval_original, dt_precision)
        save_data_now = False
        chdir(cwd)
        
    return save_data_time_interval, save_data_now