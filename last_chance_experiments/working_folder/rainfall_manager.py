# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:03:45 2023

@author: angel
"""
import numpy as np

def rainfall_data(rainfall_time_series,grid):
    
    rainfall_time = rainfall_time_series[:,0]
    rainfall_intensity = rainfall_time_series[:,1] * (2.77778 * 10 ** -7)  # rainfall in m/s
    rainfall_time_index = int(0)                                # current index for time
    rainfall_intensity_value = rainfall_intensity[rainfall_time_index] / (2.77778 * 10 ** -7)
    current_rainfall_intensity = np.full(grid.number_of_nodes, rainfall_intensity[rainfall_time_index])
    
    return rainfall_time, rainfall_intensity, rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value 
    
def update_rainfall_intensity(t, rainfall_time, rainfall_time_index, rainfall_intensity, grid, current_rainfall_intensity, rainfall_intensity_value, save_data_now):
       
    if (t >= rainfall_time[rainfall_time_index+1]):
        print('Rainfall intensity updated - Elapsed time : ', t, ' s. Current rainfall intensity :', rainfall_intensity[rainfall_time_index+1]/(2.77778 * 10 ** -7) , ' mm/hr \n' )
        rainfall_time_index += 1
        if rainfall_time_index == rainfall_time.shape[0] - 1:
            rainfall_time_index -= 1
        current_rainfall_intensity = np.full(grid.number_of_nodes, rainfall_intensity[rainfall_time_index])
        rainfall_intensity_value = rainfall_intensity[rainfall_time_index] / (2.77778 * 10 ** -7)
        save_data_now = True
       
    return rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value, save_data_now