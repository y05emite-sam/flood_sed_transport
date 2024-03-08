# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:55:01 2023

@author: angel
"""


from os import chdir
import numpy as np
from matplotlib import pyplot as plt

# Defines a function to plot results as raster
def save_raster(data,filename,grid):
    # Reshape and flip data
    reshaped_data = np.flip(np.reshape(data, [grid.number_of_node_rows, grid.number_of_node_columns]), axis=0)

    # Prepare metadata lines
    lines = [
        f'nrows {grid.number_of_node_rows}',
        f'cellsize {grid.dx}',
        f'xllcorner {grid.origin[0]}',
        f'ncols {grid.number_of_node_columns}',
        f'yllcorner {grid.origin[1]}',
        'nodata_value -9999'
    ]
    
    # Write metadata and reshaped data to file
    with open(filename, 'w') as f:
        f.write('\n'.join(lines) + '\n')
        np.savetxt(f, reshaped_data, '%.3f')
        
def plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval_original, cwd):

    plot_time_interval = round(plot_time_interval-of.dt, dt_precision)
    if plot_time_interval <= 0  or plot_now:
        chdir(output_folder)
        print('Elapsed time :',np.round(t,1),' s. Current dt =', np.round(of.dt,2),' s - Saving plot \n')
        
        # Water depth raster
        filename='depth_'+str(np.round(t,0))+'.asc'    
        data = grid.at_node["surface_water__depth"]
        save_raster(data,filename,grid)
               
        #Bed surface elevation plot
        filename='topographicElevation_'+str(np.round(t,0))+'.asc'    
        data = grid.at_node["topographic__elevation"]
        save_raster(data,filename,grid)
        
        #Bed surface variation plot
        filename = 'topographicVariation_'+str(np.round(t,0))+'.asc'    
        data = grid.at_node["topographic__elevation"] - topographic__elevation_original 
        save_raster(data,filename,grid)
    
        time, q0, q1, q2 , q3, q4, q5, q6, q7, q8, q9, q10  = np.loadtxt('output0_link_surface_water__discharge.txt', delimiter=' ', unpack=True)
        
        # Plotting the lines
        #plt.plot(time, q0, label='1500', color='blue', marker='o', linestyle='-', linewidth=2)
        #plt.plot(time, q1, label='2099', color='blue', marker='*', linestyle='-', linewidth=2)
        #plt.plot(time, q2, label='1499', color='blue', marker='s', linestyle='-', linewidth=2)
        #plt.plot(time, q3, label='900',  color='blue', marker='^', linestyle='-', linewidth=2)
        plt.plot(time, q3, label='900',  color='blue', marker='none', linestyle='-', linewidth=2)
        #plt.plot(time, q4, label='5097', color='orange', marker='o', linestyle='-', linewidth=2)
        #plt.plot(time, q5, label='5696', color='orange', marker='*', linestyle='-', linewidth=2)
        #plt.plot(time, q6, label='5096', color='orange', marker='s', linestyle='-', linewidth=2)
        #plt.plot(time, q7, label='4497', color='orange', marker='^', linestyle='-', linewidth=2)
        plt.plot(time, q7, label='4497', color='orange', marker='none', linestyle='-', linewidth=2)
        #plt.plot(time, q8, label='12291', color='green', marker='o', linestyle='-', linewidth=2)
        #plt.plot(time, q9, label='12890', color='green', marker='*', linestyle='-', linewidth=2)
        #plt.plot(time, q10, label='12290', color='green', marker='s', linestyle='-', linewidth=2)
        #plt.plot(time, q11, label='11691', color='green', marker='^', linestyle='-', linewidth=2)

        
        # Adding legend
        plt.legend()
        
        # Adding title and labels
        plt.title('Discharge Over Time')
        plt.xlabel('Time [s]')
        plt.ylabel('Discharge [m^3/s]')
        plt.xlim([0,86400])
        plt.ylim([0,80])
        
        # Showing the plot
        plt.show() 
    
        plot_now = False
        plot_time_interval = plot_time_interval_original
        chdir(cwd)

    return plot_time_interval, plot_now