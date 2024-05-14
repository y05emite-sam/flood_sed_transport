# -*- coding: utf-8 -*-
"""
Created on Fri May 10 12:30:12 2024


@author: Sam
"""
"""
Driver file for: 
Case 3b_1 - Test in large watershed - Uniform rainfall - bed evolution
"""

import copy
from clean_old_output import clean_old_output
from save_data_to_file import save_data
from save_raster import plot_results
from rainfall_manager import rainfall_data, update_rainfall_intensity
from update_time import update_time
import matplotlib.pyplot as plt
import numpy as np
from landlab.components import (FlowAccumulator, FlowDirectorSteepest,
                                SinkFiller, OverlandFlow, RiverBedDynamics)
from landlab.io import read_esri_ascii
from landlab.grid.mappers import map_mean_of_link_nodes_to_link
from numpy.ma import masked_array

(output_folder, output_path, cwd) = clean_old_output()

# Numerical simulation conditions and time control settings
zDEM = 'filled_lc3_dem.asc'  # ASCII raster DEM containing the bed surface elevation
rainfall_time_series = 'rainfall_time_series.txt'  # Rainfall time series
gsd = np.loadtxt('LC3_large_grain_size_dist1.txt')  # Grain size distribution data

t = 0
max_dt = 0.01  # Min time step
sim_max_t = 600 + max_dt
save_data_time_interval = 10
dt_precision = 1
plot_time_interval = 1800

# Load DEM and Set Up Grid
(grid, z) = read_esri_ascii(zDEM, name='topographic__elevation')
no_data_value = -9999
masked_elevation = masked_array(z, z == no_data_value)
grid.at_node['topographic__elevation'] = masked_elevation.data
grid.add_zeros('surface_water__depth', at='node')

# Depression Filling with SinkFiller (D4 routing)
sf = SinkFiller(grid, routing='D4')
sf.fill_pits()

# Flow Routing with D4 Routing
fd = FlowDirectorSteepest(grid)
fa = FlowAccumulator(grid, flow_director=fd)
fa.run_one_step()

# Store initial topographic elevation for later comparison
topographic__elevation_original = copy.deepcopy(grid.at_node['topographic__elevation'])

# OverlandFlow and RiverBedDynamics instantiation
of = OverlandFlow(grid, mannings_n=0.025, rainfall_intensity=0.0, alpha=0.25, steep_slopes=True, theta=1.0)
rainfall_time_series = np.loadtxt(rainfall_time_series)
(rainfall_time, rainfall_intensity, rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value) = rainfall_data(rainfall_time_series, grid)
of._rainfall_intensity = current_rainfall_intensity

rbd = RiverBedDynamics(grid, gsd=gsd, outlet_boundary_condition='fixedValue', bed_surf__elev_fix_node=np.zeros_like(z))

# Simulation loop
progress0 = 0
while t < sim_max_t:
    of.overland_flow(dt=max_dt)
    rbd.run_one_step()

    if t % plot_time_interval == 0 or plot_now:
        plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval, cwd)
        plot_topographic_change(grid, topographic__elevation_original, t, output_folder)

    (save_data_time_interval, save_data_now) = save_data(save_data_time_interval, of, dt_precision, save_data_now, t, output_path, grid, link_list, node_list, save_data_time_interval, cwd, rainfall_intensity_value)
    (t, plot_now, check_maximum_time, progress0) = update_time(t, of, sim_max_t, check_maximum_time, plot_now, progress0)

def plot_topographic_change(grid, original_elevation, current_time, output_folder):
    """ Plot the change in topography from the start of the simulation. """
    elevation_change = grid.at_node['topographic__elevation'] - original_elevation
    plt.figure(figsize=(10, 8))
    plt.imshow(elevation_change.reshape(grid.shape), cmap='RdYlGn', origin='lower')
    plt.colorbar(label='Elevation Change (m)')
    plt.title(f'Topographic Variation at t={current_time} s')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.savefig(f"{output_folder}/topographic_change_{int(current_time)}.png")
    plt.close()


