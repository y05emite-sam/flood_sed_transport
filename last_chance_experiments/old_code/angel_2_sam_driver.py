# -*- coding: utf-8 -*-
"""

@author: Sam
"""

""" Loads components and other libraries"""
#%reset -f

import copy
from clean_old_output import clean_old_output
from save_data_to_file import save_data
from save_raster import plot_results
from rainfall_manager import rainfall_data, update_rainfall_intensity
from update_time import update_time
import pandas as pd
import numpy as np
from landlab.components import OverlandFlow, RiverBedDynamics
from landlab.io import read_esri_ascii
from landlab.grid.mappers import map_mean_of_link_nodes_to_link

(output_folder, output_path, cwd) = clean_old_output()

""" Numerical simulation conditions and time control settings"""
# inputs filenames
zDEM = 'DEMs/lc3_dem.txt'                                # ASCII raster DEM containing the bed surface elevation
rainfall_time_series = pd.read_excel('C:/Users/Sam/MultiGrain/last_chance_experiments/climate/5min_1000yrRI_storm.xlsx')

# Update these lines with the correct column names
time_series = rainfall_time_series['Time [s]'].values
precipitation = rainfall_time_series['Precipitation [mm/hr]'].values

gsd = pd.read_excel('GSDs/LC3_grain_size_dist.xlsx', sheet_name='GSD', skiprows=0).values
                           # Check inside the txt file for more information on the file format

# Time variables
t = 0
max_dt = 0.25                          # Overland flow will use the min time step between this value and the automatically calculated. Use seconds.

# Simulation conditions settings
sim_max_t = 24*3600+max_dt          # Maximum simulation time
save_data_time_interval = 5         # Stores results every this time
dt_precision = 3                    # Avoids rounding errors
plot_time_interval = 1800           # Plots will be obtained every this seconds

# General simulation conditions
n = 0.025                                           # Manning's n
fixed_nodes_id = np.array((300,301,302,900,901,902,1500,1501,1502,2100,2101,2102,2700,2701,2702,3300,3301,3302))                    # Node ID for fixed Nodes - These will not change elevation
link_list = np.array([1500, 2099,1499, 900, 5097, 5696, 5096, 4497, 12291, 12890, 12290, 11691, 23085, 23684, 23084, 22485, 289266, 289865, 289265, 288666 ])    # Location of the links where information for postprocess will be extracted 
node_list = np.array([2701 ,6301, 144907])               # Location of the nodes where information for postprocess will be extracted 

""" OverlandFlow instantiation conditions """       
# Creates fields and instantiate the OverlandFlow component
OverlandFlow.input_var_names                                                # Gives the list of all required fields                                  
(grid, z) = read_esri_ascii(zDEM, name='topographic__elevation')             # Creates the topographic__elevation field
grid.add_zeros('surface_water__depth', at = 'node')                          # Creates the surface_water__depth field
of = OverlandFlow(grid, mannings_n=n, rainfall_intensity=0.0, alpha = 0.25, steep_slopes = False, theta = 1.0)    # instantiate the Overland flow component

(rainfall_time, rainfall_intensity, rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value) = rainfall_data(time_series, precipitation, grid)
of._rainfall_intensity = current_rainfall_intensity

""" RiverBedDynamics intantiation conditions """
# surface_water__depth and topographic__elevation were already created 
RiverBedDynamics.input_var_names #('surface_water__depth', 'surface_water__velocity', 'topographic__elevation')
RiverBedDynamics.var_mapping     # Tells where those fields have to be mapped
grid["link"]["surface_water__depth"] = map_mean_of_link_nodes_to_link(grid, "surface_water__depth") # surface_water__depth was created at nodes
grid.add_zeros('surface_water__velocity', at = 'node')
grid.add_zeros('surface_water__velocity', at = 'link')
grid.set_watershed_boundary_condition(z) # Set boundaries as closed boundaries, the outlet is set to an open boundary. 

# Creates optional fields that will be used in the simulation
# fixed_nodes defines as 1 if a node is fixed or 0 if it can varies in elevation
fixed_nodes = np.zeros_like(z)  
fixed_nodes[fixed_nodes_id] = 1

rbd = RiverBedDynamics(grid,
                       gsd = gsd,
                       outlet_boundary_condition='fixedValue',
                       bed_surf__elev_fix_node=fixed_nodes)
       
""" Defines some variables to store data"""
save_data_now = True
plot_now = True                          # Used to save the plot at time zero
check_maximum_time = True
plot_time_interval_original = copy.deepcopy(plot_time_interval)         # A copy of plot_time_interval
save_data_time_interval_original=copy.deepcopy(save_data_time_interval) # A copy of save_data_time_interval
topographic__elevation_original = copy.deepcopy(grid["node"]["topographic__elevation"]) # A copy of the original topographic__elevation

progress0 = 0 # Initializes the variable
while t < sim_max_t:
    
    # defines the velocity at previous time    
    rbd._surface_water__velocity_prev_time_link = of._grid["link"]["surface_water__discharge"] / of._grid["link"]["surface_water__depth"]

    #Calculates the rainfall intensity - variable in time
    (rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value, save_data_now) = update_rainfall_intensity(t, rainfall_time, rainfall_time_index, rainfall_intensity, grid, current_rainfall_intensity, rainfall_intensity_value, save_data_now)
    of._rainfall_intensity = current_rainfall_intensity

    # Runs overland flow for one time step      
    of.overland_flow(dt=max_dt) 

    # defines the velocity at current time        
    grid["link"]["surface_water__velocity"] = grid["link"]["surface_water__discharge"] / grid["link"]["surface_water__depth"]
    
    # Defines the time step used in RiverBedDynamics
    rbd._grid._dt = of.dt
    
    # Runs RiverBedDynamics for one time step - bed evolution in this case
    rbd.run_one_step()  # Runs riverBedDynamics for one time step
      
    ## Stores results
    (save_data_time_interval, save_data_now) = save_data(save_data_time_interval, of, dt_precision, save_data_now, t, 
                                                        output_path, grid, link_list, node_list, save_data_time_interval_original, cwd, rainfall_intensity_value)
    # Plots results
    (plot_time_interval, plot_now) = plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval_original, cwd)

    # Updating t
    (t, plot_now, check_maximum_time, progress0) = update_time(t,of,sim_max_t,check_maximum_time, plot_now, progress0)