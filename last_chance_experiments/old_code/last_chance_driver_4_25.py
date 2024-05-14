# -*- coding: utf-8 -*-
"""
@author: Sam
"""

""" Loads components and other libraries"""
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from landlab.components import (FlowAccumulator, FlowDirectorSteepest,
                                SinkFiller, ChannelProfiler, 
                                OverlandFlow, RiverBedDynamics)
from landlab.io import read_esri_ascii
from clean_old_output import clean_old_output
from save_data_to_file import save_data
from save_raster import plot_results
from rainfall_manager import rainfall_data, update_rainfall_intensity
from update_time import update_time
from numpy.ma import masked_array

# Clean old output
(output_folder, output_path, cwd) = clean_old_output()

""" DEM and GSD"""
# Inputs filenames
zDEM = r'DEMs\filled_lc1_dem.txt'  # ASCII raster DEM containing the bed surface elevation
rainfall_time_series = pd.read_excel('C:/Users/Sam/MultiGrain/last_chance_experiments/climate/5min_1yrRI_storm.xlsx')
gsd = pd.read_excel('GSDs/LC1_grain_size_dist.xlsx', sheet_name='GSD', skiprows=0).values

""" Rain stuff"""
# Extract time series and precipitation from rainfall data
time_series = rainfall_time_series['Time [s]'].values
precipitation = rainfall_time_series['Precipitation [mm/hr]'].values

# Convert DataFrame to NumPy array
rainfall_time_series_array = rainfall_time_series[['Time [s]', 'Precipitation [mm/hr]']].to_numpy()

""" Assign variables"""
# Time variables and simulation conditions
t = 0
max_dt = 0.01  # Use seconds
sim_max_t = 24*3600 + max_dt
save_data_time_interval = 5
dt_precision = 3
plot_time_interval = 1800

# General simulation conditions
n = 0.025  # Manning's n

"""DEM Stuff"""
# Load DEM and Set Up Grid 
(grid, z) = read_esri_ascii(zDEM, name='topographic__elevation')

# Assuming -9999 is used for no-data
no_data_value = -9999

# Mask the no-data values in the DEM
masked_elevation = masked_array(z, z == no_data_value)
grid.at_node['topographic__elevation'] = masked_elevation.data
grid.add_zeros('surface_water__depth', at='node')  # Creates the surface_water__depth field

# Depression Filling with SinkFiller (D4 routing) 
# Sinks were already filled in preprocessing, so this step may not be necessary
sf = SinkFiller(grid, routing='D4')  # Use D4 routing
sf.fill_pits()

# Visualization of Filled DEM
filled_z = grid.at_node['topographic__elevation']
masked_filled_z = masked_array(filled_z, filled_z == no_data_value)
plt.imshow(masked_filled_z.reshape(grid.shape), cmap='cividis', origin='lower')
plt.colorbar(label='Elevation (m)')
plt.title('Filled DEM with SinkFiller (No-Data Values Masked)')
plt.show()

""" Flow Routing with D4 Routing """
fd = FlowDirectorSteepest(grid)  # Specify D4 routing
fa = FlowAccumulator(grid, flow_director=fd)
fa.run_one_step()

# Visualization for Debugging: Flow accumulation
flow_accumulation_data = grid.at_node['drainage_area']
masked_flow_accumulation = masked_array(flow_accumulation_data, flow_accumulation_data ==
no_data_value)
plt.imshow(masked_flow_accumulation.reshape(grid.shape), cmap='cividis', origin='lower')
plt.colorbar(label='Drainage Area (m^2)')
plt.title('Flow Accumulation with D4 Routing (No-Data Values Masked)')
plt.show()

#Here change the node ID 33999 for LC3
#Here change the node ID 35680 for LC1

""" Channel Profiling """
dxy = 10  # side length of a raster model cell, or resolution [m], to make cp run

# Initialize ChannelProfiler REMEBER TO USE FOR D4
cp = ChannelProfiler(grid, 
                     #outlet_nodes=[35680], 
                     number_of_watersheds=1, 
                     main_channel_only=True,
                     minimum_channel_threshold=dxy**2)
cp.run_one_step()

# Check if ChannelProfiler has data, and plot to debug
if cp.data_structure:
    for watershed_id, segments in cp.data_structure.items():
        for segment_id, segment_data in segments.items():
            # Get the node IDs for the current segment
            node_ids = segment_data["ids"]
            
            # Extract the elevations and distances for these nodes
            elevations = grid.at_node['topographic__elevation'][node_ids]
            distances = segment_data["distances"]
            
            # Check if elevations and distances are valid
            if elevations.size > 0 and distances.size > 0:
                # Plot the profile
                plt.plot(distances, elevations, label=f'Segment {segment_id}')
            else:
                print(f"No valid elevation or distance data for segment {segment_id} in watershed {watershed_id}")

    plt.xlabel('Distance Upstream (m)')
    plt.ylabel('Elevation (m)')
    plt.title('Channel Profiles')
    plt.legend()
    plt.show()
else:
    print("No channel data found in ChannelProfiler.")

""" Nodes and links for data viz"""
# Links and nodes to record info at, NEED TO BE UODATED- THESE ARE ANGELS
link_list = np.array([1500, 2099,1499, 900, 5097, 5696, 5096, 4497, 12291, 12890, 12290 ])    # Location of the links where information for postprocess will be extracted 
node_list = np.array([2701 ,6301])               # Location of the nodes where information for postprocess will be extracted 

""" Instantiate OverlandFlow and RiverBedDynamics """
# Initialize the surface_water__velocity field with zeros as a placeholder
# This should be done before RiverBedDynamics is instantiated
grid.add_zeros('surface_water__velocity', at='link')

of = OverlandFlow(grid, mannings_n=n, rainfall_intensity=0.0, alpha=0.25, steep_slopes=False, theta=1.0)

# Call rainfall_data with the NumPy array
(rainfall_time, rainfall_intensity, rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value) = rainfall_data(rainfall_time_series, grid)
of._rainfall_intensity = current_rainfall_intensity

rbd = RiverBedDynamics(grid, gsd=gsd, outlet_boundary_condition='fixedValue', bed_surf__elev_fix_node=np.zeros_like(z))


""" Simulation Loop """
rainfall_time_index = 0  # Start at the beginning of the rainfall time series
current_rainfall_intensity = 0.0  # Initial rainfall intensity, often starts at zero
rainfall_intensity_value = 0.0    # Initial value of rainfall intensity, can be zero or first value from your data

""" Defines some variables to store data"""
save_data_now = True
plot_now = True                          # Used to save the plot at time zero
check_maximum_time = True
plot_time_interval_original = copy.deepcopy(plot_time_interval)         # A copy of plot_time_interval
save_data_time_interval_original=copy.deepcopy(save_data_time_interval) # A copy of save_data_time_interval
topographic__elevation_original = copy.deepcopy(grid["node"]["topographic__elevation"]) # A copy of the original topographic__elevation
progress0 = 0 # Initializes the variable

""" Set initial dt"""
# Check and adjust dt
if of.dt is None or np.isnan(of.dt):
    of.dt = 0.1  # Set an initial timestep if not automatically calculated
    
""" Simulation Loop """
while t < sim_max_t:
    # Debugging: Print current time step
    print(f"Current t: {t}, Current of.dt: {of.dt}")

    # Check if of.dt is None, NaN, or a negative number
    if of.dt is None:
        print("Warning: of.dt is None")
        break  # Exit the loop or handle as needed
    elif np.isnan(of.dt) or of.dt < 0:
        print("Warning: Invalid of.dt detected (of.dt: {of.dt}) at time t: {t}")
        break  # Break out of the loop or handle the situation as needed

    # Defines the velocity at previous time    
    grid.at_link['surface_water__velocity'] = grid.at_link['surface_water__discharge'] / grid.at_link['surface_water__depth']

    # Calculates the rainfall intensity - variable in time
    (rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value, save_data_now) = update_rainfall_intensity(t, rainfall_time, rainfall_time_index, rainfall_intensity, grid, current_rainfall_intensity, rainfall_intensity_value, save_data_now)
    of._rainfall_intensity = current_rainfall_intensity

    # Runs overland flow for one time step      
    of.overland_flow(dt=max_dt) 

    # Defines the velocity at current time        
    grid.at_link['surface_water__velocity'] = grid.at_link['surface_water__discharge'] / grid.at_link['surface_water__depth']
    print(f"Current t: {t}, Current of.dt: {of.dt}")

    # Defines the time step used in RiverBedDynamics
    rbd._grid._dt = of.dt
    
    # Runs RiverBedDynamics for one time step - bed evolution in this case
    rbd.run_one_step()

    # This is to check if the code made it this far
    print(f"Current t: {t}, Current of.dt: {of.dt}")

    # Stores results
    (save_data_time_interval, save_data_now) = save_data(save_data_time_interval, of, dt_precision, save_data_now, t, output_path, grid, link_list, node_list, save_data_time_interval_original, cwd, rainfall_intensity_value)

    # Plots results
    (plot_time_interval, plot_now) = plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval_original, cwd)

    # Updating t and other variables using your update_time function
    (t, plot_now, check_maximum_time, progress0) = update_time(t, of, sim_max_t, check_maximum_time, plot_now, progress0)
    # Check if t is NaN after update_time
    if np.isnan(t):
        print(f"Warning: t became NaN after update_time at time step: {of.dt}")
        break  # Break out of the loop or handle the situation as needed

    
"""
i want to brainstorm a figure to show residence time in the channel. shall i do channel steepness vs RI? 

# 4 figures total
# Figure 1. Map (with ksn), sediment size distribution, boulder size distribution,
    bed thickness distributions, maybe chi plots
    Basically this will introduce area with any relevant things (what is relevant?)
    channel steepness, sediment and boulder size, bed thickness

# Figure 2. I need some way of showing the effect of bed thicknesses on sediment and boulder size,
    and some way of showing the effect of sediment and boulder size on channel steepnesses.
    This is a range of values and so maybe can plot all the data points instead of averages

# Figure 3. Modeling, how often do storms move sediment? move boulders? 
    Plot that relates RI of "geomorphic work" to steepness
    RI on xaxis, size of sediment moved
    will any of the big things move in steep area? how much sediment motion there is?
    maybe a map which shows where bedrock is eroded and where it's not for different RI's
    put in largest 50 grain sizes for modeling
    for this storm the largest sediment that moved is this
    start with largest storm and see if sediment is mobilized 
    if sediment is moved period
    can i query to see if sediment moved in a channel section
    
# Figure 4. How could climate change "erase" this top down signal? Where does this negative feedback loop
cease to matter? is there some upperlimit to boulder size where it can resist higher intensity storms? 
    map of largest boulder size mobilized per storm? 
    
Here is under construction
    # Store and plot results
    (save_data_time_interval, save_data_now) = save_data(save_data_time_interval, of, dt_precision, save_data_now, t, 
                                                        output_path, grid, channel_links, channel_nodes, save_data_time_interval_original, cwd, rainfall_intensity_value)
    (plot_time_interval, plot_now) = plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval_original, cwd)

    # Update time
    (t, plot_now, check_maximum_time, progress0) = update_time(t, of, sim_max_t, check_maximum_time, plot_now, progress0)
    
    start with 2cm ish of water on the landscape 
    
    
   
    
"""

