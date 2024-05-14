# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 09:54:00 2024

@author: Sam
"""
"""
Driver file for: 
Case 3b_1 - Test in large watershed - Uniform rainfall - bed evolution

"""

""" Loads components and other libraries"""
#%reset -f

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
from landlab import imshow_grid


(output_folder, output_path, cwd) = clean_old_output()

""" Numerical simulation conditions and time control settings"""
# inputs filenames
zDEM = 'filled_lc1_dem.asc'                           # ASCII raster DEM containing the bed surface elevation
rainfall_time_series = 'rainfall_time_series.txt' # Rainfall time series in format #Time [s]	Precipitation [mm/hr]
gsd = np.loadtxt('LC1_grain_size_dist.txt')                             # Check inside the txt file for more information on the file format
bed_surf__gsd_loc_node = 'lc1_gsd_locations.asc'

# Time variables
t = 0
max_dt = 0.5                         # Overland flow will use the min time step between this value and the automatically calculated. Use seconds.

# Simulation conditions settings
sim_max_t = 1200+max_dt         # Maximum simulation time, the original was 24*3600+max_dt
save_data_time_interval = 10         # Stores results every this time
dt_precision = 1                    # Avoids rounding errors
plot_time_interval = 1800           # Plots will be obtained every this seconds

""" DEM Stuff """
# Load DEM and Set Up Grid 
(grid, z) = read_esri_ascii(zDEM, name='topographic__elevation')

# Needed for topographic variation plot
topographic__elevation_original = copy.deepcopy(grid.at_node['topographic__elevation'])

# Assuming -9999 is used for no-data
no_data_value = -9999

# Mask the no-data values in the DEM
masked_elevation = masked_array(z, z == no_data_value)
grid.at_node['topographic__elevation'] = masked_elevation.data
grid.add_zeros('surface_water__depth', at='node')  # Creates the surface_water__depth field



""" Depression Filling with SinkFiller (D4 routing) """
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

""" Visualization for Debugging: Flow accumulation """
flow_accumulation_data = grid.at_node['drainage_area']
masked_flow_accumulation = masked_array(flow_accumulation_data, flow_accumulation_data ==
no_data_value)
plt.imshow(masked_flow_accumulation.reshape(grid.shape), cmap='cividis', origin='lower')
plt.colorbar(label='Drainage Area (m^2)')
plt.title('Flow Accumulation with D4 Routing (No-Data Values Masked)')
plt.show()


# General simulation conditions
n = 0.025                                           # Manning's n
fixed_nodes_id = np.array((300,301,302,900,901,902,1500,1501,1502,2100,2101,2102,2700,2701,2702,3300,3301,3302))                    # Node ID for fixed Nodes - These will not change elevation
link_list = np.array([1500, 2099,1499, 900, 5097, 5696, 5096, 4497, 1291, 1890, 1220, 11691, 2085, 2384, 2384, 2245, 2866, 2865, 2865, 2866 ])    # Location of the links where information for postprocess will be extracted 
#The outlet node for LC3 is 33999
#The outlet node for LC1 is 35680
node_list = np.array([2701 ,6301, 33999])               # Location of the nodes where information for postprocess will be extracted 
rainfall_time_series = np.loadtxt(rainfall_time_series) # Loads the rainfall time series

""" OverlandFlow instantiation conditions """       
# Creates fields and instantiate the OverlandFlow component
OverlandFlow.input_var_names                                                # Gives the list of all required fields                                  
(grid, z) = read_esri_ascii(zDEM, name='topographic__elevation')             # Creates the topographic__elevation field
grid.add_zeros('surface_water__depth', at = 'node')                          # Creates the surface_water__depth field
of = OverlandFlow(grid, mannings_n=n, rainfall_intensity=0.0, alpha = 0.25, steep_slopes = True, theta = 1.0)    # instantiate the Overland flow component

(rainfall_time, rainfall_intensity, rainfall_time_index, current_rainfall_intensity, rainfall_intensity_value) = rainfall_data(rainfall_time_series,grid)
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

""" Put a thin layer of water """
if 'surface_water__depth' not in grid.at_node:
    grid.add_zeros('surface_water__depth', at='node')

# Add a thin layer of water (e.g., 1 mm)
grid.at_node['surface_water__depth'] += 0.001  # water depth in meters

# Initialize the OverlandFlow component
of = OverlandFlow(grid, mannings_n=0.025, rainfall_intensity=0.0, alpha=0.25, steep_slopes=False, theta=1.0)

# Check and adjust dt
if of.dt is None or np.isnan(of.dt):
    of.dt = 0.1  # Set an initial timestep if not automatically calculated
    
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
    
    #Bed surface variation plot
    if t % 10 == 0:  # Modulo operation to check if t is a multiple of 10
        # Plot and save the bed surface elevation variation
        ZVar = grid.at_node["topographic__elevation"] - topographic__elevation_original
        plot_name = f'Bed surface elevation variation [m] at {np.round(t, 0)} sec'
        imshow_grid(grid, ZVar, cmap='RdGy', vmin=-0.10, vmax=0.10, plot_name=plot_name)
        output = f"{output_folder}/topographicVariation_{np.round(t, 0)}.png"
        plt.savefig(output, dpi=300)
        plt.close()

    

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
    
    ORRRR
    
    Compare LC1 to LC3, find the reason for the inflection point
    
Here is under construction
    # Store and plot results
    (save_data_time_interval, save_data_now) = save_data(save_data_time_interval, of, dt_precision, save_data_now, t, 
                                                        output_path, grid, channel_links, channel_nodes, save_data_time_interval_original, cwd, rainfall_intensity_value)
    (plot_time_interval, plot_now) = plot_results(plot_time_interval, of, dt_precision, plot_now, output_folder, t, grid, topographic__elevation_original, plot_time_interval_original, cwd)

    # Update time
    (t, plot_now, check_maximum_time, progress0) = update_time(t, of, sim_max_t, check_maximum_time, plot_now, progress0)
    
    start with 2cm ish of water on the landscape 
    
    
   
    
"""
    
    