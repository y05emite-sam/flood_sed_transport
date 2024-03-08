# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 13:49:00 2024

@author: Sam
"""
from landlab.io import write_esri_ascii
import numpy as np
import matplotlib.pyplot as plt
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, SinkFillerBarnes
from landlab.io import read_esri_ascii
from numpy.ma import masked_array
from landlab import NodeStatus  # Add this import statement

# Load DEM
(grid, z) = read_esri_ascii('lc3_dem.txt', name='topographic__elevation')

# Define no-data value and mask it
no_data_value = -9999


# Find the outlet node for LC1 and LC3, LC1 is elevation 1517.548, in LC3 its 1381.671

def find_node_id(grid, value):
    """
    Find the node ID in a Landlab grid for a given cell value.

    Parameters:
    grid (RasterModelGrid): The Landlab grid.
    value (float): The value to find.

    Returns:
    int: Node ID of the cell, or None if not found.
    """
    # Find the indices of the cell with the given value
    indices = np.where(grid.at_node['topographic__elevation'] == value)

    if len(indices[0]) > 0:
        # Return the first occurrence if multiple found
        return indices[0][0]
    else:
        return None


# Replace this with the value you are looking for
cell_value = 1381.671
# 1381.671 is the value for LC3
# 1517.548 is the value for LC1

node_id = find_node_id(grid, cell_value)
if node_id is not None:
    print(f"Cell with value {cell_value} found at node ID: {node_id}")
else:
    print(f"Cell with value {cell_value} not found.")
    
grid.set_nodata_nodes_to_closed(z, -9999)

#Here change the node ID 33999 for LC3
#Here change the node ID 35680 for LC1

grid.status_at_node[33999] = NodeStatus.FIXED_VALUE


#grid.set_watershed_boundary_condition(grid, -9999)



# Flow Routing with D4 Routing
#fd = FlowDirectorSteepest(grid)
fa = FlowAccumulator(grid, flow_director='D4') # tell it im doing d4
fa.run_one_step()
grid.at_node["flow__sink_flag"][grid.core_nodes].sum()
print(grid.at_node["flow__sink_flag"][grid.core_nodes].sum())

      
# Sink Filling for D4 Flow Routing
sf = SinkFillerBarnes(grid, method='Steepest')
sf.run_one_step()
print(grid.at_node["flow__sink_flag"][grid.core_nodes].sum())

# Flow Routing with D4 Routing
#fd = FlowDirectorSteepest(grid)
fa = FlowAccumulator(grid, flow_director='D4') # tell it im doing d4
fa.run_one_step()
grid.at_node["flow__sink_flag"][grid.core_nodes].sum()
print(grid.at_node["flow__sink_flag"][grid.core_nodes].sum())

# Visualization
fig, axs = plt.subplots(2, 2, figsize=(12, 12))

# Define a function to ensure a minimum range for color scale
def ensure_min_range(data, min_range=10):
    if np.ma.is_masked(data) and data.mask.all():
        # If all data is masked, return a default range
        return 0, min_range
    else:
        min_val, max_val = np.min(data), np.max(data)
        if max_val - min_val < min_range:
            mid_val = (max_val + min_val) / 2
            return mid_val - min_range / 2, mid_val + min_range / 2
        else:
            return min_val, max_val

# Full DEM visualization
full_dem_range = ensure_min_range(z)
full_dem_img = axs[0, 0].imshow(z.reshape(grid.shape), cmap='terrain', vmin=full_dem_range[0], vmax=full_dem_range[1], origin='lower')
axs[0, 0].set_title('Full DEM')
plt.colorbar(full_dem_img, ax=axs[0, 0], orientation='vertical', fraction=0.046, pad=0.04)

# Full Flow Accumulation visualization
fa_data = grid.at_node['drainage_area']
fa_data_masked = np.log(masked_array(fa_data, fa_data == no_data_value))
full_fa_range = ensure_min_range(fa_data_masked)
full_fa_img = axs[0, 1].imshow(fa_data_masked.reshape(grid.shape), cmap='cividis', vmin=full_fa_range[0], vmax=full_fa_range[1], origin='lower')
axs[0, 1].set_title('Full Flow Accumulation')
plt.colorbar(full_fa_img, ax=axs[0, 1], orientation='vertical', fraction=0.046, pad=0.04)

#Save filled DEM
output_filled_dem_path = 'filled_lc3_dem.txt'
write_esri_ascii(output_filled_dem_path, grid, names='topographic__elevation')
