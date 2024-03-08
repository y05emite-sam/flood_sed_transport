# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 18:32:16 2024

@author: Sam
"""
import numpy as np
import matplotlib.pyplot as plt
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FlowDirectorSteepest, SinkFillerBarnes
from landlab.io import read_esri_ascii
from numpy.ma import masked_array

# Load DEM
(grid, z) = read_esri_ascii('lc3_dem.txt', name='topographic__elevation')

# Define no-data value and mask it
no_data_value = -9999
z_masked = masked_array(z, z == no_data_value)
#grid.at_node['topographic__elevation'] = z_masked.data
grid.set_nodata_nodes_to_closed(z_masked, -9999)

#grid.set_watershed_boundary_condition(grid, -9999)

"""
if set wtr shd boundary isnt workin set all -9999 to closed. manually figure out elevaiton outlet, and set that node to have an open boundary condition
"""


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

# Convert real-world coordinates to grid indices
def convert_coords_to_indices(x, y, grid):
    spacing = grid.dx  # Assuming square cells and uniform spacing
    x_index = int((x - grid.node_x[0]) / spacing)
    y_index = int((y - grid.node_y[0]) / spacing)
    return x_index, y_index

x_coord, y_coord = 528000, 3.56825
center_col, center_row = convert_coords_to_indices(x_coord, y_coord, grid)
zoom_size = 15  # Half-size for 30x30 cells

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
full_dem_range = ensure_min_range(z_masked)
full_dem_img = axs[0, 0].imshow(z_masked.reshape(grid.shape), cmap='terrain', vmin=full_dem_range[0], vmax=full_dem_range[1], origin='lower')
axs[0, 0].set_title('Full DEM')
axs[0, 0].add_patch(plt.Rectangle((center_col - zoom_size, center_row - zoom_size), 2 * zoom_size, 2 * zoom_size, fill=False, edgecolor='red', linewidth=2))
plt.colorbar(full_dem_img, ax=axs[0, 0], orientation='vertical', fraction=0.046, pad=0.04)

# Full Flow Accumulation visualization
fa_data = grid.at_node['drainage_area']
fa_data_masked = np.log(masked_array(fa_data, fa_data == no_data_value))
full_fa_range = ensure_min_range(fa_data_masked)
full_fa_img = axs[0, 1].imshow(fa_data_masked.reshape(grid.shape), cmap='cividis', vmin=full_fa_range[0], vmax=full_fa_range[1], origin='lower')
axs[0, 1].set_title('Full Flow Accumulation')
plt.colorbar(full_fa_img, ax=axs[0, 1], orientation='vertical', fraction=0.046, pad=0.04)

"""
# Zoomed-in area calculation
def zoom_in_area(center_row, center_col, zoom_size, data):
    start_row = max(center_row - zoom_size, 0)
    end_row = min(center_row + zoom_size, grid.number_of_node_rows)
    start_col = max(center_col - zoom_size, 0)
    end_col = min(center_col + zoom_size, grid.number_of_node_columns)
    return data.reshape(grid.shape)[start_row:end_row, start_col:end_col]

# Zoomed-in DEM visualization
zoomed_dem = zoom_in_area(center_row, center_col, zoom_size, z_masked)
zoomed_dem_range = ensure_min_range(zoomed_dem)
zoomed_dem_img = axs[1, 0].imshow(zoomed_dem, cmap='terrain', vmin=zoomed_dem_range[0], vmax=zoomed_dem_range[1], origin='lower')
axs[1, 0].set_title('Zoomed-in DEM')
plt.colorbar(zoomed_dem_img, ax=axs[1, 0], orientation='vertical', fraction=0.046, pad=0.04)

# Zoomed-in Flow Accumulation visualization
zoomed_fa = zoom_in_area(center_row, center_col, zoom_size, fa_data_masked)
zoomed_fa_range = ensure_min_range(zoomed_fa)
zoomed_fa_img = axs[1, 1].imshow(zoomed_fa, cmap='cividis', vmin=zoomed_fa_range[0], vmax=zoomed_fa_range[1], origin='lower')
axs[1, 1].set_title('Zoomed-in Flow Accumulation')
plt.colorbar(zoomed_fa_img, ax=axs[1, 1], orientation='vertical', fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()
"""
