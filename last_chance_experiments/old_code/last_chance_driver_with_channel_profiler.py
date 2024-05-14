# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 11:31:13 2024

@author: Sam
"""

import numpy as np
import matplotlib.pyplot as plt
from landlab import RasterModelGrid, NodeStatus
from landlab.components import FlowAccumulator, SinkFillerBarnes
from landlab.io import read_esri_ascii
from numpy.ma import masked_array

# Function to find node ID based on cell value
def find_node_id(grid, value):
    indices = np.where(grid.at_node['topographic__elevation'] == value)
    return indices[0][0] if len(indices[0]) > 0 else None

# Load DEM with error handling
try:
    grid, z = read_esri_ascii(r'DEMs\filled_lc3_dem.txt', name='topographic__elevation')
except Exception as e:
    print(f"Error reading DEM file: {e}")
    raise

# Mask no-data values and set boundary conditionst
no_data_value = -9999
grid.set_nodata_nodes_to_closed(z, no_data_value)
cell_value_lc3 = 1381.671  # Replace with the elevation value for LC3
node_id_lc3 = find_node_id(grid, cell_value_lc3)
if node_id_lc3 is not None:
    grid.status_at_node[node_id_lc3] = NodeStatus.FIXED_VALUE

# Flow routing and sink filling
fa = FlowAccumulator(grid, flow_director='D4')
sf = SinkFillerBarnes(grid, method='Steepest')

# Run components with error handling
try:
    fa.run_one_step()
    sf.run_one_step()
except Exception as e:
    print(f"Error during flow routing or sink filling: {e}")
    raise

# Visualization
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.imshow(z.reshape(grid.shape), cmap='terrain')
plt.colorbar(label='Elevation (m)')
plt.title('DEM with Sink Filling')

plt.subplot(1, 2, 2)
fa_data = grid.at_node['drainage_area']
fa_data_masked = masked_array(fa_data, fa_data == no_data_value)
plt.imshow(np.log(fa_data_masked).reshape(grid.shape), cmap='cividis')
plt.colorbar(label='Log of Drainage Area (m^2)')
plt.title('Flow Accumulation')
plt.show()
