# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 13:59:52 2022

@author: angel
"""
import numpy as np

def saveRaster(data,filename,rmg):
    data = np.flip(np.reshape(data,[rmg.number_of_node_rows,rmg.number_of_node_columns]),axis=0)
    
    line0 = 'nrows ' + str(rmg.number_of_node_rows)
    line1 = 'cellsize '+ str(rmg.dx)
    line2 = 'xllcorner ' + str(rmg.origin[0])
    line3 = 'ncols ' + str(rmg.number_of_node_columns)
    line4 = 'yllcorner ' + str(rmg.origin[1])
    line5 = 'nodata_value -9999'
    with open(filename, 'w') as f:
        f.write(line0 + '\n'+ line1 + '\n' + line2 + '\n' + line3 + '\n' + line4 + '\n'+ line5 + '\n')
    with open(filename, "ab") as f:
        np.savetxt(f, data,'%.3f')