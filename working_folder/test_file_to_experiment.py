# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 18:22:27 2022
@author: sam, angel, and mikey

experiment for lc3, later lc1 will be added, and then will become a jupyterhub
notebook tutorial for the new river bed dynamics component

"""

""" first we imporant all necessary components and other stuff"""

#%reset -f
import numpy as np
import pandas as pd
import copy
import os
import shutil
from matplotlib import pyplot as plt
from landlab.components import OverlandFlowSpatiallyVariableInputs, RiverBedDynamics
from landlab.io import read_esri_ascii
from landlab import imshow_grid

""" then we add the dem, grain size distribution info, and precipitation data. 
In the end there will be a folder where all the different precip data lives, 
but I havent gotten around it it yet. currently trying to get the model to run
and show different grain sizes in the different locations. NOTE: column 1 in the
gsd file is data for the 'shallow' channel section (above 1560m elevation), 2 
is 1530m - 1560m elevation, and 3 is steep- elevation lower than 1530. 

The two DEM's are called:
    lc3_dem.txt
    and lc1_dem.txt
    
The corresponding GSD info for LC1 is...
    LC1_grain_size_dist.xlsx
    and LC1_gsd_locations.txt
    
...and for LC3 is 
    LC3_grain_size_dist.xlsx
    and LC3_gsd_locations.txt    
"""

#This imports the DEM
watershed_dem = 'lc3_dem.txt' 
(rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')

# This stuff is precip data
rainfallFile = '5min_1000yrRI_storm_upper_bound.xlsx'
precipitation = pd.read_excel(rainfallFile)

# These import the three sediment size distributions
gsd = pd.read_excel('LC3_grain_size_dist.xlsx',sheet_name='GSD',skiprows=0).values

# And this specifies the lcation of the three GSD's
bedGSDLocationRaster = 'LC3_gsd_locations.txt'     

(rmg0, gsd_loc) = read_esri_ascii(bedGSDLocationRaster)

rmg['node']['bed_surface__grainSizeDistribution_location'] = gsd_loc   


""" this is info regarding time, and plotting stuff. remember some of the storms
will be longer than 3600 seconds, and so that will change depending on the storm
length.
"""
dtPrecision = 3             # Avoids rounding errors
max_dt = 1                  # Overland flow will use the min time step between this value and the automatically calculated. Use seconds.
tPlot = 600                  # Plots will be obtained every this seconds
storeData = 10              # Stores results every this time
tmax = 6000 + max_dt          # Maximum simulation time, adding max_dt ensures that the last time is stored

"""mmanning's n value (roughness), prolly wont be changed"""
n = 0.03                            

"""Link and node at base of DEM, where samples will be collected, i think mikey 
changed these numbers?"""
link_to_sample = 698
node_to_sample = 300

"""This removes previous figs, im not sure I'll keep it because I think it's 
deleting all the txt files (the DEMS included), and it's getting annoying 
putting them back in the working folder over and over again. UODATE: I commented
out the text file removal stuff"""
directory = os.getcwd() ; test = os.listdir( directory )

for item in test:
    if item.endswith(".png"):
        os.remove( os.path.join( directory, item ) )
   # if item.endswith(".txt"):
    #    os.remove( os.path.join( directory, item ) )  
        
"""Creates fields and instantiate the component"""
OverlandFlowSpatiallyVariableInputs.input_var_names
RiverBedDynamics.input_var_names
# (rmg, z) = read_esri_ascii(bedElevation, name='topographic__elevation')
rmg.add_zeros('bed_surface__roughness', at = 'link')
rmg.add_zeros('surface_water__depth', at = 'node')
rmg.add_zeros('rainfall__intensity', at = 'node')
#rmg['node']['bed_surface__grainSizeDistribution_location'] = np.zeros_like(z)     

"""this generates a slope map, I put it in to check and see that the DEM was 
imported correctly"""
rmg.at_node['topographic__slope'] = rmg.calc_slope_at_node(elevs='topographic__elevation')
imshow_grid(rmg,'topographic__slope');
plt.show()


""" this instatiates the two components that are needed to run the model,
i didnt comment out the two lines below, maybe mikey did?"""
of = OverlandFlowSpatiallyVariableInputs(rmg, steep_slopes=True, alpha = 0.3)
rbd = RiverBedDynamics(rmg , gsd = gsd, variableCriticalShearStress = True, bedloadEq='WilcockAndCrowe')

#z1 = z.reshape(382,469)
#print(np.where(z1==1382.996))

""" Set boundaries as closed boundaries, the outlet is set to an open boundary. 
theres some stuff here that is commented out, maybe mikey commited them out when
he set the boundary conditions?"""

#rmg.set_closed_boundaries_at_grid_edges(False, True, True, True)
rmg.set_watershed_boundary_condition_outlet_id([33999], z, nodata_value=-9999.) #[col,row] = [382,469] : 166074
outlet = rmg.set_watershed_boundary_condition(z, remove_disconnected=True, nodata_value=-9999., return_outlet_id=True) #1382.996
print(outlet)

""" Create bed and flow initial condition, remember that n is mannings roughness coefficient"""
rmg['link']['bed_surface__roughness'] = np.zeros(rmg.number_of_links) + n  

""" this is precip stuff, kinda confused about unit conversions, make sure to
go back and double check"""
precipitation = pd.read_excel(rainfallFile)
precip_time=precipitation.values[:,0]
precip_mmhr=precipitation.values[:,1]
precip_ms = precip_mmhr * (2.77778 * 10 ** -7)  # Converts mm/hr to m/s
precip_index = 0                                # current index for time

""" Defines variables to store data and run the experiment """
storeNow = True
plotNow = True                          # Used to save the plot at time zero
check_tmax = True
tPlotOrg=copy.deepcopy(tPlot)           # A copy of tPlot, for plotting purposes (i think)
storeDataOrg=copy.deepcopy(storeData)   # A copy of tPlot
outputFolder = 'output'
cwd = os.getcwd()

if os.path.exists(outputFolder):
    print('The folder') 
    print(outputFolder)
    print('Exists and it will be removed \n');
    shutil.rmtree(outputFolder)     
os.mkdir(outputFolder)

"""  runs the experiment """
t = 0                                   # Initializates the variable
while t < tmax:
    
    rbd.t = t           # Current simulation time
    
    #Calculates the rainfall intensity - variable in time
    if (t >= precip_time[precip_index+1]):
        rmg['node']['rainfall__intensity'] =  np.zeros(rmg.number_of_nodes) + precip_ms[precip_index+1]
        precip_index += 1
    else:
        rmg['node']['rainfall__intensity'] = np.zeros(rmg.number_of_nodes) + precip_ms[precip_index+1]
    
    of.overland_flow() # Runs overland flow for one time step

    rbd.run_one_step()  # Runs riverBedDynamics for one time step
    
    ## Stores results
    storeData = round(storeData-of.dt, dtPrecision)
    if (storeData <=0) or storeNow:
        os.chdir(outputFolder)
        print('Storing results at time :',np.round(t,1),' s \n')
        data = np.hstack([t,(np.abs(of._q[link_to_sample] * rmg.dx).T)])
        with open("output0_links_surface_water__discharge.txt", "ab") as f:
            np.savetxt(f, data,'%.3f')
        data = np.hstack([t,(of._h[node_to_sample].T)])
        with open("output1_node_surface_water__depth.txt", "ab") as f:
            np.savetxt(f, data,'%.3f')      
        data = np.hstack([t,np.abs(rbd._tau[link_to_sample].T)])
        with open("output2_link_surface_water__shearStress.txt", "ab") as f:
            np.savetxt(f, data,'%.3f')   
        data = np.hstack([t,rmg.at_node["topographic__elevation"][node_to_sample].T])
        with open("output3_node_topographic__elevation.txt", "ab") as f:
            np.savetxt(f, data,'%.3f') 
        data = np.hstack([t,rmg.at_link["bed_surface__medianSize"][link_to_sample].T])
        with open("output4_link_bed_surface__medianSize.txt", "ab") as f:
            np.savetxt(f, data,'%.3f')
        data = np.hstack([t,rmg.at_link['sediment_transport__bedloadRate'][link_to_sample].T])
        with open("output5_links_sediment_transport__bedloadRate.txt", "ab") as f:
            np.savetxt(f, data,'%.5f')  
        storeData = round(storeDataOrg, dtPrecision)
        storeNow = False
        os.chdir(cwd)

    tPlot = round(tPlot-of.dt, dtPrecision)
    if tPlot <= 0  or plotNow:
        os.chdir(outputFolder)
        print('Elapsed time :',np.round(t,1),' s. Current dt =',\
              np.round(of.dt,1),'. Adaptive time =',np.round(of._adaptive_dt,1),' s - Saving plot \n')
        
        # Water depth plot
        plot_name='Surface water depth [m] at ' + str(np.round(t,0)) + ' sec'
        imshow_grid(rmg, 'surface_water__depth',cmap='Blues',vmin=0,plot_name=plot_name)
        output='depth_'+str(np.round(t,0))+'.png'
        plt.savefig(output,dpi=300); plt.close()  
        
        #Bed surface elevation plot
        plot_name='Bed surface elevation [m] at ' + str(np.round(t,0)) + ' sec'
        ZBed = rmg.at_node["topographic__elevation"]
        imshow_grid(rmg, ZBed ,cmap='RdGy',plot_name=plot_name)
        output='topographicElevation_'+str(np.round(t,0))+'.png'
        plt.savefig(output,dpi=300); plt.close()  
        
        #Bed surface variation plot
        plot_name='Bed surface elevation variation [m] at ' + str(np.round(t,0)) + ' sec'
        ZVar = rmg.at_node["topographic__elevation"] - rmg.at_node['topographic__elevation_original'] 
        imshow_grid(rmg, ZVar,cmap='RdGy',plot_name=plot_name)
        output='topographicVariation_'+str(np.round(t,0))+'.png'
        plt.savefig(output,dpi=300); plt.close()    

        plotNow = False
        tPlot = tPlotOrg
        os.chdir(cwd)

    # Updating t
    if (t + of.dt > tmax) and check_tmax:
        of.dt = tmax - t
        t = tmax
        storeDataNow = True  
        plotNow = True
        check_tmax = False
    else:
        t = round(t + of.dt, dtPrecision)  
        
"""we need to make some figs... 
1) time and storm size it takes to move minimum, avg, max sediment in different channel sections
2) percent of material removed from different channel sections for different storms
3) material size removed (min, max, avg) for different different storms in different slope, channel steepness, and drainage area 
4) material deposited (min, max, avg) for different slope, chi, and drainage area 

    
    BUT first plot info about hydrographs and shear stresses produced for different storms
"""
