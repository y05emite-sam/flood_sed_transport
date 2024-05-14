# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:19:54 2023

@author: angel
"""
def update_time(t,of,sim_max_t,check_maximum_time, plot_now, progress0):
    
    if (t + of.dt > sim_max_t) and check_maximum_time:
        of.dt = sim_max_t - t
        t = sim_max_t
        plot_now = True
        check_maximum_time = False
    else:
        t += of.dt
        progress = int((t / sim_max_t) * 100)
        if progress > progress0 + 0.2:
            print("\r" + f"Progress: [{progress}%]", end="\n")
            progress0 = progress
    
    return t, plot_now, check_maximum_time, progress0

