from cmath import inf
from os import name
import matplotlib.pyplot as plt
from matplotlib.pyplot import legend
import numpy as np
import math
from scipy.signal import find_peaks
import re

def amplitude_mod_axis(A,B):
    #Calculate the amplitude modulation along a specific axis
    result=np.abs(np.abs(A+B)-np.abs(A-B))

    return result

def activating_function(field, field_grid, axis):
    length_field=len(field)
    result=[]
    for i in range(1,length_field):
        AF=0
        delta_y=(field[i,axis].real-field[i-1,axis].real)
        delta_x=(field_grid[axis,i].real-field_grid[axis,i-1].real)
        AF=delta_y/delta_x
        result.append(AF)
    return result

def activating_function_3D(field, field_grid):

    length_field=len(field)
    result=[]
    for i in range(1,length_field):
        AF=0
        for n in range(3):
            delta_y=(field[i,n].real-field[i-1,n].real)
            delta_x=(field_grid[n,i].real-field_grid[n,i-1].real)
            if abs(delta_x)>0:
                ratio=delta_y/delta_x
                AF=AF+(ratio)**2 
        result.append(math.sqrt(AF))
    return result

def loc_to_index(location,lower_limit, upper_limit, number_of_points):
    result=int(((location-lower_limit)/(upper_limit-lower_limit))*number_of_points)
    return result

def count_spikes(V_plot, threshold,beat,beat_range):
    
    #Skip the first half beat 
    time_step=V_plot[0,1]-V_plot[0,0] #remember that this is in milliseconds
    beat_steps=int((1/beat)*1000/time_step)

    check_values=V_plot[1,beat_steps*beat_range[0]:beat_steps*beat_range[1]]

    peaks, _ = find_peaks(check_values, height=threshold, distance=100)
    
    return len(peaks+1) # add 1 to count 0

def spikes_category(spikes):
    #returns the category of type of activation based on number of spikes/ wavelength at the nodes along the axon

    #spikes per envelope wavelength- split into categories: 
        # no spikes (val=0)- colour =white
        # 1 spike (val=1)- colour= green
        # 2-3 spikes (val=2)- colour = amber
        # 4+ spikes (val=4) - colour = red

    minimum_spikes=np.amin(spikes)

    if minimum_spikes>3:
        result=4
    elif minimum_spikes>1:
        result=2
    elif minimum_spikes==1:
        result=1
    else:
        result=0
        
    return result

def Get_max_abs(field):
    max_abs=0
    for i in range(len(field)):
        abs_val=math.sqrt(field[i,0].real**2+field[i,1].real**2+field[i,2].real**2)
        if abs_val>max_abs:
            max_abs=abs_val
    return max_abs

def Current_density(field):
    max_abs=0
    sum=0
    av=0
    for i in range(len(field)):
        abs_val=math.sqrt(field[i,0].real**2+field[i,1].real**2+field[i,2].real**2)
        sum=sum+abs_val
        if abs_val>max_abs:
                max_abs=abs_val
    av=sum/len(field)
    return max_abs,sum,av

def Av_interspike_interval(results,threshold, time_window):
    #Check the timestep 
    time_step=(results[0,1]-results[0,0])*1000
    window_samples=round(time_window/time_step)
    
    #Detect the spikes
    peaks, _ = find_peaks(results[1,:], height=threshold, distance=window_samples)

    
    if len(peaks)<2:
        #if less than 2 spikes detected return inf as interspike interval
        interspike=inf
    else:
        #Get average timing
        sum_intervals=0
        for i in range(len(peaks)-1):
            prev_index=peaks[i]
            index=peaks[i+1]
            sum_intervals=sum_intervals+results[0,index]-results[0,prev_index]
        
        #Get average
        interspike=sum_intervals*1000/(len(peaks)-1)
    number_of_spikes=len(peaks)
    return(interspike, number_of_spikes)

def spike_per_window_check(results,threshold, window,create_plot=False):

    #Check the timestep 
    time_step=(results[0,1]-results[0,0])*1000
    window_samples=round(window/time_step)
    n_windows=math.floor(len(results[0,:])/window_samples)-1

    if create_plot==True:
        fig, axs = plt.subplots(n_windows+1,1,figsize=(6*n_windows,12))
    output=1
    n=0

    #loop through windows 
    while (output==1 and n<=n_windows):
        #Detect the spikes
        start=int(n*window_samples)
        end=int((n+1)*window_samples)
        if create_plot==True:
            axs[n].plot(results[0,start:end],results[1,start:end])
        peaks, _ = find_peaks(results[1,start:end], height=threshold)
        print(len(peaks))
        if len(peaks)==0:
            output=0
        n+=1
            
    return output

def get_node(section_name):
    node = re.search(r"node\[(\w+)\]", section_name)
    node=int(node.group(1))
    return node