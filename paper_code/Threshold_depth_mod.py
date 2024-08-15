from __future__ import absolute_import
from __future__ import print_function
from statistics import median

from sklearn.neighbors import radius_neighbors_graph

#Import s4l modules
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.simulation.emlf as emlf
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewers
from s4l_v1.model import Vec3
import s4l_v1.units as units
from s4l_v1 import Unit
import s4l_v1.simulation.neuron as neuron

#import other modules
import sys, os
import numpy
from datetime import datetime
import itertools
import h5py 
import random
import math
import csv

#Add main file to path
sys.path.insert(1,r'C:\Users\lj386\OneDrive - University of Bath\PhD\Modelling\Sim4life\code\sim4life')

#import commonly used sim4life functions 
import s4L_functions as fun
import processing_functions as process 
import modelling_functions as modelling

#Create a folder to store the results 

#Define where the results will be saved 
parent_dir = r'C:\Users\lj386\OneDrive - University of Bath\PhD\Modelling\IFS paper\Results'

#Create a folder for today's results as a folder with the date and time
now = datetime.now()
main_folder = now.strftime("%Y_%m_%d_%H_%M")
path =os.path.join(parent_dir, main_folder)
print(path)
os.mkdir(path)

#Import arm model 
layer_names=['Skin','SAT','Muscle','Radius','Ulna','Ulnar Nerve','Median Nerve','Superficial Radial Nerve']
layer_materials=["Skin","SAT (Subcutaneous Fat)","Muscle", "Bone (Cortical)", "Bone (Cortical)",'Nerve' ,'Nerve','Nerve']

#Import Electrodes 
Electrodes=[]
Gels=[]

for i in range(2):
    Electrodes.append(model.AllEntities()['Electrode '+str(i+1)])
    Gels.append(model.AllEntities()['Gel '+str(i+1)])
pos_elec=Electrodes[0]
neg_elec=Electrodes[1]
Gel_resistivity=0.01


#Import details of median nerve axons 
model_dir = r'C:\Users\lj386\OneDrive - University of Bath\PhD\Modelling\IFS paper\Models'
# file_name='P4_median.csv'
# model_path =os.path.join(model_dir,file_name)
# file = open(model_path, "r")
# [median_diameters, median_types] = list(csv.reader(file, delimiter=","))

#import details of nerve endings axons 
file_name='P4_endings.csv'
model_path =os.path.join(model_dir,file_name)
file = open(model_path, "r")
[endings_circumference, endings_length] = list(csv.reader(file, delimiter=","))


#Setup EM simulation  
SimA_current=1 #mA
SimB_current=1 #mA
Depths=[100]

F=2000
Beat=100

thresholds_v=[]
thresholds_current=[]

#Select what extra simulations you want to run
run_all_median=0
run_all_endings=1

for Depth in Depths: 

    Nerve='Median Nerve Sensory Axon'
    Axon_name=Nerve +'_neuron'

    #Create an hdf5 file to store the data 
    save_here=os.path.join(path,(str(Depth)+".hdf5"))
    f=h5py.File(save_here,"w")

    #Createa the folders to save the different values in
    E_field=f.create_group('E field')
    Potential=f.create_group('Potential')
    Current_density=f.create_group('Current Density')
    E_field_endings=f.create_group('E field endings')
    Nerve_membrane_potential=f.create_group('Nerve Membrane Potential')
    Endings=f.create_group('Endings Potential')
    #for n in range(final_num_endings):
    #    E_field_endings.create_group(str(n))

    #Set up the simulation
    sim_settings={'Frequency':F,'Magnitude':5,'Name':'sim', 'Gel_conductivity':1/Gel_resistivity}
    sim, results =modelling.EM_Sim(sim_settings,layer_names, layer_materials,Electrodes, Gels, pos_elec, neg_elec)

    #save FEM metadata
    f.attrs['Beat']=Beat
    f.attrs['Frequency']=sim_settings['Frequency']
    f.attrs['Gel conductivity']=1/Gel_resistivity    
    
    #Run the simulation with titration 
    spike_window=math.floor((1/Beat)*1000) #ms
    Sim_duration= [0.01 ,(5/Beat)] #seconds
    Frequencies_neuron=[sim_settings['Frequency'],sim_settings['Frequency']+Beat]
    f.attrs['Stimulation duration']=Sim_duration[1]

    scale=1

    ratio_initial=fun.depth_mod_to_ratio(Depth)*sim_settings['Magnitude']
    Neuro_setup={
        'Temperature':37,'Titration': True, 'Number of Sources':2, 'Depolarisation Threshold':80.0, 'ratio':ratio_initial,
        'Duration': Sim_duration[0], 'Time step': 0.0025,'Frequency':Frequencies_neuron,'scale':scale}

    sim_name=sim_settings['Name']

    threshold_settings={
        'Long Simulation duration':Sim_duration[1], 'threshold':0, 'sim Magnitude': sim_settings['Magnitude'],
        'spike window (ms)':spike_window, 'interspike threshold (ms)':5, 'step sizes (v)': [200,100,50,25,5]}

    #save neuron metadata
    f.attrs['Neuron simulation duration']=threshold_settings['Long Simulation duration']
    f.attrs['Neuron frequencies']=str(Frequencies_neuron)
    f.attrs['interspike threshold (ms)']=threshold_settings['interspike threshold (ms)']
    
    #Get the stim threshold
    scale, titration_factor= modelling.Neuro_threshold_windows(sim_name,Axon_name, Neuro_setup,threshold_settings)
    f.attrs['Titration Voltage']=sim_settings['Magnitude']*titration_factor

    #Delete old EM simulation 
    sim.ClearResults()
    sim.ResetVoxels()
    analysis.ResetAnalysis()
    document.AllSimulations.Remove(sim)

    #Repeat the EM simulation with the new settings 
    sim_settings['Magnitude']=sim_settings['Magnitude']*scale 
    sim, results =modelling.EM_Sim(sim_settings,layer_names, layer_materials,Electrodes, Gels, pos_elec, neg_elec)

    #Repeat the neuro simulation
    Neuro_setup['scale']=1
    Sim_neuro=modelling.Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
    neuron_results=Sim_neuro.Results()

    #Save the stimulation waveform
    time,waveform=modelling.premodulated_waveform(Frequencies_neuron,threshold_settings['Long Simulation duration'],Neuro_setup['Time step'],ratio_initial)
    stimulation_waveform=[time,waveform]
    f.create_dataset('Stimulation waveform',data=stimulation_waveform)

    #Save the stimulation waveform
    time,waveform=modelling.premodulated_waveform(Frequencies_neuron,threshold_settings['Long Simulation duration'],Neuro_setup['Time step'],ratio_initial)
    stimulation_waveform=[time,waveform]
    #f.create_dataset('Stimulation waveform',data=stimulation_waveform)

    #Save the results of the main nerve 
    Axon_name='Main'
    axon_folder=Nerve_membrane_potential.create_group(Axon_name)
    point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
    for point_sensor in point_sensors:
                v=fun.extract_neuron_pointsensor(neuron_results,point_sensor,"v")
                axon_folder.create_dataset(str(point_sensor),data=v)

    Sim_neuro.ClearResults()
    document.AllSimulations.Remove(Sim_neuro)

    #Run the simulation for the rest of the fibres in the Median nerve 
    if run_all_median==1:
        for a in range(len(median_diameters)):
            #Run the neuron simulation
            Axon_name='Median_'+str(a)+'_neuron'
            axon_folder=Nerve_membrane_potential.create_group(Axon_name)
            Sim_neuro=modelling.Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
            neuron_results=Sim_neuro.Results()

            #Save the results 
            point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
            for point_sensor in point_sensors:
                        v=fun.extract_neuron_pointsensor(neuron_results,point_sensor,"v")
                        axon_folder.create_dataset(str(point_sensor),data=v)

            #Clear the simulation
            Sim_neuro.ClearResults()
            document.AllSimulations.Remove(Sim_neuro)
            analysis.ResetAnalysis()

    #Run the simulation for the rest of the fibres in the Median nerve 
    if run_all_endings==1:
        for a in range(len(endings_length)):
            #Run the neuron simulation
            Axon_name='ending_'+str(a)+'_neuron'
            ending_folder=Endings.create_group(Axon_name)
            Sim_neuro=modelling.Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
            neuron_results=Sim_neuro.Results()

            #Save the results 
            point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
            for point_sensor in point_sensors:
                        v=fun.extract_neuron_pointsensor(neuron_results,point_sensor,"v")
                        ending_folder.create_dataset(str(point_sensor),data=v)

            #Clear the simulation
            Sim_neuro.ClearResults()
            document.AllSimulations.Remove(Sim_neuro)
            analysis.ResetAnalysis()
    

    #Get the threshold values
    f.attrs['Voltage threshold']=sim_settings['Magnitude']
    f.attrs['Current threshold']=fun.Get_current(results,sim_settings['Frequency'],sim_settings['Magnitude'])
    
    #Extract overall field 
    overall_field=fun.extract_overall_field(results)
    
    #Get the electric field along the nerve 
    type="EM E(x,y,z,f0)"
    field_nerve=fun.field_along_line(overall_field, Nerve,type, 0.01)
    Field_nerve_values=fun.actual_field_values(field_nerve,type)
    along_axon_grid_A=fun.line_grid_values(field_nerve,type)
    E_field.create_dataset(str(F), data=Field_nerve_values)
    E_field.create_dataset((str(F)+'_grid'), data=along_axon_grid_A)

    #Get the potential along the nerve  
    type="EM Potential(x,y,z,f0)"
    field_nerve=fun.field_along_line(overall_field, Nerve,type, 0.01)
    Field_nerve_values=fun.actual_field_values(field_nerve,type)
    along_axon_grid_A=fun.line_grid_values(field_nerve,type)
    Potential.create_dataset(str(F), data=Field_nerve_values)
    Potential.create_dataset((str(F)+'_grid'), data=along_axon_grid_A)

    #Current density in the skin 
    type="J(x,y,z,f0)"
    current_in_skin=fun.apply_mask (overall_field,type, 'Skin')
    current_values=fun.actual_field_values(current_in_skin,type)
    Current_density.create_dataset(str(F), data=current_values)

    #Save before deleting simulations in case it crashes
    f.close()

    #Delete the simulation and results 
    sim.ClearResults()
    sim.ResetVoxels()
    document.AllSimulations.Remove(sim)
    analysis.ResetAnalysis()
    

    #delete pointsensors 
    entity_names=list(n.Name for n in list(model.AllEntities()))
    point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
    for ps in point_sensors:
        model.AllEntities()[ps].Delete()


