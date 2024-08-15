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
sys.path.insert(1,r'ADD PATH TO WHERE THE PYTHON SCRIPTS IMPORTED BELOW ARE SAVED')

#import commonly used sim4life functions 
import s4L_functions as fun
import processing_functions as process 
import modelling_functions as modellings

#Model name
model_name='P12_median.csv' # Replace with suitable name 

#Median Nerve
median_diameters=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
median_distribution=[17.09,13.51,8.87,4.90,5.30,2.91,10.46,11.39,11.26,5.56,2.78,2.52,1.19,1.19,0.53,0.53] 
Median_nerve_location=[2.52,-4.91] # Remember to update this to correspond to the actual nerve location based on the model. 
Median_radius=1
Median_length=[-100,100]

#Add Median Axons
median_axons=modelling.diameters_list(median_diameters,median_distribution,100)
median_types=modelling.sensory_vs_motor(median_axons,0.15)
modelling.axons('Median',median_axons,median_types,Median_nerve_location,Median_radius,Median_length,2)

#Add main Axon
folder=model.AllEntities()['Median Axon Trajectories']
axon=model.CreatePolyLine([Vec3(Median_nerve_location[0],Median_length[0],Median_nerve_location[1]),Vec3(Median_nerve_location[0],Median_length[1],Median_nerve_location[1])])
axon.Name='Median Nerve Sensory Axon'
folder.Add(axon)
modelling.DiscretizeAxonModel(axon.Name, 19,'sensory',folder)

axon=model.CreatePolyLine([Vec3(Median_nerve_location[0],Median_length[0],Median_nerve_location[1]),Vec3(Median_nerve_location[0],Median_length[1],Median_nerve_location[1])])
axon.Name='Median Nerve Motor Axon'
folder.Add(axon)
modelling.DiscretizeAxonModel(axon.Name, 19,'motor',folder)


#Save the results to csv
parent_dir = r'CHOOSE WHERE TO SAVE THE AXON INFORMATION'
path =os.path.join(parent_dir, model_name)
# open the file in the write mode
with open(path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(median_axons)
    writer.writerow(median_types)

f.close()



