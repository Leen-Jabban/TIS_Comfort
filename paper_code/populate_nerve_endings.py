
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import sys, os
import time
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.simulation.emlf as lf
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewerss
import inspect
import math
from s4l_v1.model import Vec3

#Add main file to path
sys.path.insert(1,r'ADD PATH TO WHERE THE PYTHON SCRIPTS IMPORTED BELOW ARE SAVED')


import s4L_functions as fun
import processing_functions as process 
import modelling_functions as modelling
import csv

model_name='P12_endings.csv'  # Replace with suitable name 

def nerve_ending_points(x1,x2,y1,y2,d1,d2):
    if x1==x2:
        x2=x2+0.000001
    if y1>=y2:

        theta=math.atan(abs((y2-y1)/(x2-x1)))
        xp1=x1+d1*math.cos(theta)
        yp1=y1-d1*math.sin(theta)
    
        alpha=math.atan(d2/d1)
        h=math.sqrt(d1**2+ d2**2)
        xp2=x1+h*math.cos(theta+alpha)
        yp2=y1-h*math.sin(theta+alpha)

    else:
        theta=math.atan(abs((y2-y1)/(x2-x1)))
        xp1=x1+d1*math.cos(theta)
        yp1=y1+d1*math.sin(theta)
    
        alpha=math.atan(d2/d1)
        h=math.sqrt(d1**2+ d2**2)
        xp2=x1+h*math.cos(alpha-theta)
        yp2=y1-h*math.sin(alpha-theta)        
        
    
    return [xp1,yp1,xp2,yp2]

# Enter the points that make up the spline

points, Translation, Rotation=modelling.get_spline_points('Spline')
points=modelling.get_real_coordinates(points,Translation,Rotation)
x_points=[]
y_points=[]
z_points=[]

for point in points:
    x_points.append(point[0])
    y_points.append(point[1])
    z_points.append(point[2])

x=x_points
y=z_points
    
d1=5
d2=1.5
L=[]
endings=[]


for l in range(len(x)-1):
    L.append(math.sqrt((x[l+1]-x[l])**2+ (y[l+1]-y[l])**2))


circumference=sum(L)
n_fibres=math.floor((circumference/d1)+1)
length=[-50,50]
end_point=100
n_layers=math.floor(abs(length[1]-length[0])/d1)+1
n_endings_tot=(n_fibres)*(n_layers)
ending_diameter=7.5

folder=model.CreateGroup('Endings Axon Trajectories')
n=0
endings_circumference=[]
endings_length=[]



for l in range(n_layers):
    s=0
    L_current=0
    xp1=x[0]
    yp1=y[0]
    start_point=length[0]+d1*l
    for i in range(n_fibres):
        if i==0:
            x1=xp1
            y1=yp1
            xp1,yp1,xp2,yp2 = nerve_ending_points(x1,x[s+1],y1,y[s+1],0.1,d2)
                
            endings.append([xp1,yp1,xp2,yp2])
            L_current=L_current+0.1
            axon=model.CreatePolyLine([Vec3(xp2,start_point,yp2),Vec3(xp1,start_point,yp1),Vec3(xp1,end_point,yp1)])
            axon_name='ending_'+ str(n)
            axon.Name = axon_name 
            folder.Add(axon)
            modelling.DiscretizeAxonModel(axon_name, ending_diameter,'sensory',folder)
            endings_circumference.append(i)
            endings_length.append(l)
            n=n+1
        else:

            if L_current+d1<= sum(L[0:s+1]):
                
                x1=xp1
                y1=yp1
                xp1,yp1,xp2,yp2 = nerve_ending_points(x1,x[s+1],y1,y[s+1],d1,d2)
                    
                endings.append([xp1,yp1,xp2,yp2])
                L_current=L_current+d1
                axon=model.CreatePolyLine([Vec3(xp2,start_point,yp2),Vec3(xp1,start_point,yp1),Vec3(xp1,end_point,yp1)])
                axon_name='ending_'+ str(n)
                axon.Name = axon_name 
                folder.Add(axon)
                modelling.DiscretizeAxonModel(axon_name, ending_diameter,'sensory',folder)
                endings_circumference.append(i)
                endings_length.append(l)
                n=n+1

            elif s<len(x)-1:
                s=s+1
                d_last=sum(L[0:s])-L_current
                new_d1=d1-d_last
                x1=x[s]
                y1=y[s]
                xp1,yp1,xp2,yp2 = nerve_ending_points(x1,x[s+1],y1,y[s+1],new_d1,d2)
                    
                endings.append([xp1,yp1,xp2,yp2])
                L_current=L_current+d1
                axon=model.CreatePolyLine([Vec3(xp2,start_point,yp2),Vec3(xp1,start_point,yp1),Vec3(xp1,end_point,yp1)])
                axon_name='ending_'+ str(n)
                axon.Name = axon_name 
                folder.Add(axon)
                modelling.DiscretizeAxonModel(axon_name,ending_diameter,'sensory',folder)
                endings_circumference.append(i)
                endings_length.append(l)

                n=n+1

#Save the results to csv

parent_dir = r'CHOOSE WHERE TO SAVE THE AXON INFORMATION'
path =os.path.join(parent_dir, model_name)
# open the file in the write mode
with open(path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(endings_circumference)
    writer.writerow(endings_length)
f.close()