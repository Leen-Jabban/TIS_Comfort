#Import s4l modules
from cmath import sin
from re import S
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.simulation.emlf as emlf
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewers
from s4l_v1.model import Vec3
import s4l_v1.units as units
from s4l_v1 import Unit
import s4l_v1.simulation.neuron as neuron
import s4l_v1.materials.database as database

import random
import numpy as np
import math
import sys, os

#Add main file to path
sys.path.insert(1,r'C:\Users\lj386\OneDrive - University of Bath\PhD\Modelling\Sim4life\code\sim4life')

#import commonly used sim4life functions 
import s4L_functions as fun
import processing_functions as process 

def diameters_list(diameters,distribution, N_fibres):

    # Get the number of fibres for each diameter 
    scaled_distribution=[element*(N_fibres/100) for element in distribution]
    rounded_distribution=[int(math.floor(element)) for element in scaled_distribution]
    check_sum=sum(rounded_distribution)

    # sort based on highest distribution first (for the next step)
    rounded_distribution,diameters=zip(*sorted(zip(rounded_distribution,diameters),reverse=True))
    
    # Change back to a list 
    rounded_distribution=list(rounded_distribution)
    diameters=list(diameters)

    # A way to get to the right number. Not great, but those numbers are an estimate anyway
    n=0
    while check_sum < N_fibres:
        rounded_distribution[n]=rounded_distribution[n]+1
        n=n+1
        check_sum=sum(rounded_distribution)

    # Create a list of fibre diameters
    result=[] 

    for i in range(len(diameters)):
        n=0
        while n< rounded_distribution[i]:
            result.append(diameters[i])
            n=n+1
          
    random.shuffle(result)
    return result

def layered_cylinder(names,thickness,material,length):
    current=0
    for l in range(len(names)):
        layer = model.CreateSolidTube( Vec3(0,0,0),Vec3(0,0,length),thickness[l]+current,current)
        layer.MaterialName=material[l]
        layer.Name = names[l]
        current+=thickness[l]

def axons(nerve_name,Diameters,type,centre,R,length,length_axis):
    folder=model.CreateGroup(nerve_name+' Axon Trajectories')
    a=0
    while a<len(Diameters):
        x=2*R*(-0.5+np.random.rand())
        y=2*R*(-0.5+np.random.rand())
        if (x*x+y*y<R*R):
            if length_axis==3:
                axon=model.CreatePolyLine([Vec3(centre[0]+x,centre[1]+y,length[0]),Vec3(centre[0]+x,centre[1]+y,length[1])])
            elif length_axis==2:
                axon=model.CreatePolyLine([Vec3(centre[0]+x,length[0],centre[1]+y),Vec3(centre[0]+x,length[1],centre[1]+y)])
            axon.Name=nerve_name+'_'+str(a)
            folder.Add(axon)
            DiscretizeAxonModel(axon.Name, Diameters[a],type[a],folder)
            a+=1

def sensory_vs_motor(Diameters, perc_motor):
    num_motor=len(Diameters)*perc_motor
    type=[]
    for a in range(len(Diameters)):
        if a<=num_motor:
            type.append('motor')
        else:
            type.append('sensory')
    return type

def DiscretizeAxonModel(Axon_name, Diameter, type,folder):
    axon_entity = model.AllEntities()[Axon_name]

    if type=='motor':
        model_properties=model.MotorMrgNeuronProperties()
    elif type=='sensory':
        model_properties=model.SensoryMrgNeuronProperties()
    else:
        model_properties=model.MotorNeuronProperties()
    
    model_properties.AxonDiameter=Diameter
    discretized_axon = model.CreateAxonNeuron(axon_entity,model_properties)
    discretized_axon.Name = Axon_name +'_neuron'
    folder.Add(discretized_axon)


def fibre_distribution(layers_thickness,section_length, distribution):
    arm_radius=sum(layers_thickness)
    circumference=arm_radius*math.pi*2
    area=circumference*section_length
    tot_fibres=distribution*area
    fibres_circum=math.sqrt(circumference*tot_fibres/section_length)
    fibres_length=tot_fibres/fibres_circum
    num_fibres=[fibres_circum, fibres_length]
    return num_fibres

def distribute_nerve_endings(num_fibres, section_length,arm_length,R_tot,skin, SAT, dermis,diameters):
    folder=model.CreateGroup('Endings Axon Trajectories')
    start_point=(arm_length+section_length)/2
    d_l=section_length/(num_fibres[1]-1)
    d_theta=2*math.pi/(num_fibres[0])
    l=start_point
    theta=0
    a=0

    while l>=start_point-section_length:
        theta=0
        while theta<2*math.pi:
            points=[]
            #Ending point 
            R=R_tot-skin+dermis*random.random()
            x=R*math.sin(theta)
            y=R*math.cos(theta)
            points.append(Vec3(x,y,l))

            #mid-point 
            R=R_tot-skin-SAT*0.25*random.random()
            x=R*math.sin(theta)
            y=R*math.cos(theta)
            points.append(Vec3(x,y,l))

            #deep point
            points.append(Vec3(x,y,0))  
            
            # Create the axon
            axon=model.CreatePolyLine(points)
            axon.Name='ending_'+str(a)
            folder.Add(axon)
            DiscretizeAxonModel(axon.Name, diameters[a],'sensory',folder)
            a+=1    
            theta+=d_theta

        
        l=l-d_l
        return a

def place_electrodes(number,radius,length_spacing, arm_radius, arm_length, thickness):
    section_length=(number[1]-1)*length_spacing
    starting_point=(arm_length-section_length)/2
    gels=model.CreateGroup('Gels')
    electrodes=model.CreateGroup('Electrodes')
    identifiers=[]
    for r in range(number[0]):
        angle=r* (2*math.pi/number[0])

        x_in=(arm_radius-thickness[0])*math.sin(angle)
        y_in=(arm_radius-thickness[0])*math.cos(angle)
        x_mid=(arm_radius+thickness[1])*math.sin(angle)
        y_mid=(arm_radius+thickness[1])*math.cos(angle)
        x_out=(arm_radius+thickness[1]+thickness[2])*math.sin(angle)
        y_out=(arm_radius+thickness[1]+thickness[2])*math.cos(angle)

        for l in range(number[1]):
            length=length_spacing*l
            z=starting_point+length

            identifier='('+str(r)+','+ str(l) + ')'
            identifiers.append(identifier)
            gel=model.CreateSolidCylinder(Vec3(x_in,y_in,z), Vec3(x_mid,y_mid,z),radius)
            gel.Name= 'Gel'+identifier
            gels.Add(gel)

            electrode=model.CreateSolidCylinder( Vec3(x_mid,y_mid,z),Vec3(x_out,y_out,z),radius)
            electrode.Name='Electrode'+identifier
            electrodes.Add(electrode)
    return identifiers

def EM_Sim(setup,Entities, Materials,Electrodes, Gels, pos_elec, neg_elec):

    #Create the simulation 
    sim = emlf.ElectroQsOhmicSimulation()
    sim.Name = setup['Name']

    #Set simulation frequncy
    sim.SetupSettings.Frequency = setup['Frequency']

    #Loop through the different material settings 
    Arm=[]
    for e in range(len(Entities)):
        entitiy=model.AllEntities()[Entities[e]]
        Arm.append(entitiy)
        material_settings = emlf.MaterialSettings()
        components = [entitiy]
        mat = database["IT'IS 4.1"][Materials[e]]
        sim.LinkMaterialWithDatabase(material_settings, mat)
        sim.Add(material_settings, components)

    #Add Gels material settings 
    material_settings = emlf.MaterialSettings()
    component=Gels
    material_settings.Name='Gels'
    material_settings.ElectricProps.Conductivity = setup['Gel_conductivity'], Unit("S/m")
    sim.Add(material_settings, component)

    #Set boundary conditions
    boundary_pos = sim.AddBoundarySettings(pos_elec)
    boundary_pos.DirichletValue = setup['Magnitude']
    boundary_pos.Name = 'Dirichlet +'

    boundary_neg = sim.AddBoundarySettings(neg_elec)
    boundary_neg.DirichletValue = -1*setup['Magnitude']
    boundary_neg.Name = 'Dirichlet -'

    # Grid
    sim.GlobalGridSettings.ManualDiscretization = True
    sim.GlobalGridSettings.MaxStep = np.array([1, 1, 1]), units.MilliMeters
    sim.GlobalGridSettings.Resolution = np.array([0.5, 0.5, 0.5]), units.MilliMeters
    sim.GlobalGridSettings.ManualPadding = True
    sim.GlobalGridSettings.BottomPadding = sim.GlobalGridSettings.TopPadding = np.array([10, 10, 10]), units.MilliMeters

    manual_grid_settings = sim.AddManualGridSettings(Arm)
    manual_grid_settings.MaxStep = np.array([0.65, 0.65, 0.65]), units.MilliMeters
    manual_grid_settings.Resolution = np.array([0.65, 0.65, 0.65]), units.MilliMeters

    anual_grid_settings = sim.AddManualGridSettings(Gels)
    manual_grid_settings.MaxStep = np.array([0.65, 0.65, 0.65]), units.MilliMeters
    manual_grid_settings.Resolution = np.array([0.65, 0.65, 0.65]), units.MilliMeters

    
    # Voxels
    # Removing AutomaticVoxelerSettings Automatic Voxeler Settings
    automatic_voxeler_settings = [x for x in sim.AllSettings if isinstance(x, emlf.AutomaticVoxelerSettings) and x.Name == "Automatic Voxeler Settings"][0]
    sim.RemoveSettings(automatic_voxeler_settings)

    # Adding a new ManualVoxelerSettings for the Arm 
    manual_voxeler_settings = emlf.ManualVoxelerSettings()
    manual_voxeler_settings.Name = "Arm"
    manual_voxeler_settings.Priority = 2
    sim.Add(manual_voxeler_settings, Arm)

    # Adding a new ManualVoxelerSettings for the Electrodes 
    manual_voxeler_settings = emlf.ManualVoxelerSettings()
    components = Electrodes
    manual_voxeler_settings.Name = "Electrodes"
    manual_voxeler_settings.Priority = 3
    sim.Add(manual_voxeler_settings, components)

    # Adding a new ManualVoxelerSettings for the Gels 
    manual_voxeler_settings = emlf.ManualVoxelerSettings()
    components = Gels
    manual_voxeler_settings.Name = "Gels"
    manual_voxeler_settings.Priority = 0
    sim.Add(manual_voxeler_settings, components)

    # Solver settings
    sim.SolverSettings.PredefinedTolerances = sim.SolverSettings.PredefinedTolerances.enum.High

    # Update the materials with the new frequency parameters
    sim.UpdateAllMaterials()

    # Update the grid with the new parameters
    sim.UpdateGrid()

    #Add simulation to document 
    document.AllSimulations.Add( sim )
    sim.CreateVoxels()
    sim.RunSimulation()
    results = sim.Results()
    return sim, results

def current_sim(setup,Entities, Materials,Electrodes, Gels, pos_elec, neg_elec):

    #Run the simulation with a 1v setup
    Current=setup['Magnitude']
    setup['Magnitude']=1
    sim, results =EM_Sim(setup,Entities, Materials,Electrodes, Gels, pos_elec, neg_elec)

    #Get the current 
    Current_sim=fun.Get_current(results,setup['Frequency'],2)
    scale=Current_sim/Current

    #Delete the previous simulation
    document.AllSimulations.Remove(sim)

    #Repeat the simulation with the correct voltage
    setup['Magnitude']=scale
    sim, results=EM_Sim(setup,Entities, Materials,Electrodes, Gels, pos_elec, neg_elec)

    return sim, results

def get_largest_sensory_axon(diameters,type):
    max_diameter=max(diameters)
    index=index=np.where(np.array(diameters) == max_diameter)[0]
    axons=[]
    for i in index:
        if type[i]=='sensory':
            axons.append(i)
    return axons

def Neuro_sim(sim_name,Axon_name, setup):

    # Creating the simulation
    simulation = neuron.Simulation()
    simulation.Name = Axon_name +' Sim'

    duration=setup['Duration']
   
    # Mapping the components and entities
    entities = model.AllEntities()
    axon = entities[Axon_name]

    # Setup Settings
    setup_settings = simulation.SetupSettings
    setup_settings.PerformTitration = True
    setup_settings.DepolarizationDetection = setup_settings.DepolarizationDetection.enum.Threshold
    setup_settings.DepolarizationThreshold = setup['Depolarisation Threshold']

    # Simulation links

    # Will Titration be performed?
    setup_settings = simulation.SetupSettings
    setup_settings.PerformTitration = setup['Titration']

    # Neuron Settings
    neuron_settings = neuron.AutomaticAxonNeuronSettings()
    neuron_settings.Temperature = setup['Temperature']
    simulation.Add(neuron_settings,axon)

    simulation.LoadModel()

    # Adding a new SourceSettings
    for n in range(setup['Number of Sources']):
        field= document.AllSimulations[sim_name[n]]
        source= simulation.AddSource(field, "Overall Field")
        source.Name = "Field_"+ str(n)
        source.PulseType = source.PulseType.enum.Sinusoidal
        source.FrequencySine = setup['Frequency'][n], units.Hz
        source.AmplitudeSine = setup['scale'][n]
        source.NumberOfHalfPeriodsSine =int(setup['Frequency'][n]*duration)*2 
        del field

    #Add point sensor
    sectionNames = simulation.GetSectionNames(axon)
    for n in sectionNames:
        #Only compute at nodes
        is_node=n.startswith('node')
        if is_node==True: 
            point_sensor_settings = simulation.AddPointSensor(axon)
            point_sensor_settings.SectionName = n
            point_sensor_settings.RecordEExtracellular = True
            point_sensor_settings.RecordIMembrane = True
            point_sensor_settings.RecordV= True
            point_sensor_settings.RelativePosition = 0.5

    # Editing SolverSettings "Solver"
    solver_settings = simulation.SolverSettings
    solver_settings.Duration = duration, units.Seconds
    solver_settings.TimeStep=setup['Time step'] , units.MilliSecond

    # Add the simulation to the UI
    document.AllSimulations.Add( simulation )
    simulation.RunSimulation()

    return simulation

def Neuro_threshold(sim_name,Axon_name, Neuro_setup,threshold_settings):
    
    threshold=threshold_settings['threshold']
    spike_window=threshold_settings['spike window (ms)']
    interspike_threshold=threshold_settings['interspike threshold (ms)']
    step_sizes=threshold_settings['step sizes (v)']
    sim_magnitude=threshold_settings['sim Magnitude']
    minimum_spikes=math.floor((Neuro_setup['Duration']*1000/interspike_threshold))
    
    Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
    neuron_results=Sim_neuro.Results()

    #Get the tiration factor 
    titration_factor=fun.Get_titration_factor(neuron_results).item()
    titration_voltage=sim_magnitude*titration_factor

    #Delete the simulation
    Sim_neuro.ClearResults()
    document.AllSimulations.Remove(Sim_neuro)

    #delete points sensors 
    entity_names=list(n.Name for n in list(model.AllEntities()))
    point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
    for ps in point_sensors:
        model.AllEntities()[ps].Delete()
    del neuron_results

    #Repeat simulation with scaled input
    Neuro_setup['Duration']=threshold_settings['Long Simulation duration']
    Neuro_setup['scale'][0]=math.floor(titration_voltage)/sim_magnitude
    Neuro_setup['Titration']=False
    

    for step_size in step_sizes:

        #Next simulation
        Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
        neuron_results=Sim_neuro.Results()

        #Get the time between intervals at the last node 
        point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
        last_node=point_sensors[0] #Is there a better way to do this? 
        v=fun.extract_neuron_pointsensor(neuron_results,last_node,"v")

        #Get the interspike interval
        interspike, number_spikes=process.Av_interspike_interval(v,threshold,spike_window)

        while (interspike>interspike_threshold) or (number_spikes<= minimum_spikes) :

            #Get the next step
            current_voltage=Neuro_setup['scale'][0]*sim_magnitude
            scale=(current_voltage+step_size)/sim_magnitude
            Neuro_setup['scale'][0]=scale

            #Delete simulation
            Sim_neuro.ClearResults()
            document.AllSimulations.Remove(Sim_neuro)

            #delete pointsensors 
            entity_names=list(n.Name for n in list(model.AllEntities()))
            point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
            for ps in point_sensors:
                model.AllEntities()[ps].Delete()

            del Sim_neuro, neuron_results

            #Next simulation
            Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
            neuron_results=Sim_neuro.Results()

            #Get the time between intervals at the last node 
            point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
            last_node=point_sensors[0] #Is there a better way to do this? 
            v=fun.extract_neuron_pointsensor(neuron_results,last_node,"v")

            #Get the interspike interval (in seconds)
            threshold=0
            spike_window=0.001*2 #in seconds
            interspike, number_spikes=process.Av_interspike_interval(v,threshold,spike_window)
        
        current_voltage=Neuro_setup['scale'][0]*sim_magnitude
        Neuro_setup['scale'][0]=(current_voltage-step_size)/sim_magnitude

        #Delete simulation
        Sim_neuro.ClearResults()
        document.AllSimulations.Remove(Sim_neuro)

        #delete pointsensors 
        entity_names=list(n.Name for n in list(model.AllEntities()))
        point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
        for ps in point_sensors:
            model.AllEntities()[ps].Delete()

        del Sim_neuro, neuron_results

    return scale, titration_factor

def Neuro_threshold_IFS(sim_name,Axon_name, Neuro_setup,threshold_settings):
    
    threshold=threshold_settings['threshold']
    spike_window=threshold_settings['spike window (ms)']
    interspike_threshold=threshold_settings['interspike threshold (ms)']
    step_sizes=threshold_settings['step sizes (v)']
    sim_magnitude=threshold_settings['sim Magnitude']
    minimum_spikes=threshold_settings['Minimum number of Spikes']
    
    Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
    neuron_results=Sim_neuro.Results()

    #Get the tiration factor 
    titration_factor=fun.Get_titration_factor(neuron_results).item()
    titration_voltage=sim_magnitude*titration_factor

    #Delete the simulation
    Sim_neuro.ClearResults()
    document.AllSimulations.Remove(Sim_neuro)

    #delete points sensors 
    entity_names=list(n.Name for n in list(model.AllEntities()))
    point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
    for ps in point_sensors:
        model.AllEntities()[ps].Delete()
    del neuron_results

    #Repeat simulation with scaled input
    Neuro_setup['Duration']=threshold_settings['Long Simulation duration']
    for n in range(Neuro_setup['Number of Sources']):
        Neuro_setup['scale'][n]=math.floor(titration_voltage)/sim_magnitude
    Neuro_setup['Titration']=False
    current_voltage=titration_voltage
    

    for step_size in step_sizes:

        #Next simulation
        Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
        neuron_results=Sim_neuro.Results()

        #Get the time between intervals at the last node 
        point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
        last_node=point_sensors[0] #Is there a better way to do this? 
        v=fun.extract_neuron_pointsensor(neuron_results,last_node,"v")

        #Get the interspike interval
        interspike, number_spikes=process.Av_interspike_interval(v,threshold,spike_window)

        while ((interspike<interspike_threshold) or (number_spikes<minimum_spikes)) :

            #Get the next step
            current_voltage=Neuro_setup['scale'][0]*sim_magnitude
            scale=(current_voltage+step_size)/sim_magnitude
            for n in range(Neuro_setup['Number of Sources']):
                Neuro_setup['scale'][n]=scale

            #Delete simulation
            Sim_neuro.ClearResults()
            document.AllSimulations.Remove(Sim_neuro)

            #delete pointsensors 
            entity_names=list(n.Name for n in list(model.AllEntities()))
            point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
            for ps in point_sensors:
                model.AllEntities()[ps].Delete()

            del Sim_neuro, neuron_results

            #Next simulation
            Sim_neuro=Neuro_sim(sim_name,Axon_name,Neuro_setup)
            neuron_results=Sim_neuro.Results()

            #Get the time between intervals at the last node 
            point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
            last_node=point_sensors[0] #Is there a better way to do this? 
            v=fun.extract_neuron_pointsensor(neuron_results,last_node,"v")

            #Get the interspike interval (in seconds)
            threshold=0
            spike_window=0.001*2 #in seconds
            interspike, number_spikes=process.Av_interspike_interval(v,threshold,spike_window)
        
        current_voltage=Neuro_setup['scale'][0]*sim_magnitude
        for n in range(Neuro_setup['Number of Sources']):
            Neuro_setup['scale'][n]=(current_voltage-step_size)/sim_magnitude

        #Delete simulation
        Sim_neuro.ClearResults()
        document.AllSimulations.Remove(Sim_neuro)

        #delete pointsensors 
        entity_names=list(n.Name for n in list(model.AllEntities()))
        point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
        for ps in point_sensors:
            model.AllEntities()[ps].Delete()

        del Sim_neuro, neuron_results

    return scale, titration_factor

def premodulated_waveform(f,t_max,dt, ratio):
    dt_s=dt/1000
    time_s=np.arange(0,t_max,dt_s)
    time_ms=time_s*1000
    waveform=ratio[0] * np.sin(2 * np.pi * f[0] * time_s)+ ratio[1] * np.sin(2 * np.pi * f[1] * time_s)

    return time_ms, waveform

def Neuro_sim_premodulated(sim_name,Axon_name, setup):

    # Creating the simulation
    simulation = neuron.Simulation()
    simulation.Name = Axon_name +' Sim'

    duration=setup['Duration']
   
    # Mapping the components and entities
    entities = model.AllEntities()
    axon = entities[Axon_name]

    # Setup Settings
    setup_settings = simulation.SetupSettings
    setup_settings.PerformTitration = True
    setup_settings.DepolarizationDetection = setup_settings.DepolarizationDetection.enum.Threshold
    setup_settings.DepolarizationThreshold = setup['Depolarisation Threshold']

    # Simulation links

    # Will Titration be performed?
    setup_settings = simulation.SetupSettings
    setup_settings.PerformTitration = setup['Titration']

    # Neuron Settings
    neuron_settings = neuron.AutomaticAxonNeuronSettings()
    neuron_settings.Temperature = setup['Temperature']
    simulation.Add(neuron_settings,axon)

    simulation.LoadModel()

    # Adding a new SourceSettings
    time,waveform=premodulated_waveform(setup['Frequency'],duration,setup['Time step'],setup['ratio'])

    field= document.AllSimulations[sim_name]
    source= simulation.AddSource(field, "Overall Field")
    source.Name = "Field_"+ str(0)
    source.PulseType = source.PulseType.enum.Vectors
    
    source.TimeVector = time
    source.ModulationVector = waveform* setup['scale']
    del field

    #Add point sensor
    sectionNames = simulation.GetSectionNames(axon)
    for n in sectionNames:
        #Only compute at nodes
        is_node=n.startswith('node')
        if is_node==True: 
            point_sensor_settings = simulation.AddPointSensor(axon)
            point_sensor_settings.SectionName = n
            point_sensor_settings.RecordEExtracellular = True
            point_sensor_settings.RecordIMembrane = True
            point_sensor_settings.RecordV= True
            point_sensor_settings.RelativePosition = 0.5

    # Editing SolverSettings "Solver"
    solver_settings = simulation.SolverSettings
    solver_settings.Duration = duration, units.Seconds
    solver_settings.TimeStep=setup['Time step'] , units.MilliSecond

    # Add the simulation to the UI
    document.AllSimulations.Add( simulation )
    simulation.RunSimulation()

    return simulation

def Neuro_threshold_windows(sim_name,Axon_name, Neuro_setup,threshold_settings): 
    
    
    threshold=threshold_settings['threshold']
    spike_window=threshold_settings['spike window (ms)']
    interspike_threshold=threshold_settings['interspike threshold (ms)']
    step_sizes=threshold_settings['step sizes (v)']
    sim_magnitude=threshold_settings['sim Magnitude']
    
    Sim_neuro=Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
    neuron_results=Sim_neuro.Results()

    #Get the tiration factor 
    titration_factor=fun.Get_titration_factor(neuron_results).item()
    titration_voltage=sim_magnitude*titration_factor

    #Delete the simulation
    Sim_neuro.ClearResults()
    document.AllSimulations.Remove(Sim_neuro)

    #delete points sensors 
    entity_names=list(n.Name for n in list(model.AllEntities()))
    point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
    for ps in point_sensors:
        model.AllEntities()[ps].Delete()
    del neuron_results

    #Repeat simulation with scaled input
    Neuro_setup['Duration']=threshold_settings['Long Simulation duration']
    Neuro_setup['scale']=math.floor(titration_voltage)/sim_magnitude
    Neuro_setup['Titration']=False
    current_voltage=titration_voltage

    for step_size in step_sizes:

        #Next simulation
        Sim_neuro=Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
        neuron_results=Sim_neuro.Results()

        #Get the time between intervals at the last node 
        point_sensors = [sensor for sensor in neuron_results.keys() if '@node' in sensor]
        #Look at the last node 
        last_node=0
        for node in point_sensors:
            node_number=process.get_node(node)
            if node_number>last_node:
                last_node=node_number

        node_before_last=str(last_node-1)
        saved_node=next((s for s in point_sensors if node_before_last in s), None)
        v=fun.extract_neuron_pointsensor(neuron_results,saved_node,"v")

        #Get the interspike interval
        output=process.spike_per_window_check(v,threshold,spike_window)

        while output==0 :

            #Get the next step
            current_voltage=Neuro_setup['scale']*sim_magnitude
            scale=(current_voltage+step_size)/sim_magnitude
            Neuro_setup['scale']=scale

            #Delete simulation
            Sim_neuro.ClearResults()
            document.AllSimulations.Remove(Sim_neuro)

            #delete pointsensors 
            entity_names=list(n.Name for n in list(model.AllEntities()))
            point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
            for ps in point_sensors:
                model.AllEntities()[ps].Delete()

            del Sim_neuro, neuron_results

            #Next simulation
            Sim_neuro=Neuro_sim_premodulated(sim_name,Axon_name,Neuro_setup)
            neuron_results=Sim_neuro.Results()

            #Get the time between intervals at the last node 
            v=fun.extract_neuron_pointsensor(neuron_results,saved_node,"v")

            #Get the interspike interval (in seconds)
            output=process.spike_per_window_check(v,threshold,spike_window)
        
        current_voltage=Neuro_setup['scale']*sim_magnitude
        Neuro_setup['scale']=(current_voltage-step_size)/sim_magnitude

        #Delete simulation
        Sim_neuro.ClearResults()
        analysis.ResetAnalysis()
        document.AllSimulations.Remove(Sim_neuro)

        #delete pointsensors 
        entity_names=list(n.Name for n in list(model.AllEntities()))
        point_sensors = [sensor for sensor in entity_names if '@node' in sensor]
        for ps in point_sensors:
            model.AllEntities()[ps].Delete()

        del Sim_neuro, neuron_results

    return scale, titration_factor


def get_spline_points(spline_name):
    spline = model.AllEntities()[spline_name]
    points = [p.Value for p in spline.Parameters.Points]
    Translation=spline.Transform.Translation
    Rotation=spline.Transform.DecomposeRotation

    return points,Translation,Rotation


def get_real_coordinates(points, Translation, Rotation):
    
    real_coordinates=[]
    for n in range(len(points)):
        #start with x axis rotation
        cos_theta_x=math.cos(Rotation[0])
        sin_theta_x=math.sin(Rotation[0])
        x_rotation_matrix=np.array([[1,0,0],[0,cos_theta_x,-sin_theta_x],[0,sin_theta_x,cos_theta_x]])
        point=np.array(points[n])
        new_points=np.matmul(x_rotation_matrix,point)
        
        #y axis roration
        cos_theta_y=math.cos(Rotation[1])
        sin_theta_y=math.sin(Rotation[1])
        y_rotation_matrix=np.array([[cos_theta_y,0,sin_theta_y],[0,1,0],[-sin_theta_y,0,cos_theta_y]])
        new_points=np.matmul(y_rotation_matrix,new_points)
        
        #z axis roration
        cos_theta_z=math.cos(Rotation[2])
        sin_theta_z=math.sin(Rotation[2])
        z_rotation_matrix=np.array([[cos_theta_z,-sin_theta_z,0],[sin_theta_z,cos_theta_z,0],[0,0,1]])
        new_points=np.matmul(z_rotation_matrix,new_points)
        
        real_coordinates.append(new_points+ Translation)
                
    return real_coordinates

def AxonPropertyValues(entity,index):
    # credit to Harshita Tangella March 2023

    # 0 - AxonDiameter; 1 - NodeNodeSeparation; 2 - NumberOfMyelinLamellae; 3 - NodeLength; 4 - NodeDiameter
    # 5 - MYSALength; 6 - MYSADiameter; 7 - MYSAPeriaxonalSpaceWidth; 8 - FLUTLength; 9 - FLUTDiameter
    # 10 - FLUTPeriaxonalSpaceWidth; 11 - STINLength; 12 - STINDiameter; 13 - STInPeriaxonalSpaceWidth
    properties = model.GetAxonProperties(entity)
    parameter = properties[index]
    parameter_value = parameter.Value
    return parameter_value