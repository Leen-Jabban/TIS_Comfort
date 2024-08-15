from s4l_v1 import ReleaseVersion
from s4l_v1 import Unit
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.materials.database as database
import s4l_v1.simulation.emlf as emlf
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewers
import s4l_v1.units as units
import numpy as np 
import itertools



#Define a function to run the simulation 

def apply_mask (Results,type, mask):
	
	#select desired masking location
	selected_entity = model.AllEntities()[mask]
	entities_from_model = all_entities_within_group(selected_entity)

	# Adding a new FieldMaskingFilter
	inputs = [Results.Outputs[type]]
	field_masking_filter = analysis.core.FieldMaskingFilter(inputs=inputs)
	field_masking_filter.SetAllMaterials(False)
	field_masking_filter.SetEntities(entities_from_model)
	field_masking_filter.UpdateAttributes()
	document.AllAlgorithms.Add(field_masking_filter)
	return field_masking_filter

def sum_fields (A,B,type):
	field_A=actual_field_values(A,type)
	field_B=actual_field_values(B,type)
	
	sum=np.add(field_A,field_B)
	return sum

def actual_field_values (field,output_name):
	E_field = field.Outputs[output_name]
	E_field.Update()
	raw_data =(E_field.Data.Field(0))
	
	return raw_data

def line_grid_values(field,output_name):
	E_field=field.Outputs[output_name]
	E_field.Update()
	num_of_points=E_field.Data.Grid.NumberOfPoints
	x=[]
	y=[]
	z=[]
	for point in range(num_of_points):
		this_point= E_field.Data.Grid.GetPoint(point)
		x=np.append(x,this_point[0])
		y=np.append(y,this_point[1])
		z=np.append(z,this_point[2])
	
	raw_data=np.stack([x,y,z])
	return raw_data

def crop_field (Results, upper, lower,name):
	##Note: upper and lower have to be numpy arrays created using numpy.array
	## Adding a new FieldCropFilter
	inputs = [Results.Outputs["EM E(x,y,z,f0)"]]
	field_crop_filter = analysis.core.FieldCropFilter(inputs=inputs)
	field_crop_filter.Name = name
	field_crop_filter.LowerExtent = upper
	field_crop_filter.UpperExtent = lower
	field_crop_filter.UpdateAttributes()
	document.AllAlgorithms.Add(field_crop_filter)

	return field_crop_filter
		
def all_entities_within_group(entity_group):
   
	if isinstance(entity_group, model.EntityGroup):
		return list(itertools.chain.from_iterable(
		all_entities_within_group(e) for e in entity_group.Entities))
	else:
		return [entity_group]

def extract_overall_field(Results):
	em_sensor_extractor = Results["Overall Field"]
	em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"
	em_sensor_extractor.Update()
	document.AllAlgorithms.Add(em_sensor_extractor)
	return em_sensor_extractor

def scale_field(Results, scale_factor):
	# Adding a new EmSensorExtractor
	em_sensor_extractor = Results["Overall Field"]
	em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"
	em_sensor_extractor.Update()
	document.AllAlgorithms.Add(em_sensor_extractor)

	# Adding a new UserDefinedFieldNormalizer
	inputs = [em_sensor_extractor.Outputs["EM E(x,y,z,f0)"]]
	user_defined_field_normalizer = analysis.field.UserDefinedFieldNormalizer(inputs=inputs)
	user_defined_field_normalizer.Target.Value = scale_factor
	user_defined_field_normalizer.UpdateAttributes()
	user_defined_field_normalizer.Update()
	document.AllAlgorithms.Add(user_defined_field_normalizer)
	
	# print('x')
	return user_defined_field_normalizer

def get_tan(field, nerve,max_edge_length):

	# Adding a new ModelToGridFilter
	inputs = []
	model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
	model_to_grid_filter.Name = nerve
	model_to_grid_filter.Entity = model.AllEntities()[nerve]
	model_to_grid_filter.MaximumEdgeLength = max_edge_length, units.Meters
	model_to_grid_filter.UpdateAttributes()
	document.AllAlgorithms.Add(model_to_grid_filter)

	# Adding a new FieldInterpolationFilter
	inputs = [field.Outputs["EM E(x,y,z,f0)"], model_to_grid_filter.Outputs["Line"]]
	field_interpolation_filter = analysis.core.FieldInterpolationFilter(inputs=inputs)
	field_interpolation_filter.UpdateAttributes()
	document.AllAlgorithms.Add(field_interpolation_filter)

	# Adding a new TangentialFieldEvaluator
	inputs = [field_interpolation_filter.Outputs["EM E(x,y,z,f0)"]]
	tangential_field_evaluator = analysis.field.TangentialFieldEvaluator(inputs=inputs)
	tangential_field_evaluator.UpdateAttributes()
	document.AllAlgorithms.Add(tangential_field_evaluator)

	return tangential_field_evaluator

def field_along_line(field_results, line_name,type, max_edge):

	# Adding a new ModelToGridFilter
	inputs = []
	model_to_grid_filter = analysis.core.ModelToGridFilter(inputs=inputs)
	model_to_grid_filter.Name = line_name
	model_to_grid_filter.Entity = model.AllEntities()[line_name]
	model_to_grid_filter.MaximumEdgeLength = max_edge, units.MilliMeters
	model_to_grid_filter.UpdateAttributes()
	document.AllAlgorithms.Add(model_to_grid_filter)

	# Adding a new FieldInterpolationFilter
	inputs = [field_results.Outputs[type], model_to_grid_filter.Outputs["Line"]]
	field_interpolation_filter = analysis.core.FieldInterpolationFilter(inputs=inputs)
	field_interpolation_filter.Name = 'field_along_'+ line_name
	field_interpolation_filter.UpdateAttributes()
	document.AllAlgorithms.Add(field_interpolation_filter)

	return field_interpolation_filter

def extract_neuron_pointsensor(Results, ps_name, output_type):

	sensor_extractor = Results[ps_name]
	output=sensor_extractor.Outputs[output_type]
	document.AllAlgorithms.Add(sensor_extractor)
	output.Update()
	axis=output.Data.Axis
	values =(output.Data.GetComponent(0))
	raw_data=np.stack([axis,values])
	return raw_data

def extract_spike_times(sim):
	
	# Adding a new SensorExtractor
	sensor_extractor = sim["Action Potential Sensor"]
	document.AllAlgorithms.Add(sensor_extractor)

	# Adding a new ActionPotentialEvaluator
	inputs = [sensor_extractor.Outputs["Action Potential"]]
	action_potential_evaluator = analysis.neuron_evaluators.ActionPotentialEvaluator(inputs=inputs)
	action_potential_evaluator.UpdateAttributes()
	document.AllAlgorithms.Add(action_potential_evaluator)

	# Adding a new DataTableHTMLViewer
	inputs = [action_potential_evaluator.Outputs["Action Potential"]]
	data_table_html_viewer = analysis.viewers.DataTableHTMLViewer(inputs=inputs)
	data_table_html_viewer.UpdateAttributes()
	document.AllAlgorithms.Add(data_table_html_viewer)

	input=action_potential_evaluator.Outputs["Action Potential"]
	values=input.Data.ToList()
	return values

def Get_current(simulation,frequency, voltage):
	
	# Adding a new EmSensorExtractor
	em_sensor_extractor = simulation["Overall Field"]
	em_sensor_extractor.FrequencySettings.ExtractedFrequency = u"All"
	document.AllAlgorithms.Add(em_sensor_extractor)

	# Adding a new SarStatisticsEvaluator
	inputs = [em_sensor_extractor.Outputs["EM E(x,y,z,f0)"]]
	sar_statistics_evaluator = analysis.em_evaluators.SarStatisticsEvaluator(inputs=inputs)
	sar_statistics_evaluator.Snapshot = u"{}".format(frequency)
	sar_statistics_evaluator.UpdateAttributes()
	document.AllAlgorithms.Add(sar_statistics_evaluator)
	
	#Extract the value of the power loss overall 
	sar_statistics_evaluator.Update()
	Rows=sar_statistics_evaluator.Outputs["SAR Statistics"].Data.RowKeys()
	y=Rows.index('All Regions')
	Columns=sar_statistics_evaluator.Outputs["SAR Statistics"].Data.ColumnKeys()
	x=Columns.index('Total Loss, W')

	all_SAR_data=sar_statistics_evaluator.Outputs["SAR Statistics"].Data.ToList()

	Power=(all_SAR_data[y][x])

	#Calculate the current in mA 
	Current=Power/(voltage*0.5*1000)
	return Current

def Get_titration_factor(neuron_results):
	# Adding a new SensorExtractor
	Titration_sensor = neuron_results["Titration Sensor"]
	document.AllAlgorithms.Add(Titration_sensor)

	# Adding a new TitrationEvaluator
	inputs = [Titration_sensor.Outputs["Titration"]]
	titration_evaluator = analysis.neuron_evaluators.TitrationEvaluator(inputs=inputs)
	titration_evaluator.UpdateAttributes()
	document.AllAlgorithms.Add(titration_evaluator)
	titration_evaluator.Update()
	titration_factor_data=titration_evaluator.Outputs["Titration factor"]
	titration_factor_data.Update()
	titration_factor=titration_factor_data.Data.GetComponent(0)[0]
	return titration_factor

def depth_mod_to_ratio(depth):
	A=(200-depth)/200
	B=1-A
	ratio=[A, B]
	return ratio