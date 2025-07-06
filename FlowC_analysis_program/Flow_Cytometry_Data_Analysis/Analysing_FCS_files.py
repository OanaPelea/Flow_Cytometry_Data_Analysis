#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 13:00:13 2020

@author: oanapelea
"""
import FlowCal
import pylab
import fcsparser
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from functions_fcs_file_analysis import define_gating_strategy_Blanc_condition
from functions_fcs_file_analysis import define_gating_strategy_REP_condition_iBlue_ECFP
from functions_fcs_file_analysis import output_log_x_values_and_min_and_max_values
from functions_fcs_file_analysis import selecting_cells_from_given_gate_and_input_file
from functions_FlowC_analysis import guide_colour_code
from functions_FlowC_analysis import plot_graph
import matplotlib.patches as patches
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse

gate_position_multiplication_factor=20 

##########################################################################################
##########################################################################################
input_file='input_FCS_files_and_lists/0_Files_inputs_for_fcs_analysis_program/00_example_data.csv'

output_file_individual_analysis=open('output_numerical_data/00_example_data_numerical_info.csv', 'a') # will need to change append 

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#####do not change
intermediary_filename='output_numerical_data/0_intermediary_output.csv'
output_file_all_data_analysis=open('output_numerical_data/Centralised_FlowC_data.csv', 'a') #delete test

##########################################################################################
##########################################################################################
###start intermediary comments
intermediary_output_file=open(intermediary_filename, 'w')
param_reporter_colour='405-450/50-A' ## '405-450/50-A'
#########################################################################################   


input_files_replicate_1=[]
input_files_replicate_2=[]
input_files_replicate_3=[]

names_conditions=[]
full_condition_details_and_indexes=[]

for line in open(input_file).readlines():
   if (",,,,,,,,,,,,,,,,," not in line) and ("Condition_name" not in line):
     #print line
     f=line.split(',')
     if len(f[1])>0:
       condition_name=str(f[1])
       names_conditions.append(condition_name)
       full_condition_details_and_indexes.append(f[0: 19])

     total_number_replicates=0
     if len(str(f[9]))>2:
      input_files_replicate_1.append(str(f[9]))
      total_number_replicates=total_number_replicates+1

     if len(str(f[10]))>2:
      input_files_replicate_2.append(str(f[10]))
      total_number_replicates=total_number_replicates+1

     if len(str(f[11]))>2:
      input_files_replicate_3.append(str(f[11]))
      total_number_replicates=total_number_replicates+1
       
input_files_replicate=[input_files_replicate_1,input_files_replicate_2, input_files_replicate_3]

#determine_number_guides 
index_g=0
full_gRNA_names=[]
gRNA_colours=[]
gRNA_names=[]
gRNA_start_indexes=[]
list_counts_number_conditions_for_individual_gRNAs=[]
index=0
while index < len(full_condition_details_and_indexes):
    condition=full_condition_details_and_indexes[index]
    CTS_used=condition[3]    
    
    if ('Blanc' in condition) and ('REP' in full_condition_details_and_indexes[index+1]) :
       next_condition=full_condition_details_and_indexes[index+1]
       CTS_used_next_condition=next_condition[3]
       gRNA_name=CTS_used_next_condition[0:2]#change
       gRNA_with_c=gRNA_name+'c'
       if gRNA_with_c in next_condition:
            gRNA_name=gRNA_with_c
       
    else:
        gRNA_name=CTS_used[0:2]#change
        gRNA_with_c=gRNA_name+'c'
        if gRNA_with_c in condition:
            gRNA_name=gRNA_with_c
    full_gRNA_names.append(gRNA_name)
    gRNA_colours.append(guide_colour_code(gRNA_name))
    if gRNA_name not in gRNA_names:
      gRNA_names.append(gRNA_name)
      gRNA_start_indexes.append(index) 
    index=index+1    
             
for gRNA in gRNA_names:
    guide_occurence=full_gRNA_names.count(gRNA)
    list_counts_number_conditions_for_individual_gRNAs.append(guide_occurence)
                                                              
#Other input parameters
HEK_cell_fraction=0.5
single_cell_gate_fraction=0.95

#Checking if input is correct

if total_number_replicates==1:
    if len(input_files_replicate_1)==len(names_conditions):
        print "Input files and condition names have desired format"
    else:
        print "!!!!Check input file and condition names"
        
if total_number_replicates==2:
    if (len(input_files_replicate_1)==len(names_conditions)) and (len(input_files_replicate_2)==len(names_conditions)):
        print "Input files and condition names have desired format"
    else:
        print "!!!!Check input file and condition names"

if total_number_replicates==3:
    if (len(input_files_replicate_1)==len(names_conditions)) and (len(input_files_replicate_2)==len(names_conditions)) and (len(input_files_replicate_3)==len(names_conditions)):
        print "Input files and condition names have desired format"
    else:
        print "!!!!Check input file and condition names"
#angles of rotation
        
        
#output lists
list_number_HEK293T_cells_all_replicates=[]
list_number_single_cells_all_replicates=[]
list_number_transfected_cells_all_replicates=[]
list_number_percentage_activated_transfected_all_replicates=[]
list_median_ECFP_transfected_cells_all_replicates=[]

#global variables
global_HEK_centre=[0,0]
global_HEK_angle=0
global_HEK_axes=[0,0]

global_single_cell_centre=[0,0]
global_single_cell_angle=0
global_single_cell_axes=[0,0]

#running code

index_number_replicates=0    

while index_number_replicates < total_number_replicates:
 files_corresponding_replicate=input_files_replicate[index_number_replicates]
 
 list_number_HEK293T_cells=[]
 list_number_single_cells=[]
 list_number_transfected_cells=[]
 list_number_percentage_activated_transfected=[]
 list_median_ECFP_transfected_cells=[]
 
 index_guide_number=0
 
 while index_guide_number < len(gRNA_start_indexes):
  
  gRNA_name=gRNA_names[index_guide_number]
  number_occurence_conditions_with_particular_gRNA=list_counts_number_conditions_for_individual_gRNAs[index_guide_number]

  #Determine blanc file 
  index_blanc=gRNA_start_indexes[index_guide_number]
  while index_blanc<len(names_conditions):
    if ('Blanc' in full_condition_details_and_indexes[index_blanc]) and (gRNA_name in full_gRNA_names[index_blanc]):
       Blanc_filename_and_path=str(files_corresponding_replicate[index_blanc])
       
       final_index_blanc=index_blanc

    index_blanc=index_blanc+1

  #Determine REP file
  index_REP=gRNA_start_indexes[index_guide_number]
  while index_REP<len(names_conditions):
    if ('REP' in full_condition_details_and_indexes[index_REP]) and (gRNA_name in full_gRNA_names[index_REP]):
       REP_filename_and_path=str(files_corresponding_replicate[index_REP])     
       print REP_filename_and_path
    index_REP=index_REP+1
   
  #Extracting transfection and minimal activation gate from the blanc
  blanc_max_iBlue_and_ECFP_values=define_gating_strategy_Blanc_condition(Blanc_filename_and_path, HEK_cell_fraction, single_cell_gate_fraction, '640-670/14-A', param_reporter_colour)
  transfection_gate_cuttoff=blanc_max_iBlue_and_ECFP_values[0]
  initial_value_to_assess_for_activation_gate=float(blanc_max_iBlue_and_ECFP_values[1])
  print initial_value_to_assess_for_activation_gate
  
  #Extracting max ECFP activation in the REP using blanc function
  max_value_to_assess_for_activation_gate=gate_position_multiplication_factor*(define_gating_strategy_Blanc_condition(REP_filename_and_path, HEK_cell_fraction, single_cell_gate_fraction, '640-670/14-A', param_reporter_colour)[1])
 
  #Figuring out HEK293T, single cells, iBlue+ and ECFP+ gates in the REP file
  gating_values_reporter_file=define_gating_strategy_REP_condition_iBlue_ECFP(REP_filename_and_path, HEK_cell_fraction, single_cell_gate_fraction, transfection_gate_cuttoff, initial_value_to_assess_for_activation_gate, max_value_to_assess_for_activation_gate)
  HEK293T_contour=gating_values_reporter_file[0]
  single_cell_contour=gating_values_reporter_file[1]
  iBlue_gate_position=transfection_gate_cuttoff
  ECFP_gate_position=gating_values_reporter_file[8]  
  percentage_activated_ECFP_cells_in_reporter_condition=gating_values_reporter_file[7] 
  
  #contout HEK
  HEK_x_list=[]
  HEK_y_list=[]
  for element in HEK293T_contour:
      HEK_x_list.append(float(element[0]))
      HEK_y_list.append(float(element[1]))
      
  #contout single cells   
  single_cells_x_list=[]
  single_cells_y_list=[]
  for element in single_cell_contour:
      single_cells_x_list.append(float(element[0]))
      single_cells_y_list.append(float(element[1]))
      
  #####output_log_x_values_and_min_and_max_values for HEK293T cells
  HEK_293T_output=output_log_x_values_and_min_and_max_values(HEK293T_contour)
  HEK_min_x_value=HEK_293T_output[0]
  HEK_max_x_value=HEK_293T_output[1]
  HEK_min_y_value=HEK_293T_output[2]
  HEK_max_y_value=HEK_293T_output[3]
  HEK_list_log_x_values=HEK_293T_output[4]
  HEK_list_float_min_y_corresp_x=HEK_293T_output[5]
  HEK_list_float_max_y_corresp_x=HEK_293T_output[6]
  HEK_centre=HEK_293T_output[7]
  HEK_angle=HEK_293T_output[8]
  HEK_axes=HEK_293T_output[9]
  
  if (HEK_centre!=[0,0]) and  (HEK_angle!=0) and (HEK_axes !=[0,0]) and (float(HEK_centre[0])>50000) and (float(HEK_centre[1])>50000):
      global_HEK_centre=HEK_centre
      global_HEK_angle=HEK_angle
      global_HEK_axes=HEK_axes
      print  "Use local HEK cell data"
      
  else: 
      HEK_centre=global_HEK_centre
      HEK_angle=global_HEK_angle
      HEK_axes=global_HEK_axes 
      print  "Use global HEK cell data"


      
  print  "HEK local", HEK_centre, HEK_angle, HEK_axes
  print  "HEK global", global_HEK_centre, global_HEK_angle, global_HEK_axes
  
  #####output_log_x_values_and_min_and_max_values for single cells
  single_cell_output=output_log_x_values_and_min_and_max_values(single_cell_contour)
  single_cell_min_x_value=single_cell_output[0]
  single_cell_max_x_value=single_cell_output[1]
  single_cell_min_y_value=single_cell_output[2]
  single_cell_max_y_value=single_cell_output[3]
  single_cell_list_log_x_values=single_cell_output[4]
  single_cell_list_float_min_y_corresp_x=single_cell_output[5]
  single_cell_list_float_max_y_corresp_x=single_cell_output[6]
  single_cell_centre=single_cell_output[7]
  single_cell_angle=single_cell_output[8]
  single_cell_axes=single_cell_output[9]
  
  if (single_cell_centre!=[0,0]) and  (single_cell_angle!=0) and (single_cell_axes !=[0,0]):
     global_single_cell_centre=single_cell_centre
     global_single_cell_angle=single_cell_angle
     global_single_cell_axes=single_cell_axes
     print "Use local single cell data"
     
  if (single_cell_centre==[0,0]) and  (single_cell_angle==0) and (single_cell_axes ==[0,0]):
      single_cell_centre=global_single_cell_centre
      single_cell_angle=global_single_cell_angle
      single_cell_axes=global_single_cell_axes
      print "Use global single cell data"
            
  print "Single cells local:", single_cell_centre, single_cell_angle, single_cell_axes
  print "Single cells global:",global_single_cell_centre, global_single_cell_angle, global_single_cell_axes
  
  index_individual_file=0
    
  while index_individual_file < len(names_conditions):
    if (final_index_blanc<=index_individual_file) and (index_individual_file<(final_index_blanc+number_occurence_conditions_with_particular_gRNA)):
      file_to_analyse=files_corresponding_replicate[index_individual_file]
      print "File to analyse: ", file_to_analyse           
      print "Transfection condition: ", names_conditions[index_individual_file]
      
      #selecting HEK293T cells
      transformed_input_data = FlowCal.io.FCSData(file_to_analyse)
      transformed_input_data = FlowCal.transform.to_rfi(transformed_input_data)
      HEK_desired_parameter_1='FSC-A' 
      HEK_desired_parameter_2='SSC-A'
      
      #Plotting cells before selection:
      initial_all_data = FlowCal.gate.density2d(transformed_input_data, channels=[HEK_desired_parameter_1, HEK_desired_parameter_2], gate_fraction=1) 
      FlowCal.plot.density2d(initial_all_data, channels=[HEK_desired_parameter_1, HEK_desired_parameter_2], mode='scatter', xscale='linear', yscale='linear')
      plt.plot(HEK_x_list, HEK_y_list, c='black', marker='.', markersize=0.0001)
      plt.show()
      
      selected_HEK293T_cells_output=selecting_cells_from_given_gate_and_input_file(transformed_input_data, HEK_desired_parameter_1, HEK_desired_parameter_2, HEK_min_x_value, HEK_max_x_value, HEK_min_y_value, HEK_max_y_value, HEK_list_log_x_values, HEK_list_float_min_y_corresp_x, HEK_list_float_max_y_corresp_x, HEK_centre,HEK_angle,HEK_axes)
      selected_HEK293T_cells=selected_HEK293T_cells_output[0]
      
      #Plotting cells after selection:
      selected_all_data = FlowCal.gate.density2d(selected_HEK293T_cells, channels=[HEK_desired_parameter_1, HEK_desired_parameter_2], gate_fraction=1) 
      FlowCal.plot.density2d(selected_all_data, channels=[HEK_desired_parameter_1, HEK_desired_parameter_2], mode='scatter', xscale='linear', yscale='linear')
      plt.plot(HEK_x_list, HEK_y_list, c='black', marker='.', markersize=0.0001)
      plt.show()
      
      number_HEK293T_cells=selected_HEK293T_cells_output[1]
      list_number_HEK293T_cells.append(number_HEK293T_cells)
      
      print "No HEK cells:", number_HEK293T_cells
      
      #selecting single cells
            
      single_cell_desired_parameter_1='SSC-A'
      single_cell_desired_parameter_2='SSC-H'
      
      #Plotting cells before selection:
      initial_single_cell_data = FlowCal.gate.density2d(selected_HEK293T_cells, channels=[single_cell_desired_parameter_1, single_cell_desired_parameter_2], gate_fraction=1) 
      FlowCal.plot.density2d(initial_single_cell_data, channels=[single_cell_desired_parameter_1, single_cell_desired_parameter_2], mode='scatter', xscale='linear', yscale='linear')
      plt.plot(single_cells_x_list, single_cells_y_list, c='black', marker='.', markersize=0.0001)
      plt.show()

      selected_single_cells_output=selecting_cells_from_given_gate_and_input_file(selected_HEK293T_cells, single_cell_desired_parameter_1, single_cell_desired_parameter_2, single_cell_min_x_value, single_cell_max_x_value, single_cell_min_y_value, single_cell_max_y_value, single_cell_list_log_x_values, single_cell_list_float_min_y_corresp_x, single_cell_list_float_max_y_corresp_x, single_cell_centre, single_cell_angle, single_cell_axes)
      selected_single_cells=selected_single_cells_output[0]
      
      #Plotting cells after selection:
      selected_single_cell_data = FlowCal.gate.density2d(selected_single_cells, channels=[single_cell_desired_parameter_1, single_cell_desired_parameter_2], gate_fraction=1) 
      FlowCal.plot.density2d(selected_single_cell_data, channels=[single_cell_desired_parameter_1, single_cell_desired_parameter_2], mode='scatter', xscale='linear', yscale='linear')
      plt.plot(single_cells_x_list, single_cells_y_list, c='black', marker='.', markersize=0.0001)
      plt.show()
      
      number_single_cells=selected_single_cells_output[1]
      list_number_single_cells.append(number_single_cells)     
      
      print "No single cells:", number_single_cells
      
      #gating transfected cells
      g_transfected_cells= FlowCal.gate.high_low(selected_single_cells, channels='640-670/14-A', low=iBlue_gate_position)
      number_transfected_cells=len(g_transfected_cells)
      list_number_transfected_cells.append(number_transfected_cells)
      
      #selecting activated cells
      g_activated_cells= FlowCal.gate.high_low(g_transfected_cells, channels=param_reporter_colour, low=ECFP_gate_position)
      number_activated_cells=len(g_activated_cells) 
      if number_transfected_cells ==float(0):
          list_number_percentage_activated_transfected.append(0)
      else: 
          list_number_percentage_activated_transfected.append(float(float(number_activated_cells)/float(number_transfected_cells)*100))
      
      #gating the two fluorophores at single-cell level
      FlowCal.plot.density2d(selected_single_cells, channels=['640-670/14-A', param_reporter_colour], mode='scatter', xscale='log', yscale='log',xlim=(5, 110000), ylim=(5, 110000), xlabel='iBlue transfection control', ylabel='ECFP reporter')
      
      print index_individual_file
      plot_filename=str(full_condition_details_and_indexes[index_individual_file][13])+"_"+str(full_condition_details_and_indexes[index_individual_file][0])+"_"+str(full_condition_details_and_indexes[index_individual_file][1])+'.png'
      plot_name="output_figure_plots/"+plot_filename
      
      #drawing the y gate
      x_values=[]
      y_values=[]
      index=0
      while index <250000:
        x=index
        x_values.append(x)
        y=ECFP_gate_position
        y_values.append(y)
        index=index+200
      plt.plot(x_values, y_values, c='black', marker='.', markersize=0.0001)

      #drawing the x gate
      x2_values=[]
      y2_values=[]
      index2=0
      while index2 <250000:
        y2_values.append(index2)
        x2_values.append(iBlue_gate_position)
        index2=index2+200
      plt.plot(x2_values, y2_values, c='black', marker='.', markersize=0.0001)
      plt.savefig(plot_name, dpi=300, bbox_inches='tight')
      plt.show()
      
      
      #Extracting information about ECFP intensity
      s_transformed = FlowCal.transform.to_rfi(g_transfected_cells, channels=param_reporter_colour)
      ECFP_positive_cells=s_transformed[:,[param_reporter_colour]]
      #Change for mean rather then median
      median_ECFP_intensity_transfected_cells=np.median(ECFP_positive_cells)
      list_median_ECFP_transfected_cells.append(median_ECFP_intensity_transfected_cells)
      
      print 'NO HEK293T cells: ', number_HEK293T_cells
      print 'NO single cells: ', number_single_cells
      print 'NO transfected cells: ', number_transfected_cells
      if number_transfected_cells==0:
         print '%age activated transfected cells: ', 0
      else:
         print '%age activated transfected cells: ', float(float(number_activated_cells)/float(number_transfected_cells)*100)
      print 'Median ECFP: ', median_ECFP_intensity_transfected_cells
         
    index_individual_file=index_individual_file+1
  
  index_guide_number=index_guide_number+1
  
 list_number_HEK293T_cells_all_replicates.append(list_number_HEK293T_cells)
 list_number_single_cells_all_replicates.append(list_number_single_cells)
 list_number_transfected_cells_all_replicates.append(list_number_transfected_cells)
 list_number_percentage_activated_transfected_all_replicates.append(list_number_percentage_activated_transfected)
 list_median_ECFP_transfected_cells_all_replicates.append(list_median_ECFP_transfected_cells)
  
 index_number_replicates=index_number_replicates+1

output_line_1='@Index_condition,Condition_name,Plasmid_1, Plasmid_2, Plasmid_3, Plasmid_4, Plasmid_5, Plasmid_6, Plasmid_7, filename_repeat_1, filename_repeat_2, filename_repeat_3, Replicate_type, transfection_date_repeat_1,transfection_date_repeat_2, transfection_date_repeat_3, gRNA_colour, percentage_activation_R1, percentage_activation_R2, percentage_activation_R3, median_ECFP_R1, median_ECFP_R2, median_ECFP_R3, average_percentage_activated_cells,std_percentage_activated_cells,average_med_ECFP,std_med_ECFP,Norm_SPA_average_percentage_activated_cells,Norm_SPA_std_percentage_activated_cells,Norm_SPA_average_med_ECFP,Norm_SPA_std_med_ECFP,ON/OFF_ratios_percentage_activation,ON/OFF_ratios_median_ECFP,,,,\n' 
intermediary_output_file.write(output_line_1)

index=0
while index < len(full_condition_details_and_indexes):
    print 'index',index
    condition_name=full_condition_details_and_indexes[index]
    gRNA_colour=gRNA_colours[index]
    if total_number_replicates ==1:
      percentage_activation_R1=list_number_percentage_activated_transfected_all_replicates[0][index]
      median_ECFP_R1=list_median_ECFP_transfected_cells_all_replicates[0][index]
      percentage_activation_R2='NA'
      median_ECFP_R2='NA'
      percentage_activation_R3='NA'
      median_ECFP_R3='NA'
    
    if total_number_replicates ==2:
      percentage_activation_R1=list_number_percentage_activated_transfected_all_replicates[0][index]
      median_ECFP_R1=list_median_ECFP_transfected_cells_all_replicates[0][index]
      percentage_activation_R2=list_number_percentage_activated_transfected_all_replicates[1][index]
      median_ECFP_R2=list_median_ECFP_transfected_cells_all_replicates[1][index]
      percentage_activation_R3='NA'
      median_ECFP_R3='NA'
    
    if total_number_replicates==3:
      percentage_activation_R1=list_number_percentage_activated_transfected_all_replicates[0][index]
      median_ECFP_R1=list_median_ECFP_transfected_cells_all_replicates[0][index]
      percentage_activation_R2=list_number_percentage_activated_transfected_all_replicates[1][index]
      median_ECFP_R2=list_median_ECFP_transfected_cells_all_replicates[1][index]
      percentage_activation_R3=list_number_percentage_activated_transfected_all_replicates[2][index]
      median_ECFP_R3=list_median_ECFP_transfected_cells_all_replicates[2][index]

    output=condition_name[0]+','+condition_name[1]+','+condition_name[2]+','+condition_name[3]+','+condition_name[4]+','+condition_name[5]+','+condition_name[6]+','+condition_name[7]+','+condition_name[8]+','+condition_name[9]+','+condition_name[10]+','+condition_name[11]+','+condition_name[12]+','+condition_name[13]+','+condition_name[14]+','+condition_name[15]+','+gRNA_colour+','+str(percentage_activation_R1)+","+str(percentage_activation_R2)+","+str(percentage_activation_R3)+","+str(median_ECFP_R1)+","+str(median_ECFP_R2)+","+str(median_ECFP_R3)+",,,\n"
    intermediary_output_file.write(output)
    
    index=index+1
    

intermediary_output_file.close()


normalisation_vector_percentage_activation=[]
normalisation_vector_median_ECFP=[]  
OFF_average_percentage_activation_value=0
OFF_average_median_ECFP=0
spa_index=0

for line in open(intermediary_filename).readlines():
  #print line
  if "percentage_activation_R1" in line:
      output=str(line)
      output_file_individual_analysis.write(output)
      output_file_all_data_analysis.write(output)
  elif ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," in line:
      print "do nothing" 
  else:
    f=line.split(',')
    info_to_maintain=str(f[0:23])[1:-1]
    percentage_activation_replicates=f[17:20]
    median_ECFP_replicates=f[20:23]

    list_replicates_percentage_activation_replicates=[]
    for element in percentage_activation_replicates:
        if element!='NA':
            list_replicates_percentage_activation_replicates.append(float(element))
    
    list_replicates_median_ECFP_replicates=[]
    for element2 in median_ECFP_replicates:
        if element2!='NA':
            list_replicates_median_ECFP_replicates.append(float(element2))
    
    array_replicates_percentage_activation=np.array(list_replicates_percentage_activation_replicates)
    array_median_ECFP_replicates=np.array(list_replicates_median_ECFP_replicates)
    
    ########determine first 4 important parameters
    average_percentage_activation=np.mean(array_replicates_percentage_activation)
    stdev_percentage_activation=np.std(array_replicates_percentage_activation)
    
    average_median_ECFP=np.mean(array_median_ECFP_replicates)
    stdev_median_ECFP=np.std(array_median_ECFP_replicates)
    
#    print  average_percentage_activation, stdev_percentage_activation, average_median_ECFP, stdev_median_ECFP  
    
    #determine spacer and rep condition
    vector_percentage_activation=[]
    vector_median_ECFP=[]
    
    if (f[1]=='REP') or (f[1]=='Blanc'):
        spa_index=0
        for element in list_replicates_median_ECFP_replicates:
            vector_percentage_activation.append(1)
            vector_median_ECFP.append(1)
        normalisation_vector_percentage_activation=np.array(vector_percentage_activation)
        normalisation_vector_median_ECFP=np.array(vector_median_ECFP)
        
    if ('SPA' in line) and (f[1]=='SPA'):
            normalisation_vector_percentage_activation=array_replicates_percentage_activation
            normalisation_vector_median_ECFP=array_median_ECFP_replicates
            spa_index=1
    print f[1], spa_index
    #determining the other parameters
    if (('REP' in line) or ('Blanc' in line)) and (spa_index==0):
       norm_average_percentage_activation='NA'
       norm_stdev_percentage_activation='NA'
    
       norm_average_median_ECFP='NA'
       norm_stdev_median_ECFP='NA' 
    
    else:  
      array_normalised_percentage_activated_cells=(array_replicates_percentage_activation*100)/normalisation_vector_percentage_activation
      array_normalised_median_ECFP=(array_median_ECFP_replicates*100)/normalisation_vector_median_ECFP
   
      norm_average_percentage_activation=np.mean(array_normalised_percentage_activated_cells)
      norm_stdev_percentage_activation=np.std(array_normalised_percentage_activated_cells)
    
      norm_average_median_ECFP=np.mean(array_normalised_median_ECFP)
      norm_stdev_median_ECFP=np.std(array_normalised_median_ECFP) 


    output_2=str(average_percentage_activation)+","+str(stdev_percentage_activation)+","+str(average_median_ECFP)+","+str(stdev_median_ECFP)
    output_3=str(norm_average_percentage_activation)+","+str(norm_stdev_percentage_activation)+","+str(norm_average_median_ECFP)+","+str(norm_stdev_median_ECFP)
   
    if ('OFF' in line) and f[1]=="OFF":
      OFF_average_percentage_activation_value=average_percentage_activation
      OFF_average_median_ECFP=average_median_ECFP
      ON_OFF_ratio_percentage_activation='NA'
      ON_OFF_ratio_median_ECFP='NA'
    elif 'ON' in line and f[1]=="ON":
      ON_OFF_ratio_percentage_activation=float(average_percentage_activation/OFF_average_percentage_activation_value)
      ON_OFF_ratio_median_ECFP=float(average_median_ECFP/OFF_average_median_ECFP)
    else:
      ON_OFF_ratio_percentage_activation='NA'
      ON_OFF_ratio_median_ECFP='NA'
      
    output_4=str(ON_OFF_ratio_percentage_activation)+","+str(ON_OFF_ratio_median_ECFP)+',,,,\n'
    
    output_to_print=(info_to_maintain+','+output_2+','+output_3+','+output_4)
    output_file_individual_analysis.write(output_to_print)
    output_file_all_data_analysis.write(output_to_print)
    
output_file_individual_analysis.close()
output_file_all_data_analysis.close()
