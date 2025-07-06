#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:49:27 2020

@author: oanapelea
"""

import FlowCal
import pylab
import fcsparser
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import math
import pandas as pd
from matplotlib.patches import Ellipse
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

def define_gating_strategy_Blanc_condition(Blanc_filename_and_path, HEK_cell_fraction, single_cell_gate_fraction, fluorophore_1, fluorophore_2):
  
    input_file_path=Blanc_filename_and_path
    
    #plotting all data
    all_data = FlowCal.io.FCSData(input_file_path)
    all_data = FlowCal.transform.to_rfi(all_data)
    
    #gating HEK293T cells gate
    g_HEK293T = FlowCal.gate.density2d(all_data, channels=['FSC-A', 'SSC-A'], gate_fraction=HEK_cell_fraction, full_output=True)

    #gating  single cells
    g_single_cell_HEK293T=FlowCal.gate.density2d(g_HEK293T[0], channels=['SSC-A', 'SSC-H'], gate_fraction=single_cell_gate_fraction, full_output=True)

    #determining max fluorescence intensity for fluorophore 1 and fluorophore 2
    untransfected_HEK293T =FlowCal.gate.density2d(g_single_cell_HEK293T[0], channels=[fluorophore_1, fluorophore_2], gate_fraction=0.99, full_output=True )
    
    contour_blanc=untransfected_HEK293T[2][0]
    list_values_fluorophore_1=[]
    list_values_fluorophore_2=[]
    for element in contour_blanc:
      if float(element[0])>0:
        list_values_fluorophore_1.append(float(element[0]))
      if float(element[1])>0:
        list_values_fluorophore_2.append(float(element[1]))
        
    transfection_gate_cuttoff_fluorophore_1=max(list_values_fluorophore_1)
    transfection_gate_cuttoff_fluorophore_2=max(list_values_fluorophore_2)
        
    return(transfection_gate_cuttoff_fluorophore_1, transfection_gate_cuttoff_fluorophore_2)
   
def define_gating_strategy_REP_condition_iBlue_ECFP(REP_filename_and_path, HEK_cell_fraction, single_cell_gate_fraction, transfection_gate_cuttoff, initial_value_to_assess_for_activation_gate, max_value_to_assess_for_activation_gate):
    input_file_path=REP_filename_and_path
    
    max_desired_activation_gate_cutt_off=float(max_value_to_assess_for_activation_gate)
    print initial_value_to_assess_for_activation_gate, max_desired_activation_gate_cutt_off
    #plotting all data
    all_data = FlowCal.io.FCSData(input_file_path)
    all_data = FlowCal.transform.to_rfi(all_data)
    number_measured_events=len(all_data)
    
    #gating HEK293T cells gate
    g_HEK293Ti = FlowCal.gate.density2d(all_data, channels=['FSC-A', 'SSC-A'], gate_fraction=1)
    g_HEK293T = FlowCal.gate.density2d(all_data, channels=['FSC-A', 'SSC-A'], gate_fraction=HEK_cell_fraction, full_output=True)
    FlowCal.plot.density2d(g_HEK293Ti, channels=['FSC-A', 'SSC-A'], mode='scatter', xscale='linear', yscale='linear')

    list_x_val_contour_coordinates=[]
    list_y_val_contour_coordinates=[]
    contour=g_HEK293T[2][0]
    for element in contour:
        list_x_val_contour_coordinates.append(float(element[0]))
        list_y_val_contour_coordinates.append(float(element[1]))
    plt.plot(list_x_val_contour_coordinates, list_y_val_contour_coordinates, c='black', marker='.', markersize=0.0001)
    plt.show()

    number_HEK293T_cells=len(g_HEK293T[0])
    
    HEK_cell_contour_details=contour #important output
    
    #gating  single cells
    single_cell_HEK293T=FlowCal.gate.density2d(g_HEK293T[0], channels=['SSC-A', 'SSC-H'], gate_fraction=1)
    FlowCal.plot.density2d(single_cell_HEK293T, channels=['SSC-A', 'SSC-H'], mode='scatter', xscale='linear', yscale='linear')

    g_single_cell_HEK293T=FlowCal.gate.density2d(g_HEK293T[0], channels=['SSC-A', 'SSC-H'], gate_fraction=single_cell_gate_fraction, full_output=True)
    contour_single_cells=g_single_cell_HEK293T[2][0]

    list_x_val_contour_coordinates_single_cells=[]
    list_y_val_contour_coordinates_single_cells=[]

    for element2 in contour_single_cells:
        list_x_val_contour_coordinates_single_cells.append(float(element2[0]))
        list_y_val_contour_coordinates_single_cells.append(float(element2[1]))
    plt.plot(list_x_val_contour_coordinates_single_cells, list_y_val_contour_coordinates_single_cells, c='black', marker='.', markersize=0.0001)
    plt.show()

    number_single_cells=len(g_single_cell_HEK293T[0])
    single_cell_contour_details= contour_single_cells#important output
    
    #gating transfected cells
    cuttoff_transfected_cells=transfection_gate_cuttoff

    #gating transfected cells
    g_transfected_cells= FlowCal.gate.high_low(g_single_cell_HEK293T[0], channels='640-670/14-A', low=cuttoff_transfected_cells)

    number_transfected_cells=len(g_transfected_cells)
    
    #gating_activated_cells

    desired_interval_fraction_activated_cells_in_reporter=[float(0.10000),float(0.1500)]
    initial_value_of_the_gate=initial_value_to_assess_for_activation_gate

    intermediate_g_activated_cells= FlowCal.gate.high_low(g_transfected_cells, channels='405-450/50-A', low=initial_value_of_the_gate)
    FlowCal.plot.density2d(intermediate_g_activated_cells,channels=['640-670/14-A', '405-450/50-A'], mode='scatter', xscale='log', yscale='log')


    intermediate_number_activated_cells=len(intermediate_g_activated_cells)
    intermediate_fraction_activated_cells_in_the_reporter=100*float(intermediate_number_activated_cells)/float(number_transfected_cells)

    #decrease gate position case
    if (intermediate_fraction_activated_cells_in_the_reporter < desired_interval_fraction_activated_cells_in_reporter[0]):
        index4=initial_value_of_the_gate
        desired_activated_cells_gate_position_comment="will need to decrease gate position"
        print desired_activated_cells_gate_position_comment
        intermediate2_fraction_activated_cells_in_the_reporter=intermediate_fraction_activated_cells_in_the_reporter
        while index4 > 0.0:
            new_gate_position2=float(initial_value_of_the_gate-index4)
            intermediate2_g_activated_cells= FlowCal.gate.high_low(g_transfected_cells, channels='405-450/50-A', low=new_gate_position2)
            number_activated_cells_int2=float(len(intermediate2_g_activated_cells))
            intermediate2_fraction_activated_cells_in_the_reporter=float(100*(float(number_activated_cells_int2)/float(number_transfected_cells)))
            if (intermediate2_fraction_activated_cells_in_the_reporter >= desired_interval_fraction_activated_cells_in_reporter[0]) and (intermediate2_fraction_activated_cells_in_the_reporter <= desired_interval_fraction_activated_cells_in_reporter[1]):
               desired_ECFP_gate_position= new_gate_position2   
               number_activated_cells=number_activated_cells_int2
               fraction_activated_cells=float(intermediate2_fraction_activated_cells_in_the_reporter)
            index4=index4-0.1
            
    #equality case
    if (intermediate_fraction_activated_cells_in_the_reporter >=desired_interval_fraction_activated_cells_in_reporter[0]) and (intermediate_fraction_activated_cells_in_the_reporter <=desired_interval_fraction_activated_cells_in_reporter[1]):
        g_activated_cells=intermediate_g_activated_cells
        number_activated_cells=len( g_activated_cells)
        desired_activated_cells_gate_position_comment='keep gate as it is'
        print desired_activated_cells_gate_position_comment
        fraction_activated_cells=intermediate_fraction_activated_cells_in_the_reporter
        desired_ECFP_gate_position=initial_value_of_the_gate

    #increase gate position case
    if (intermediate_fraction_activated_cells_in_the_reporter > desired_interval_fraction_activated_cells_in_reporter[1]):
        new_gate_position3=initial_value_of_the_gate
        desired_activated_cells_gate_position_comment="will need to increase gate position"
        print desired_activated_cells_gate_position_comment
        intermediate3_fraction_activated_cells_in_the_reporter=intermediate_fraction_activated_cells_in_the_reporter
        while new_gate_position3 < max_desired_activation_gate_cutt_off:
            new_gate_position3=new_gate_position3+0.1
            intermediate3_g_activated_cells= FlowCal.gate.high_low(g_transfected_cells, channels='405-450/50-A', low=new_gate_position3)
            number_activated_cells_int3=float(len(intermediate3_g_activated_cells))
            intermediate3_fraction_activated_cells_in_the_reporter=float(100*(float(number_activated_cells_int3)/float(number_transfected_cells)))

            if (intermediate3_fraction_activated_cells_in_the_reporter >= desired_interval_fraction_activated_cells_in_reporter[0]) and (intermediate3_fraction_activated_cells_in_the_reporter <= desired_interval_fraction_activated_cells_in_reporter[1]):
               desired_ECFP_gate_position=new_gate_position3   
               number_activated_cells=number_activated_cells_int3
               fraction_activated_cells=float(intermediate3_fraction_activated_cells_in_the_reporter)
            new_gate_position3=new_gate_position3+0.5


    #gating the two fluorophores at single-cell level
    FlowCal.plot.density2d(g_single_cell_HEK293T[0], channels=['640-670/14-A', '405-450/50-A'], mode='scatter', xscale='log', yscale='log',xlim=(8, 110000), ylim=(7, 110000), xlabel='iBlue transfection control', ylabel='ECFP reporter')

    #drawing the y gate
    x_values=[]
    y_values=[]
    index=0
    while index <250000:
        x=index
        x_values.append(x)
        y=desired_ECFP_gate_position
        y_values.append(y)
        index=index+1
    plt.plot(x_values, y_values, c='black', marker='.', markersize=0.0001)


    #drawing the x gate
    x2_values=[]
    y2_values=[]
    index2=0
    while index2 <250000:
        y2_values.append(index2)
        x2_values.append(cuttoff_transfected_cells)
        index2=index2+1
    plt.plot(x2_values, y2_values, c='black', marker='.', markersize=0.0001)
    plt.show()

    #Extracting information about ECFP intensity
    s_transformed = FlowCal.transform.to_rfi(g_transfected_cells, channels='405-450/50-A')

    ECFP_positive_cells=s_transformed[:,['405-450/50-A']]
    median_ECFP_intensity_transfected_cells=np.median(ECFP_positive_cells)
    
    return(HEK_cell_contour_details, single_cell_contour_details,number_measured_events ,number_HEK293T_cells, number_single_cells, number_transfected_cells, number_activated_cells, fraction_activated_cells,  desired_ECFP_gate_position, median_ECFP_intensity_transfected_cells)

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])
    

def output_log_x_values_and_min_and_max_values(gate_contour):
    list_x_values=[]
    list_y_values=[]
    unique_list_x_log_values=[]
    unique_list_y_log_values=[]
    my_points=[]
    for point in gate_contour:
        x_value=int(point[0])
        y_value=int(point[1])
        list_x_values.append(x_value) 
        list_y_values.append(y_value)  
        my_points.append([x_value, y_value])
        
        log_x_value=float(format(np.log(x_value), '.2f'))
        log_y_value=float(format(np.log(y_value), '.2f'))
        if log_x_value not in unique_list_x_log_values:
          unique_list_x_log_values.append(log_x_value)
        if log_y_value not in unique_list_y_log_values:
          unique_list_y_log_values.append(log_y_value)        
        
    overall_min_x_value=min(list_x_values)
    overall_max_x_value=max(list_x_values)
    overall_min_y_value=min(list_y_values)
    overall_max_y_value=max(list_y_values)
    
    xc=0
    yc=0
    a=0
    b=0
    theta=0

    a_points = np.array(my_points)      
    ell = EllipseModel()
    ell.estimate(a_points)
    
    if ell.estimate(a_points) ==True:
      xc, yc, a, b, theta = ell.params

      print("center = ",  [xc, yc])
      print("angle of rotation = ",  theta)
      print("axes = ", (a,b))
    
    centre=[xc, yc]
    angle= theta
    axes=[a,b]      
      
    list_log_x_values=[]
    index=min(unique_list_x_log_values)
    while index <= max(unique_list_x_log_values):
        x_value_to_consider=index
        list_log_x_values.append(float(format((x_value_to_consider), '.2f')))
        index=index+0.01
        
    list_log_x_with_documented_y_values=[]
    list_min_documented_y=[]
    list_max_documented_y=[]
    for unique_log_x_unique_value in list_log_x_values:
        list_y_values_given_x=[]
        for point in gate_contour:         
          x_value=int(point[0])
          log_x_value=float(format(np.log(x_value), '.2f'))
          y_value=int(point[1])
          
          if unique_log_x_unique_value == log_x_value:
              list_y_values_given_x.append(y_value) 
              
        if len(list_y_values_given_x)   >0:
          y_min=min(list_y_values_given_x)
          y_max=max(list_y_values_given_x)
          list_log_x_with_documented_y_values.append(unique_log_x_unique_value)
          list_min_documented_y.append(y_min)
          list_max_documented_y.append(y_max)
   
    list_float_min_y_corresp_x=[]
    list_float_max_y_corresp_x=[]                                    
    index=0
    while index < len(list_log_x_values):
        x_value_of_choice=list_log_x_values[index]
        if x_value_of_choice in list_log_x_with_documented_y_values:
            x_val_index=list_log_x_with_documented_y_values.index(x_value_of_choice)
            min_y_value=list_min_documented_y[x_val_index]
            max_y_value=list_max_documented_y[x_val_index]
            
        if x_value_of_choice not in list_log_x_with_documented_y_values:
          if (1<=index) and (index<=(len(list_log_x_values)-1)):
              previous_x_value=list_log_x_values[index-1]
              next_x_value=list_log_x_values[index+1]
              min_list=[]
              max_list=[]
              if previous_x_value in list_log_x_with_documented_y_values:
                 x_prev_val_index=list_log_x_with_documented_y_values.index(previous_x_value)
                 min_y_value=list_min_documented_y[x_prev_val_index]
                 max_y_value=list_max_documented_y[x_prev_val_index]
                 min_list.append(min_y_value)
                 max_list.append(max_y_value)
              if next_x_value in list_log_x_with_documented_y_values:
                 x_next_val_index=list_log_x_with_documented_y_values.index(next_x_value)
                 min_y_value=list_min_documented_y[x_next_val_index]
                 max_y_value=list_max_documented_y[x_next_val_index]
                 min_list.append(min_y_value)
                 max_list.append(max_y_value)
                 
              min_y_value=np.mean(min_list)
              max_y_value=np.mean(max_list)
              
          else:
              min_y_value=0
              max_y_value=0
              
        list_float_min_y_corresp_x.append(min_y_value)
        list_float_max_y_corresp_x.append(max_y_value)
        index=index+1
    
    return(overall_min_x_value, overall_max_x_value, overall_min_y_value, overall_max_y_value, list_log_x_values, list_float_min_y_corresp_x, list_float_max_y_corresp_x, centre, angle, axes)

def selecting_cells_from_given_gate_and_input_file(transformed_input_data, desired_parameter_1, desired_parameter_2, min_x_value, max_x_value, min_y_value, max_y_value, list_log_x_values, list_float_min_y_corresp_x, list_float_max_y_corresp_x, centre, angle, axes):

    #plotting_all_data
    all_data = transformed_input_data
      
    x_truncated_data=FlowCal.gate.high_low(all_data, channels=[desired_parameter_1], high=max_x_value, low=min_x_value, full_output=False)
    y_truncated_data=FlowCal.gate.high_low(x_truncated_data, channels=[desired_parameter_2], high=max_y_value, low=min_y_value, full_output=False)
    selected_data=FlowCal.gate.ellipse(y_truncated_data, channels=[desired_parameter_1, desired_parameter_2],  center=centre, a=axes[0], b=axes[1], theta=angle, log=False, full_output=False)

    return(selected_data, len(selected_data))
    