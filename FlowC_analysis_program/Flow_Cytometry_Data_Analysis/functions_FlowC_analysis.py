#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:47:41 2019

@author: oanapelea
"""
import numpy as np
# import matplotlib.pyplot as plt
import matplotlib  
matplotlib.use('TkAgg')   
import matplotlib.pyplot as plt  

output_subdirectory="output_figure_barrgraphs/"

#def changing_guide_names_to_standard(guide_name):
#    new_guide_name=guide_name
#    if guide_name=='P1':
#        new_guide_name='G1'
#    if guide_name=='P2':
#        new_guide_name='G2'
#    if guide_name=='G2':
#        new_guide_name='G4'
#    if guide_name=='G3':
#        new_guide_name='G3'
#    if guide_name=='G5':
#        new_guide_name='G5'              
#    return(new_guide_name)
    
def guide_colour_code(guide_name):
    colour_guide='black'
    if ("G1" == guide_name) or ("G1c" == guide_name):
        colour_guide='paleturquoise'
    if ("G2" == guide_name) or ("G2c" == guide_name):
        colour_guide='rebeccapurple'   
    if ("G3" == guide_name) or ("G3c" == guide_name):
        colour_guide='navy'
    if ("G4" == guide_name) or ("G4c" == guide_name):
        colour_guide='blue'
    if ("G5" == guide_name) or ("G5c" == guide_name):
        colour_guide='mediumorchid'
    if ("P1" == guide_name) or ("P1c" == guide_name):
        colour_guide='lightskyblue'
    if ("P2" == guide_name) or ("P2cc" == guide_name):
        colour_guide='dodgerblue'  
        
    if ("E2" == guide_name) or ("E2c" == guide_name):
        colour_guide='black'
    if ("E3" == guide_name) or ("E3c" == guide_name):
        colour_guide='teal'
    if ("E4" == guide_name) or ("E4c" == guide_name):
        colour_guide='cyan'
    if ("E6" == guide_name) or ("E6c" == guide_name):
        colour_guide='darkturquoise'
    if ("E7" == guide_name) or ("E7c" == guide_name):
        colour_guide='cadetblue'
    if ("E8" == guide_name) or ("E8c" == guide_name):
        colour_guide='powderblue'
        
    if ("T4a" == guide_name):
        colour_guide='deeppink'
    if ("T4b" == guide_name):
        colour_guide='fuchsia'
    if ("T4c" == guide_name):
        colour_guide='violet'
        
    return(colour_guide)
  
def plot_graph(plot_title, list_condition_names, y_axis_title, list_average_values, list_stDEVs, list_guide_colours, plot_y_max_lim):
    
    final_plot_title_summary=''
    output_summary_image_name=plot_title+".svg"

    x_pos = np.arange(len(list_condition_names))
    names=list_condition_names
    CTEs = list_average_values
    error = list_stDEVs
    desired_colors=list_guide_colours
    rotation_angle=55
    
    number_conditions=len(list_average_values)
    if number_conditions <=6:
        size=8
       
    if (number_conditions > 5) and (number_conditions <10):
        size=6

    if number_conditions >= 10:
        size=4

    
    fig, ax = plt.subplots()
    ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=1, color=desired_colors, capsize=float(size))
    ax.set_ylabel(y_axis_title)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, rotation=rotation_angle)
    ax.set_title(final_plot_title_summary)
    #ax.yaxis.grid(color='black',linestyle='--', linewidth=0.1)

    plt.tight_layout()
    plt.ylim(0, plot_y_max_lim)

    plt.savefig(output_subdirectory+output_summary_image_name, format='svg', dpi=1200) 
    
    plt.show()
    
    output_message="Check your plot:"+ output_subdirectory+output_summary_image_name +" \n"
    return(output_message)

def figure_out_intermediary_output_file_data(given_parameter_to_plot):  
    if given_parameter_to_plot=="median CFP intensity":
        desired_file="intermediary_outputs/VII_unprocessed_data_median_CFP.csv"
    if given_parameter_to_plot=="normalised median CFP intensity":
        desired_file="intermediary_outputs/IV_normalised_data_median_CFP.csv"
    if given_parameter_to_plot=="proportion activated cells":
        desired_file="intermediary_outputs/IV_normalised_data_median_CFP.csv"
    if given_parameter_to_plot=="normalised proportion activated cells":
        desired_file="intermediary_outputs/V_normalised_data_proportion_activated_cells.csv"
    if given_parameter_to_plot=="activation score":
        desired_file="intermediary_outputs/IX_unprocessed_data_activation_score.csv"
    if given_parameter_to_plot=="normalised activation score":
        desired_file="intermediary_outputs/VI_normalised_data_activation_score.csv"
    return(desired_file)


