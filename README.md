# Flow_Cytometry_Data_Analysis

Flow_Cytometry_Data_Analysis is a pipeline for analysing flow cytometry .fcs files to quantify reporter activation in HEK293T cells. It applies gating strategies, calculates activation metrics, and generates outputs for further analysis. This code is associated with the publication:

# Specific Modulation of CRISPR Transcriptional Activators through RNA-Sensing Guide RNAs in Mammalian Cells and Zebrafish Embryos
Oana Pelea, Sarah Mayes, Quentin R.V. Ferry, Tudor A. Fulga, Tatjana Sauka-Spengler
eLife. https://doi.org/10.7554/eLife.87722.2

Scripts

# Analysing_FCS_files.py
Main pipeline script. Reads input .csv files specifying conditions and replicates, processes corresponding .fcs data, applies gating strategies, and outputs numerical summaries and plots.

# functions_fcs_file_analysis.py

Contains functions to:

    Define gating strategies for Blanc and REP conditions
    Select cells within gates
    Extract contour data and gate coordinates

# functions_FlowC_analysis.py

Contains functions to:

    Assign guide colour codes
    Generate bar plots from output data

# Example input files

Folder input_FCS_files_and_lists/0_Files_inputs_for_fcs_analysis_program includes an example input .csv file (00_example_data.csv) together with example .fcs files.

Please use the same format for your own experiments. The input file has three columns for replicate filenames (filename_repeat_1, filename_repeat_2, and filename_repeat_3). Fill in the appropriate data names and file paths for up to three experimental replicates per condition – each column can take the input filename for one replicate.

Optionally, you can also modify condition names as well as the exact plasmid names transfected in each experiment to match your specific setup.

# Prerequisites

Requirement	Version	Notes
Python	2.7.x	Code won’t run under Python ≥3.0 without edits.
FlowCal	–	For flow cytometry data parsing and plotting
numpy	–	
scipy	–	
pandas	–	
matplotlib	–	
scikit-image	–	
fcsparser	–	
