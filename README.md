# Learning experiment
> Syngjoo Choi, Sanjeev Goyal, Frederic Moisan, Yu Yang Tony To

The objective of this project is to understand the effects of network types on information aggregation and social learning. <br />
This repository includes all the coding done for the project, including Inputs, DeGroot simulations, Experiment analysis, Regressions, Learning-rule analysis and Outputs. <br />

## How to use
0. Ensure you can run R
1. Download all the documents from "R", "input" folder
2. Begin with running 00->05 .R files in order
   - "00-Simulation.R" takes input from "networks", "signals", "functions" and output simulation results "Simulation_xxx.pdf"
   - "01-process_data.R" takes input from "raw_data", "00-simulation-envir.RData", and output cleaned data summaries "01-process_data-envir.RData"
   - "02-analysis.R" takes prev. as input, to perform analysis and output important dataframes "02-analysis-envir.RData"
   - "03-analysis_output.R" takes prev. as input, to create output of analysis into "output"
   - "04-regression.R" takes input from "02-analysis-envir.RData" to run regressions, and output "Reg_output_xxx.txt"
   - "05-degroot.R" takes input from "02-analysis-envir.RData" to do learning rule analysis, and output into "output"
3. Read "output" folder for graphs, tables, csv

## Content
"R" folder includes:
- R files:  
- "RData" folder: workspace saved from R files
- "function" folder: functions used for simulation and measuring network statistics

"input" folder includes:
- "networks" folder: csv files of network graphs, images of networks
- "raw_data" folder: May2021 experimental data
- "signals" folder: 24 sets of signals used

"output" folder includes:
- "dataframe" folder: dataframe outputted from R files
- "latex" folder: tables for printing to LaTeX file "Output_prints.tex" 
- outputs from simulations
- outputs from analysis
- outputs from regressions
