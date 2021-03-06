### Resolving network communities in the human brain
-----------------------------------------------
#### what's in this repo?
This repository contains simulation and analysis code for the project "Improving Resolution of Dynamic Communities in Human Brain Networks through Targeted Node Removal." 

The goal is to improve the resolution of network communities and reveal information about the differences between dynamic network states, by removing targeted nodes from the network both before and after community detection and comparing the results.

The project tests and evaluates the technique on synthetic Kuramoto networks of oscillators with a known underlying communities. Then it applies the technique to networks of activity in the human brain, revealing new ways to compare data sets across different cognitive states.

#### organization
The main analysis files in the repository include:

 - mini_fig.m : prototype of small network schematic (Figs. 1 and 2)  
 - kuramoto_examples.m : example code for running kuramoto simulations and saving results  
 - kuramoto_statplots.m : script to create figures and test hypotheses from saved kuramoto oscillator examples (Figs. 3-5)  
 - brain_cd_analysis.m : script to compare brain network flexibility before and after removal of visual cortex (Figs. 6, 7, and 9)  
 - power_region_analysis.m : script for computing differences in brain functional system recruitment before and after removal of visual cortex (Figs. 8 and S1)

Other files include:

 - directories of utility functions for...  
     - community detection (CD_Utils/)  
     - running kuramoto simulations (Kuramoto/)  
     - plotting and statistical analysis (Analysis_Utils/)
     - brain region labeling (Region_Info/)  
     - basic MATLAB functionality (MATLAB_Utils/)  
