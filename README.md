# HCVHIVWithinHostModel
This Repository contains the code for our within-host model of HCV-HIV dynamics

There 3 scripts, which call various functions as described below

####Figures_2_4_6_Script
This script creates and plots Figures 2, 4, and 6 from our Epidemics paper. It uses a data file, SnoeckData_AllOutcomes.txt, and calls the following other functions and scripts: paramsets, VLplot, VLplotIC, VLplotDecay, SimulateTreatmentDurations, SimulateTreatmentDurationsDecay, adjust_delta

####Figure3script
This script creates and plots Figure 3 from our Epidemics paper. It calls paramsets and bifurc_ThIc.

####Figure5script
This script creates and plots Figure 5 from our Epidemics paper. It calls paramsets and TrtEffvsDurPlot.

####paramsets
This script populates all the parameter sets used in creating figures 2-6.

####VLplot, VLplotIC, VLplotDecay
These functions take as input symbolic equilibrium formulae and parameter sets, and generate initial population vectors. Using these vectors, they run the model and get numerical solutions for the ODEs. They output time and population vectors. VLplot using the standard version of the model, VLplotIC uses the QSS approximation, and VLPlotDecay uses the standard version with treatment decay. They call, respectively, the functions hcvhivThIcVc, hcvhivThIc, hcvhivThIcVcDecay









