# HCVHIVWithinHostModel
This Repository contains the code for our within-host model of HCV-HIV dynamics

Functions and scripts are described below

####adjust_delta
This function adjusts delta such that delta/(1+alpha*CD4) = delta, in order to make sure that the total clearance term is within a realistic parameter range. It takes as input a parameter set and the symbolic equilibrium formula for CD4 count.

####clear_test
This function assesses whether SVR has been achieved for a given treatment efficacy and duration. It takes as input a parameter set, a treatment duration, and a treatment efficacy, and runs the model using VLplot. It outputs 1 if SVR is achieved (i.e. if viral load is below 1 at the end of the run) and -1 if it is not.

####Figures_2_4_6_Script
This script creates and plots Figures 2, 4, and 6 from our Epidemics paper. It uses a data file, SnoeckData_AllOutcomes.txt, and calls the following other functions and scripts: paramsets, VLplot, VLplotIC, VLplotDecay, SimulateTreatmentDurations, SimulateTreatmentDurationsDecay, adjust_delta.

####Figure3script
This script creates and plots Figure 3 from our Epidemics paper. It calls paramsets, bifurc_ThIc and adjust_delta.

####Figure5script
This script creates and plots Figure 5 from our Epidemics paper. It calls paramsets, TrtEffvsDurPlot and adjust_delta.

####hcvhivThIcVc, hcvhivThIc, hcvhivThIcVcDecay
These functions contain the model ordinary differential equations (ODEs) for the standard, the quasi-steady state approximation, and standard with decay of treatment efficacy versions of the model, respectively. They take as input a time point, a population vector, and a parameter set, and output the evalution of the ODEs.

####paramsets
This script populates all the parameter sets used in creating figures 2-6.

####SimulateTreatmentDurations,SimulateTreatmentDurationsDecay
These functions simulate viral dynamics for a given treatment efficacy and a short and long treatment duration under 3 different HIV scenarios (HIV negative, reduced CD4 and severely immunocompromised). They take as input symbolic equilibrium formulae, a parameter set, durations, CD4 counts and HIV viral load. They call adjust_delta and VLplot or VLplotDecay, respectively,  and output time and population vectors for each scenario. SimulateTreatmentDurationsDecay performs these simulations using the assumption of treatment efficacy decay over time.

####TrtEffvsDurPlot
This function calculates and plots minimum treatment efficacy and duration pairs necessary for sustained virologic response under a range of HIV viral loads/CD4 counts. It takes as input a parameter set and an axis handle, and plots the graph on the given axis handle. It calls clear_test and adjust_delta.

####VLplot, VLplotIC, VLplotDecay
These functions take as input symbolic equilibrium formulae and parameter sets, and generate initial population vectors. Using these vectors, they run the model and get numerical solutions for the ODEs. They output time and population vectors. VLplot using the standard version of the model, VLplotIC uses the QSS approximation, and VLPlotDecay uses the standard version with treatment decay. They call adjust_delta and, respectively, the functions hcvhivThIcVc, hcvhivThIc, hcvhivThIcVcDecay, numerically evaluating them with the matlab ODE solver function ode45.












