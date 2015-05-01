%Figure 5 Script
clear
clc
clear syms
%% Load parameter sets
paramsets

%% Run and plot Patient 1
% create parent figure
pat12 = figure('units','normalized','outerposition',[0 0 1 1]);
%get subplot handle
s1 = subplot(2,1,1,'Parent',pat12);
% run Efficacy vs duration function for patient 1
TrtEffvsDurPlot(paramsALL(ind_bistable),s1);

%% Run and plot Patient 2
% get subplot handle
s2 = subplot(2,1,2,'Parent',pat12);
% run Efficacy vs duration function for patient 1

TrtEffvsDurPlot(paramsALL(ind_nonbist),s2);

