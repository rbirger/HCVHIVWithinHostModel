%Figures for paper
%% Set Default Figure values
clear
clc

set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize', 11)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 11)
%% Load parameter sets
paramsets;
%% Figure 2: Recreate viral load figures from Snoeck et al 2010 
%clear
% Initialize variables.
filename = 'SnoeckData_AllOutcomes.txt';
delimiter = '\t';
startRow = 2;

% Format string for each line of text:
formatSpec = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
SnoeckDataAllOutcomes = [dataArray{1:end-1}];

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

breakthrough = SnoeckDataAllOutcomes(find(~isnan(SnoeckDataAllOutcomes(:,1))),1:2);
nullResponse = SnoeckDataAllOutcomes(find(~isnan(SnoeckDataAllOutcomes(:,3))),3:4);
svResponse = SnoeckDataAllOutcomes(find(~isnan(SnoeckDataAllOutcomes(:,5))),5:6);
partialResponse = SnoeckDataAllOutcomes(find(~isnan(SnoeckDataAllOutcomes(:,7))),7:8);
relapse = SnoeckDataAllOutcomes(find(~isnan(SnoeckDataAllOutcomes(:,9))),9:10);



%% Figure 2: Set up equations and obtain analytical solutions:

%Initiate symbolics variables
syms dTc dIc dVc dTh dVh Tcstar Vcstar Thstar Vhstar

%Write equations as strings

dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -beta_c*Tc*Vc';
dIc = 'beta_c*Tc*Vc + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta_c*(1+alph*Th)*Ic';
dVc = 'p*Ic -c*Vc';
dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*VL*Th';

%Solve for equilibria
[Tcstar, Icstar,  Vcstar, Thstar] = solve(dTc, dIc, dVc, dTh,  'Tc','Ic','Vc','Th');

% Assign index for which equilibrium solution to use
equil_index = 1;

%% Figure 2: Run with different parameter sets


% Null Responder
% load Null responder parameters
temp_params = paramsALL(ind_NullResp);
%Treatment duration and efficacy
temp_params.tdur = 48*7;
temp_params.epsilon = .5;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

[t_outNull, P_outNull] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params );

% Partial response
% load Partial Responder parameters
temp_params = paramsALL(ind_Part);
%Treatment duration and efficacy
temp_params.tdur = 48*7;
temp_params.epsilon = .625;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

[t_outPart, P_outPart] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params);

% Relapse
% load Relapse parameters
temp_params = paramsALL(ind_Relapse);
%Treatment duration and efficacy
temp_params.tdur = 48*7;
temp_params.epsilon = .595;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

[t_outRelapse, P_outRelapse] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params);

% Breakthrough
% load Breakthrough parameters
temp_params = paramsALL(ind_Breakthrough);
%Treatment duration and efficacy
temp_params.tdur = 48*7;
temp_params.epsilon = .5915;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

[t_outBreak, P_outBreak] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params);

% SVR (Patient II)
% load SVR parameters
temp_params = paramsALL(ind_nonbist);
%Treatment duration and efficacy
temp_params.epsilon = .7725;
temp_params.tdur = 48*7;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
[t_outSVR, P_outSVR] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params);

%Plot Figure 2: First plot Snoeck Data
h = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
set(gca,'ColorOrder', gray(6))

for i = [5 3 7 1 9]
    subplot(2,2,1)
    hold on
    set(gca,'ColorOrder', gray(6))
    plot(SnoeckDataAllOutcomes(:,i), SnoeckDataAllOutcomes(:,i+1),'-*','LineWidth',2)
    hold on
end

%plot Detection threshold line and hour
plot(SnoeckDataAllOutcomes(:,1),2*ones(size(SnoeckDataAllOutcomes(:,1))),'k--')
legend1 =legend('SVR' ,'Null Response','Partial Response','Breakthrough','Relapse', 'Location','SouthEast');
set(legend1,...
    'Position',[0.392376077586207 0.586647727272727 0.108826809342422 0.115757335316051],...
    'FontSize',14.4);
% Create textarrow
annotation(h,'textarrow',[0.172377992897239 0.145148356054531],...
    [0.596563852813854 0.622159090909091],'String',{'Detection Threshold'},...
    'FontSize',15);
set(gca, 'FontSize', 16,'YLim',[1 9])
xlabel('Weeks')
ylabel('log_{10} HCV Viral Load')
title('a)','FontWeight','normal')

%Plot Figure 2: Model outputs
subplot(2,2,2)
hold on
set(gca,'ColorOrder', gray(6))
plot(t_outSVR/7,log10(P_outSVR(:,3)),'LineWidth',2)
plot(t_outNull/7,log10(P_outNull(:,3)),'LineWidth',2)
plot(t_outPart/7,log10(P_outPart(:,3)),'LineWidth',2)
plot(t_outBreak/7,log10(P_outBreak(:,3)),'LineWidth',2)
plot(t_outRelapse/7,log10(P_outRelapse(:,3)),'LineWidth',2)

plot(SnoeckDataAllOutcomes(:,1),2*ones(size(SnoeckDataAllOutcomes(:,1))),'k--')
annotation(h,'textarrow',[0.673616680032076 0.654370489174018],...
    [0.603693181818183 0.620738636363636],'String',{'Detection Threshold'},...
    'FontSize',16);
set(gca, 'FontSize', 16,'YLim',[1 9],'XLim',[0 80])
legend2 = legend('SVR' ,'Null Response','Partial Response','Breakthrough','Relapse', 'Location','SouthEast');
set(legend2,...
    'Position',[0.839850015036086 0.588068181818182 0.108826809342422 0.115757335316051],...
    'FontSize',14.4);

xlabel('Weeks')
ylabel('log_{10} HCV Viral Load')
title('b)','FontWeight','normal')

%%  Figure 2: Neuman data plus biphasic model scenarios

% Input Neumann Patient Data
neumann_patient_low = [0.183,5.959;0.666,5.235; 0.508,4.905;0.969,4.841;1.138, 4.544;1.960,4.827;2.920, 4.571;4.017, 4.426;5.005, 4.426;6.983, 4.213;9.006, 4.101];
neumann_patient_high = [0.225, 5.648;0.346, 4.775;0.521, 4.575;0.815, 4.234;0.957, 3.991;1.382, 3.701;1.972, 3.186;2.968, 3.305;5.005, 2.670;7.007, 2.689;8.998, 2.016];

% Load parameters for Patient A 
temp_params = paramsALL(ind_Neumann1);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment Efficacy and Duration
temp_params.tdur = 48*7;
temp_params.epsilon = .78;
[t_outSVR1, P_outSVR1] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index,temp_params);


% Load parameters for Patient B 
temp_params = paramsALL(ind_Neumann2);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment Efficacy and Duration
temp_params.tdur = 48*7;
temp_params.epsilon = .63;

[t_outSVR2, P_outSVR2] = VLplot(Tcstar, Icstar,  Vcstar, Thstar,equil_index, temp_params);


% Plot Neumann Data and model outputs
subplot(2,2,3)
hold on
set(gca,'ColorOrder', gray(5))
plot(neumann_patient_low(:,1),neumann_patient_low(:,2),'*','MarkerSize',10)
plot(neumann_patient_high(:,1),neumann_patient_high(:,2),'o','MarkerSize',10)


plot(t_outSVR1,log10(P_outSVR1(:,3)),'LineWidth',3)
plot(t_outSVR2,log10(P_outSVR2(:,3)),'LineWidth',3)
set(gca, 'FontSize', 16,'YLim',[1 9],'XLim',[0 14])
ylabel('log_{10} HCV Viral Load')
xlabel('Days')
legend('Patient A','Patient B', 'Model Output A','Model Output B','Location','NorthEast')
title('c)','FontWeight','normal')
%% Figure 2: Quasi-steady state approximation
% Set up equations and obtain analytical solutions

% Initiate symbolics variables
syms dTc dIc dTh dVh Tcstar Vcstar Thstar Vhstar

%Write equations as strings

dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -(p/c)*beta_c*Tc*Ic';
dIc = '(p/c)*beta_c*Tc*Ic + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta_c*(1+alph*Th)*Ic';
dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*VL*Th';

%Solve for equilibria
[Tcstar, Thstar, Icstar] = solve(dTc, dTh, dIc,  'Tc','Th','Ic');

% Assign index for which equilibrium solution to use
equil_index = 4; 

%% Figure 2: Run parameter sets with QSS equilibria


%Null response
temp_params = paramsALL(ind_NullResp);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment efficacy and duration
temp_params.tdur =  48*7;
temp_params.epsilon = .5;

[t_outNull, P_outNull] = VLplotIC(Tcstar,   Icstar, Thstar,equil_index, temp_params );

% Partial response
temp_params = paramsALL(ind_Part);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment efficacy and duration
temp_params.epsilon = .625;
temp_params.tdur =  48*7;

[t_outPart, P_outPart] = VLplotIC(Tcstar,   Icstar, Thstar,equil_index, temp_params);

% Relapse
temp_params = paramsALL(ind_Relapse);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment efficacy and duration
temp_params.tdur =  48*7;
temp_params.epsilon = .593;
[t_outRelapse, P_outRelapse] = VLplotIC(Tcstar,  Icstar, Thstar,equil_index, temp_params);


% Breakthrough
temp_params = paramsALL(ind_Breakthrough);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment efficacy and duration
temp_params.tdur =  48*7;
temp_params.epsilon = .5905;
[t_outBreak, P_outBreak] = VLplotIC(Tcstar, Icstar, Thstar,equil_index, temp_params);

% SVR
temp_params= paramsALL(ind_nonbist);
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);
% Treatment efficacy and duration
temp_params.tdur = 48*7;
temp_params.epsilon = .77;
[t_outSVR, P_outSVR] = VLplotIC(Tcstar,  Icstar, Thstar,equil_index, temp_params);

% Plot last panel of Figure 2 with QSS outputs
subplot(2,2,4)
hold on
set(gca,'ColorOrder', gray(6))
plot(t_outSVR/7,log10(P_outSVR(:,2)),'LineWidth',2)
plot(t_outNull/7,log10(P_outNull(:,2)),'LineWidth',2)
plot(t_outPart/7,log10(P_outPart(:,2)),'LineWidth',2)
plot(t_outBreak/7,log10(P_outBreak(:,2)),'LineWidth',2)
plot(t_outRelapse/7,log10(P_outRelapse(:,2)),'LineWidth',2)

% Plot detection threshold
plot(SnoeckDataAllOutcomes(:,1),2*ones(size(SnoeckDataAllOutcomes(:,1))),'k--')

% Create textarrow
annotation(h,'textarrow',[0.670682730923695 0.624899598393574],...
    [0.120738636363636 0.151988636363636],'String',{'Detection Threshold'},...
    'FontSize',16);
set(gca, 'FontSize', 16,'YLim',[1 9],'XLim',[0 80])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
legend4 = legend('SVR' ,'Null Response','Partial Response','Breakthrough','Relapse');
set(legend4,...
    'Position',[0.835385131835938 0.110795454545455 0.106021118164063 0.115757335316051],...
    'FontSize',14.4);

title('d)','FontWeight','normal')
%% Figure 4: show how SVR dynamics change with reduced immune system

% Initiate symbolics variables
syms dTc dVc dTh dVh Tcstar Vcstar Thstar Vhstar

% Write equations as strings

dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -beta_c*Tc*Vc';
dIc = 'beta_c*Tc*Vc + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta_c*(1+alph*Th)*Ic';
dVc = 'p*Ic -c*Vc';
dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*VL*Th';

% Assign index for which equilibrium solution to use
equil_index = 1;

% Solve for equilibria
[Tcstar, Icstar,  Vcstar, Thstar] = solve(dTc, dIc, dVc, dTh,  'Tc','Ic','Vc','Th');
%%
% Set up figure colors 
set(groot,'defaultAxesColorOrder',(1/255)*[ 230 97 1 ;253 184 99;94 60 153; 178 171 210 ],...
      'defaultAxesLineStyleOrder','-|--|:')

%%%%%%%%%%%%%
% Patient I %
%%%%%%%%%%%%%

% Load bistable parameter set (Patient I)  
temp_params = paramsALL(ind_bistable);
% Treatment efficacy
temp_params.epsilon = .8;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

% Set up durations to test
shortdur = 24*7;
longdur = 48*7;
% Set up Immune system parameters: normal CD4 recruitment rate, low CD4 recruitment rate, high HIV viral load 
norm_sh = 9;
low_sh = 5;
highVL = 1000000;
% Simulate varying treatment durations
[t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurations(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);


%Plot first panel of figure 4
scrsz = get(0,'ScreenSize');
figure('Position', scrsz);
subplot(1,2,1)
hold on
%set(gca, 'ColorOrder', gray(5))
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l1/7, log10(P_out_l1(:,3)),'LineWidth',4)
plot(t_out_l2/7, log10(P_out_l2(:,3)),'--','LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l2(1,4))),' for ',num2str(longdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(longdur/7),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 80])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
title(['a) Patient I, ', num2str(100*temp_params.epsilon),' % efficacy Treatment'],'FontWeight','Normal')

%%%%%%%%%%%%%
% Patient II %
%%%%%%%%%%%%%

% Load parameter set
temp_params = paramsALL(ind_nonbist);
% Treatment efficacy
temp_params.epsilon = .7725; 
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, 3);


% Simulate varying treatment durations
[t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurations(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);



subplot(1,2,2)
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l1/7, log10(P_out_l1(:,3)),'LineWidth',4)
plot(t_out_l2/7, log10(P_out_l2(:,3)),'--','LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l2(1,4))),' for ',num2str(longdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(longdur/7),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 80])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
title(['b) Patient II, ',num2str(100*temp_params.epsilon),' % efficacy Treatment'],'FontWeight','Normal')

%% Figure 6: Show what happens with high and declining efficacy


% Set up figure colors
set(groot,'defaultAxesColorOrder',(1/255)*[ 230 97 1 ;94 60 153; 178 171 210 ],...
      'defaultAxesLineStyleOrder','-|--|:')

%%%%%%%%%%%%%
% Patient I %
%%%%%%%%%%%%%

% Load parameter set
temp_params = paramsALL(ind_bistable);
% Treatment Efficacy
temp_params.epsilon = .95;
% Calcultate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

% Set up treatment durations
shortdur = 18.49*7;
longdur = 18.49*7;
% Set up Immune system parameters: normal CD4 recruitment rate, low CD4 recruitment rate, high HIV viral load 
norm_sh = 9;
low_sh = 5;
highVL = 1000000;
[t_out_n1,P_out_n1, t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurations(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);


% Start plotting figure 6: panel 1
scrsz = get(0,'ScreenSize');
figure('Position', scrsz);
subplot(2,2,1)
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l1/7, log10(P_out_l1(:,3)),'LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(round(shortdur/7)),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(round(longdur/7)),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(round(longdur/7)),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 30])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
%legend(['CD4 =',num2str(round(P_outSVRneg(1,4)))],['CD4 =',num2str(round(P_outSVRpos(1,4)))],['CD4 =',num2str(round(P_outSVRposhigh(1,4)))],'Location','NorthEast')
title(['a) Patient I, ',num2str(100*temp_params.epsilon),' % efficacy Treatment'],'FontWeight','Normal')


[t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurationsDecay(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);


% Plot panel 3 of Figure 6

subplot(2,2,3)
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l2/7, log10(P_out_l1(:,3)),'LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(round(shortdur/7)),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(round(longdur/7)),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(round(longdur/7)),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 30])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
title(['c) Patient I, ',num2str(100*temp_params.epsilon),' % efficacy Treatment with Decay'],'FontWeight','Normal')

%%%%%%%%%%%%%%
% Patient II %
%%%%%%%%%%%%%%

% Load parameters
temp_params = paramsALL(ind_nonbist);
% Treatment efficacy
temp_params.epsilon = .95;
% Calculate adjusted delta
temp_params.delta_c = adjust_delta(temp_params.delta_c_orig,Thstar, temp_params, equil_index);

% Set up treatment durations
shortdur = 12*7;
longdur = 12*7;
% Set up Immune system parameters: normal CD4 recruitment rate, low CD4 recruitment rate, high HIV viral load 
norm_sh = 9;
low_sh = 5;
highVL = 1000000;

[t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurations(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);

% Plot panel 2 of Figure 6
subplot(2,2,2)
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l2/7, log10(P_out_l1(:,3)),'LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(longdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(longdur/7),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 30])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
title(['b) Patient II, ',num2str(100*temp_params.epsilon),' % efficacy Treatment'],'FontWeight','Normal')

[t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurationsDecay(Tcstar,Icstar,Vcstar,Thstar,equil_index,temp_params, shortdur, ...
    longdur, norm_sh, low_sh, highVL);

%Plot panel 4 of figure 6

subplot(2,2,4)
hold on
plot(t_out_n1/7, log10(P_out_n1(:,3)),'LineWidth',3)
plot(t_out_l2/7, log10(P_out_l1(:,3)),'LineWidth',3)
plot(t_out_h1/7, log10(P_out_h1(:,3)),'LineWidth',3)
plot(t_out_n1/7, 2*ones(size(t_out_n1)),'k--')        
legend(['CD4 = ', num2str(round(P_out_n1(1,4))),' for ',num2str(shortdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_l1(1,4))),' for ',num2str(longdur/7),' wks'],...
    ['CD4 = ', num2str(round(P_out_h1(1,4))),' for ',num2str(longdur/7),' wks'])
     
set(gca, 'FontSize', 20,'YLim',[1.8 9], 'XLim',[0 30])
ylabel('log_{10} HCV Viral Load')
xlabel('Weeks')
title(['d) Patient II, ',num2str(100*temp_params.epsilon),' % efficacy Treatment with Decay'],'FontWeight','Normal')



