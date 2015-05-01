% Script to create Bifurcation figure
clear
clear syms
clc
%% Set Default Figure values
clear
clc

set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize', 11)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 11)
%% Set up Equations with Quasi-steady state approximation
% Initiate symbolics variables

syms dTc dIc dTh  Tcstar Icstar Thstar 

% Write equations as strings (Quasi-steady state equations)

dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -(p/c)*beta_c*Tc*Ic';
dIc = '(p/c)*beta_c*Tc*Ic + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta_c*(1+alph*Th)*Ic';
dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*VL*Th';

%Solve for equilibria
[Tcstar, Thstar, Icstar] = solve(dTc, dTh, dIc,  'Tc','Th','Ic');

% Assign index for which equilibrium solution to use
equil_index = 4;

%% Load parameter sets
paramsets

% Load Bistable parameter set (Patient I)
temp_params = paramsALL(ind_bistable);
temp_params = rmfield(temp_params, 'delta_c_orig');

paramnames = fieldnames(temp_params);
for i = 1:numel(paramnames)
    eval(char(strcat(paramnames(i),'= temp_params.',paramnames{i},';')))
end


%% Set up values of delta, the bifurcation parameter

%get max and min value of delta_x by adjusting values for d and 3 (min and
%max values of delta/(1+alpha*H0)

delta_max_num = adjust_delta(3,Thstar, temp_params, equil_index);
delta_min_num = adjust_delta(temp_params.dc,Thstar, temp_params, equil_index);

% set how many values of delta to try
delta_res = 1200;

% set up HIV viral loads to use in bifurcation
trial_VLs = [0, 9e3, 2.5e4, 6e4, 2e5, 1e6];

delta_range = zeros(length(trial_VLs),delta_res);
delta_max = repmat(delta_max_num,size(trial_VLs));
delta_min = repmat(delta_min_num,size(trial_VLs));

% Get corresponding CD4 counts for trial viral loads

CD4s_bifurc = zeros(length(trial_VLs),2);
for i = 1:length(trial_VLs)
    for j = [1 2]
        VL = trial_VLs(i);
        if j == 1
            delta_c = delta_min(i);
        else
            delta_c = delta_max(i);
        end
        CD4s_bifurc(i,j) = eval(char(Thstar(equil_index)));
    end
end

% Break up delta ranges
for i = 1:length(delta_max)
    delta_range(i,:)= linspace(delta_min(i),delta_max(i),delta_res);
end

delta_range_check = (1+alph*CD4s_bifurc(1))*delta_range(1,:);

%back-calculate max delta
delta_index = min(find(delta_range_check>3))-1;
delta_range = delta_range(:,1:delta_index);
%% Run bifurcation  function

% Pre-allocate matrices for Infected cell counts and R0s
IL = zeros(length(delta_range),2, length(trial_VLs));
R0 = zeros(length(delta_range), length(trial_VLs));

% Run function for each delta
for i = 1: length(trial_VLs)
    [IL(:,:,i), R0(:,i)] = bifurc_ThIc(trial_VLs(i),delta_range(i,:), temp_params);
end

% save infected cell counts
IL_old = IL;
% convert infected cell counts into viral load
IL =real((temp_params.p/temp_params.c)*IL);

%% Plot bifurcation figure



%Adjust Equilibria vectors to remove zeros for plotting
F = IL;
F(~F) = nan;

% Calculate R0s for each delta (from original HIV-negative patient)
delta_inv = 1./delta_range;
Tc0 = (Tcmax/(2*r1))*((r1-dc)+sqrt((r1-dc)^2+4*sc*r1/Tcmax));
R0_orig =  ((p/c)*beta_c*Tc0 +r2*(1-Tc0/Tcmax))./(delta_range*(1+alph*CD4s_bifurc(1)));

%Get rounded CD4 counts for legend
CD4s = round(CD4s_bifurc);

%create figure
figure('units','normalized','outerposition',[0 0 1 1])
figure1 = gcf;
hold all
set(gca, 'ColorOrder',bone(7))
axes1 = gca;

% Plot stable equilibria
plot(R0_orig', squeeze(F(:,1,:)),'-','LineWidth',3);
plot(0,0)
hold on
% Plot unstable equilibria
plot(R0_orig', squeeze(F(:,2,:)),'--','LineWidth',3);
xlabel('Within-Host R_0','FontSize',22)
ylabel('HCV Viral Load (copies per mL)','FontSize',22)
set(axes1,'XLim',[0 1.05],'YLim',[0 6e6],'FontSize',22)
legend(['CD4 = ', num2str(CD4s(1)), sprintf(', HIV VL = 0')],...
    ['CD4 = ', num2str(CD4s(2)), sprintf(', HIV VL = %1.1e', trial_VLs(2))],...
    ['CD4 = ', num2str(CD4s(3)), sprintf(', HIV VL = %1.1e', trial_VLs(3))],...
    ['CD4 = ', num2str(CD4s(4)), sprintf(', HIV VL = %1.1e', trial_VLs(4))],...
    ['CD4 = ', num2str(CD4s(5)), sprintf(', HIV VL = %1.1e', trial_VLs(5))],...
    ['CD4 = ', num2str(CD4s(6)), sprintf(', HIV VL = %1.1e', trial_VLs(6))],'Location','East')



% Create textarrow
annotation(figure1,'textarrow',[R0_orig(6,length(R0_orig))/1.05+.05 R0_orig(6,length(R0_orig))/1.05+.075],...
   [F(min(find(IL(:,1:2,1)==0))-1,1,1)/max(max(max(F))) F(min(find(IL(:,1:2,1)==0))-1,1,1)/max(max(max(F)))-.1],'TextEdgeColor','none','FontSize',20,...
    'String',{'\delta = ' round(1000*delta_range(6,length(delta_range)))/1000 });

% Create textarrow
annotation(figure1,'textarrow',[R0_orig(1,min(find(IL(:,1:2,1)==0)))/1.05+.1 R0_orig(1,min(find(IL(:,1:2,1)==0)))/1.05+.05],...
   [F(min(find(IL(:,1:2,1)==0))-1,1,1)/max(max(max(F)))-.11 F(min(find(IL(:,1:2,1)==0))-1,1,1)/max(max(max(F)))-.12],'TextEdgeColor','none','FontSize',20,...
    'String',{'\delta = ' round(1000*delta_range(1,min(find(IL(:,1:2,1)==0))))/1000});

% Create textarrow
annotation(figure1,'textarrow',[0.7 0.66],...
    [0.37 0.34],'TextEdgeColor','none','FontSize',20,...
    'String',{'Unstable'});

% Create textarrow
annotation(figure1,'textarrow',[0.84 0.8],...
    [0.72 0.76],'TextEdgeColor','none','FontSize',20,...
    'String',{'Stable'});


