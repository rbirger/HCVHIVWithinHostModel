function [IL,R0] = bifurc_ThIc(HIV_VL, delta_range,params)
%Inputs: HIV viral load to generate bifurcation figure for corresponding CD4
%count, range of deltas, parameters

%Outputs: Equilibrium values of Infected Cell count, R0

syms dTc dIc dTh dVh Tcstar Vcstar Thstar Vhstar  %symbolics variables

%Write equations as strings

dTc = 'sc + r1*Tc*(1-(Tc +Ic)/Tcmax) - dc*Tc -(p/c)*beta_c*Tc*Ic';
dIc = '(p/c)*beta_c*Tc*Ic + r2*Ic*(1-(Tc +Ic)/Tcmax)-delta*(1+alph*Th)*Ic';
dTh = 'sh*(1+gamm*Ic) - dh*Th - beta_h*HIV_VL*Th';
%dVh = 'beta_h*(k/e)*Vh*Th-delta_h*Vh';

%Solve for equilibria
[Tcstar, Thstar, Icstar] = solve(dTc, dTh, dIc,  'Tc','Th','Ic');

% assign parameter values to names for evaluation
paramnames = fieldnames(params);
for i = 1:numel(paramnames)
    eval(char(strcat(paramnames(i),'= params.',paramnames{i},';')))
end

%Solve for initial value of Healthy Target cells
Tc0 = (Tcmax/(2*r1))*((r1-dc)+sqrt((r1-dc)^2+4*sc*r1/Tcmax));

%Set up empty vectors for equilibria values
eq_vals = zeros(length(delta_range),2); %Hepatocytes
IL = zeros(length(delta_range),2); %Free virus
R0 = zeros(length(delta_range),1); %Within-host HCV R0

%Calculate equilibrium values for each value of delta for the given HIV VL
for i = 1:length(delta_range)
    % get current delta
    delta = delta_range(i);
    % get CD4 count
    Th = eval(char(Thstar(4)));
    % calculate adjusted delta
    W = delta_range(i)*(1+alph*Th);
    
    % Create vector with all possible values of Healthy Target Cell
    Tc = [0:Tcmax];
    % Solve for all values of f1
    f1 = sc + r1*Tc.*(1-((Tc)/Tcmax)) - dc*Tc ;
    % Solve for all values f1
    f2 =(r1+(p/c)*beta_c*Tcmax)*(1/r2)*(((p/c)*beta_c-r2/Tcmax)*Tc+ r2-W).*Tc;
    
    % Get differences between f2 and f1
    f12 = f2-f1;
    
    % Get sign of differences to find crossings
    signum = sign(f12);
    % set sign of exact data zeros to positive
    signum(f12==0) = 1; 
    % get zero crossings by diff ~= 0 
    eqs = Tc(diff(signum)~=0); 
    
    % Get parabola crossings and value of infected cell count at each
    % crossing
    if ~isempty(eqs) && max(size(eqs))<=2
        eq_vals(i,:) = eqs;
        eq_vals(find(eq_vals(:,1)==eq_vals(:,2)),2) = 0;
        IL(i,:) = (Tcmax/r2)*(((p/c)*beta_c-r2/Tcmax)*eq_vals(i,:)+r2-W);
    end
    
    % calculate corresponding within-host R0
    R0(i) = ((p/c)*beta_c*Tc0 +r2*(1-Tc0/Tcmax))/(delta_range(i)*(1+alph*sh/dh));
end

IL(find(eq_vals ==0))=0;
IL(find(IL<0))=0;
R0 = real(R0);
