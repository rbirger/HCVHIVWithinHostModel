function [t_out, P_out] = VLplotDecay(Tcstar, Icstar,  Vcstar, Thstar,eq_index, params)
% Inputs: Equilibrium formulae of Target Cells, Infected Cells, Free  CD4 cells

% Output time and population vectors for standard version of model with treatment efficacy decay
% given a set of parameters

% Get parameter names, and read in values for evaluation

paramnames = fieldnames(params);
for i = 1:numel(paramnames)
    eval(char(strcat(paramnames(i),'= params.',paramnames{i},';')))
end

%Get numerical values of equilibria with input parameters

Tcstar_vec = eval(char(Tcstar(eq_index)));
Icstar_vec = eval(char(Icstar(eq_index)));
Vcstar_vec = eval(char(Vcstar(eq_index)));
Thstar_vec = eval(char(Thstar(eq_index)));

% Set up initial population vector
P0 = [ Tcstar_vec; Icstar_vec; Vcstar_vec;Thstar_vec];

% set up time span
t0 = 0;
tf = 96*7;
% If there are negative values in P0, don't run and output zeros

if ~isempty(find(P0<0))
    t_out = [0:tf] ;
    P_out = -99*ones(tf+1,4);
else
    %run odesolver
    options = odeset('NonNegative',[1:4]);
    [t_out, P_out] = ode23(@hcvhivThIcVcDecay, [t0:tf], P0,options, params);
end
        
        

