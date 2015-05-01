function [t_out, P_out] = VLplotIC(Tcstar,  Icstar, Thstar,eq_index,params)
% Inputs: Equilibrium formulae  of Target Cells, Infected Cells, CD4 cells

% Output time and population vectors for QSS approximation, given a set of
% parameters

% Get parameter names, and read in values for evaluation
paramnames = fieldnames(params);
for i = 1:numel(paramnames)
    eval(char(strcat(paramnames(i),'= params.',paramnames{i},';')))
end

%Get numerical values of equilibria with input parameters

Tcstar_vec = eval(char(Tcstar(eq_index)));
Icstar_vec = eval(char(Icstar(eq_index)));
Thstar_vec = eval(char(Thstar(eq_index)));

% Set up initial population vector
P0 = [ Tcstar_vec;  Icstar_vec; Thstar_vec];

% Set up time span
t0 = 0;
tf = 96*7;
% If there are negative values in P0, don't run and output zeros
if ~isempty(find(P0<0))
    t_out = [0:tf] ;
    P_out = zeros(tf+1,3);
else
    %run odesolver
    options = odeset('NonNegative',[1:3]);
    [t_out, P_out] = ode23(@hcvhivThIc, [t0:tf], P0,options, params);
end

% Convert Infected Cell count to viral load using p/c. While treatment is
% on, conversion is (1-epsilon)p/c

if ~(length(t_out)<tf)
    P_out((1:params.tdur),2) = (params.p*(1-params.epsilon)/params.c)*P_out((1:params.tdur),2);
    P_out((params.tdur+1:tf+1),2) = (params.p/params.c)*P_out((params.tdur+1:tf+1),2);
else
     P_out = zeros(length(t_out),3);
end

