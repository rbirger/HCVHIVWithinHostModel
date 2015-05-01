function dvdt = hcvhivThIc(t,v,params)

%ODEs for Quasi-Steady State version of Model

Tc = v(1); %Healthy HCV Target Cells
Ic = v(2); %Infected Target Cells
Th = v(3); %CD4 Cells


%Turn on treatment if t is less that treatment duration, off if t is
%greater
if t <params.tdur 
    eps = params.epsilon;
else
    eps = 0;
end

% Cure boundary implementation: turn of p if Infected Target Cell count is
% less than 1
if Ic <1
    params.p = 0;
end

dvdt = zeros(size(v));
dvdt(1) = params.sc +params.r1*Tc*(1-(Tc +Ic)/params.Tcmax) - params.dc*Tc -(params.p*(1-eps)/params.c)*params.beta_c*Tc*Ic;
dvdt(2) = params.beta_c*(params.p*(1-eps)/params.c)*Tc*Ic + params.r2*Ic*(1-(Tc +Ic)/params.Tcmax)-params.delta_c*(1+params.alph*Th)*Ic;
dvdt(3) = params.sh*(1+params.gamm*Ic) - params.dh*Th-params.beta_h*params.VL*Th;
