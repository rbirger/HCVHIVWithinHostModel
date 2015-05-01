function dvdt = hcvhivThIcVcDecay(t,v, params)

%ODEs for Standard Version of Model with treatment decay

Tc = v(1); %Healthy HCV Target cells
Ic = v(2); %Infected HCV Target cells
Vc = v(3); %HCV free virus
Th = v(4); %CD4 cells

% assign decay coefficient
epsdecay = -0.0025;  

%Turn on treatment if t is less that treatment duration, off if t is
%greater
if t <params.tdur
    %calculate treatment efficacy as a function of time
    eps = params.epsilon*exp(epsdecay*t);
else
    eps = 0;
end

% Cure boundary implementation: turn of p if Infected Target Cell count is
% less than 1
if Ic<1
    params.p = 0;
end
dvdt = zeros(size(v));
dvdt(1) = params.sc +params.r1*Tc*(1-(Tc +Ic)/params.Tcmax) - params.dc*Tc -params.beta_c*Tc*Vc;
dvdt(2) = params.beta_c*Tc*Vc + params.r2*Ic*(1-(Tc +Ic)/params.Tcmax)-params.delta_c*(1+params.alph*Th)*Ic;
dvdt(3) = params.p*(1-eps)*Ic - params.c*Vc;
dvdt(4) = params.sh*(1+params.gamm*Ic) - params.dh*Th-params.beta_h*params.VL*Th;
