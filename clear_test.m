function cleared = clear_test(eff,dur, Tcstar, Icstar,  Vcstar,Thstar,eq_index,params)

% Function to run the model with certain treatment and efficacy combos to
% see if SVR is achieved

%outputs 1 for SVR, -1 for nonresponse

% Assign treatment duration and efficacy
params.tdur = dur;
params.epsilon = eff;

% Run model
[t_out,P_out] = VLplot(Tcstar, Icstar,  Vcstar,Thstar,eq_index,params);

% check final viral load
if P_out(length(t_out),3) < 1 && ~isequal(P_out(length(t_out),3),-99)
    cleared =1;
else
    cleared = -1;
end

    