function [t_out_n1,P_out_n1,t_out_l1,...
    P_out_l1,t_out_l2, P_out_l2,t_out_h1, P_out_h1] = ...
    SimulateTreatmentDurationsDecay(Tcstar,Icstar, Vcstar,Thstar,eq_index,params, shortdur, longdur, norm_sh, low_sh, highVL)

   %Function to simulate varying treatment durations for a given initial
   %efficacy that decays over time under varying CD4 conditions
   %Inputs: Equilibria formulae for Target cells, Infected Target Cells,
   %Free virus, CD4 count, equilibrium index, parameter set, short
   %duration, long duration, CD4 recruitment rate for normal patient, CD4
   %recruitment rate for treated HIV patient, HIV viral load for untreated
   %patient
   %Outputs: time and population vectors for each of the patients: normal
   %for short duration, treated HIV positive patient for short and long
   %durations, HIV positive for long durations
  
        % Set up parameters to test short duration in normal patient
        params.sh = norm_sh;
        params.tdur = shortdur;
        params.VL = 0;
        params.delta_c = adjust_delta(params.delta_c_orig,Thstar, params, eq_index);

        [t_out_n1, P_out_n1] = VLplotDecay(Tcstar, Icstar,  Vcstar, Thstar,eq_index,params);

        % Set up parameters to test short duration in reduced CD4 patient
        % (Patient on ART with partial CD4 recovery)

        params.sh = low_sh;
        params.tdur = shortdur;
        [t_out_l1, P_out_l1] = VLplotDecay(Tcstar, Icstar,  Vcstar, Thstar,eq_index,params);
        
        % Set up parameters to test long duration in reduced CD4 patient 
        params.tdur = longdur;
      
        [t_out_l2, P_out_l2] = VLplotDecay(Tcstar, Icstar,  Vcstar, Thstar,eq_index,params);

        % Set up parameters to test long duration in high HIV viral load patient 
        
        %longer duration
        params.sh = norm_sh;
        params.VL = highVL;
        params.tdur = longdur;
        
        [t_out_h1, P_out_h1] = VLplotDecay(Tcstar, Icstar,  Vcstar, Thstar,eq_index,params);

