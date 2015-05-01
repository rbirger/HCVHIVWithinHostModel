function delta_c = adjust_delta(delta_c_orig,Thstar, params, equil_index)
% Function that takes a value of delta between d and 3 and re-calculates it
% such that delta/(1+alpha*Th) = original delta
%Inputs: original delta, equilibrium formula for Th, parameters,
%equilibrium index
%Output: Adjusted delta


params.delta_c = delta_c_orig;

%get names of parameters
paramnames = fieldnames(params);

% feed parameter values into names for evaluation
for i = 1:numel(paramnames)
    eval(char(strcat(paramnames(i),'= params.',paramnames{i},';')))
end

% Set up delta_c
% calculate appropriate equilibrium value of Th (CD4s)
Thstar_vec = eval(char(Thstar(equil_index)));

% calulate new delta
delta_c = delta_c_orig/(1+params.alph*Thstar_vec);
