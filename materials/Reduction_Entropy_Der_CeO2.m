function ds = Reduction_Entropy_Der_CeO2(delta)
% Derivative of entropy of reduction as a function of nonstoichiometry 
% extent from Bulfin et al. (2016)
% Input: delta in [-]
% Output: s in [J/mol-K]
R = 8.31446261815324;   % Universal gas constant [J/mol-K]
delta_m = 1/2.9;
ds = -R./(delta.*delta_m-delta.^2);
end