function s = Reduction_Entropy_LSM40(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent from
% Bulfin et al. (2016)
% Input: delta in [-]
% Output: s in [J/mol-K]
s = polyval([28793.4727581100	-18846.7503038220	4956.39158700211	-661.580690655142	147.805209002488],delta);
end