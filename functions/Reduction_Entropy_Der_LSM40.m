function ds = Reduction_Entropy_Der_LSM40(delta)
% Derivative of entropy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: ds in [J/mol-K]
ds = polyval(polyder([28793.4727581100	-18846.7503038220	4956.39158700211	-661.580690655142	147.805209002488]),delta);
end