function ds = Reduction_Entropy_Der_CeZr20(delta)
% Derivative of entropy of reduction as a function of nonstoichiometry extent from
% Bulfin et al. (2016)
% Input: delta in [-]
% Output: ds in [J/mol-K]
ds = polyval(polyder([-654575246.005027	484617338.809207	-143521773.042047	21752793.3562150	-1805338.38584696	81990.9071596598	-1992.16529093549	183.801571687312]),delta);
end