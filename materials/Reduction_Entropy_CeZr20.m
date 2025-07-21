function s = Reduction_Entropy_CeZr20(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent from
% Bulfin et al. (2016)
% Input: delta in [-]
% Output: s in [J/mol-K]
s = polyval([-654575246.005027	484617338.809207	-143521773.042047	21752793.3562150	-1805338.38584696	81990.9071596598	-1992.16529093549	183.801571687312],delta);
end