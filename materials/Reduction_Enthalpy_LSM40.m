function h = Reduction_Enthalpy_LSM40(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
h = polyval([18526.9125312381	-4866.74542772799	1131.38738218696	0.572640604439397	261.475684571276],delta).*1e3;
end