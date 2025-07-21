function dh = Reduction_Enthalpy_Der_LSM40(delta)
% Derivative of enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: dh in [J/mol]
dh = polyval(polyder([18526.9125312381	-4866.74542772799	1131.38738218696	0.572640604439397	261.475684571276]),delta).*1e3;
end