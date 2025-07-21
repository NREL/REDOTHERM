function dh = Reduction_Enthalpy_Der_CeO2(delta)
% Derivative of the enthalpy of reduction as a function of nonstoichiometry
% extent
% Input: delta in [-]
% Output: dh/ddelta in [J/mol]
dh = (-1158+2.*1790.*delta+3.*23368.*delta.^2-4.*64929.*delta.^3).*1e3;
end