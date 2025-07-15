function dh = Reduction_Enthalpy_Der_LCMA(delta)
% Derivatire of the enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
dh = polyval(polyder([-4261.68344016110	3151.98240960251	-815.995941174202	-98.0106100884883	297.756258853530]),delta).*1e3;
end