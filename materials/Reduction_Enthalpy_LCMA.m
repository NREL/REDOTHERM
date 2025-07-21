function h = Reduction_Enthalpy_LCMA(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
h = polyval([-4261.68344016110	3151.98240960251	-815.995941174202	-98.0106100884883	297.756258853530],delta).*1e3;
end