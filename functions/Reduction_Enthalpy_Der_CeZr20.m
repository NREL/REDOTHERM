function dh = Reduction_Enthalpy_Der_CeZr20(delta)
% Derivative of enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: dh in [J/mol]
dh = polyval(polyder([-22398297.0416170	19072499.7381167	-6314600.83891391	1026326.98304749	-87433.6492450955	4216.49909987392	251.221224172867]),delta).*1e3;
end