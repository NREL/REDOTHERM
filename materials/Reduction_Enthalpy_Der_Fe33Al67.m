function h = Reduction_Enthalpy_Der_Fe33Al67(delta)
% Derivative of the enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
h = polyval(polyder([388.6984	-243.1832	50.8679	334.5232]),delta).*1e3;
end