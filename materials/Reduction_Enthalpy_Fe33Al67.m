function h = Reduction_Enthalpy_Fe33Al67(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: h in [J/mol]
% Original fit by Janna
% h = polyval([388.6984	-243.1832	50.8679	334.5232],delta).*1e3;
% Method 1b
h = polyval([-184.104	63.5255	-95.886	-111.2472],delta).*1e3;
end