function s = Reduction_Entropy_Fe33Al67(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: s in [J/mol-K]
% Original fit by Janna
% s = polyval([6045.9950 -3307.9447 927.7121  61.1391],delta);
% Method 1b
s = polyval([-2270.2 1117.503 -346.298  -19.536],delta);
end