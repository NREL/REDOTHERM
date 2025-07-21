function s = Reduction_Entropy_Der_Fe33Al67(delta)
% Derivative of reduction entropy as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: s in [J/mol-K]
s = polyval(polyder([6045.9950 -3307.9447 927.7121  61.1391]),delta);
end