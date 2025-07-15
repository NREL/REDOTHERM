function ds = Reduction_Entropy_Der_LCMA(delta)
% Derivative of entropy of reduction as a function of nonstoichiometry extent
% Input: delta in [-]
% Output: ds in [J/mol-K]
ds = polyval(polyder([-254063.095418054	206629.588039536	-64774.6632670582	9665.93962602571	-1025.26399686063	189.514703171167]),delta);
end