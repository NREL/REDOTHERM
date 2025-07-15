function s = Reduction_Entropy_LCMA(delta)
% Enthalpy of reduction as a function of nonstoichiometry extent from
% Bulfin et al. (2016)
% Input: delta in [-]
% Output: s in [J/mol-K]
s = polyval([-254063.095418054	206629.588039536	-64774.6632670582	9665.93962602571	-1025.26399686063	189.514703171167],delta);
end