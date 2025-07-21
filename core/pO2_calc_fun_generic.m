function p = pO2_calc_fun_generic(T,delta,dH_fun,ds_fun,dphi_fun)
% This function calculates the O2 partial pressure of a nonstoichiometric
% redox material based on enthalpy and entropy of reduction (generic)
% The values of the enthalpy and entropy can be calculated as a function of
% the nonstoichiometry extent a priori
% Input:
% T - Temperature [K]
% delta - nonstoichiometry extent
% dH - enthalpy of reduction [J/mol]
% ds_th - entropy of reduction [J/mol-K]
% dphi_fun - d(phi)/d(delta) function handle
% Output:
% p - Oxygen partial pressire [Pa]
R = 8.3144598; % Universal Gas constant [J/mol-K]
dH = dH_fun(delta);
ds_th = ds_fun(delta);
p = 1e5.*(exp((dH-T.*ds_th)./(R.*T.*dphi_fun(delta)))).^2;
end