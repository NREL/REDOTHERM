function p_der = pO2_der_calc_fun_generic(T,delta,dH_fun,ds_fun,dH_ddelta_fun,ds_ddelta_fun,dphi_fun)
% This function calculates the derivative of the O2 partial pressure of a 
% nonstoichiometric redox material based on enthalpy and entropy of
% reduction (generic).
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
dH_ddelta = dH_ddelta_fun(delta);
ds_ddelta = ds_ddelta_fun(delta);
p_der = 2.*1e5.*(exp(-(dH-T.*ds_th)./(R.*T.*dphi_fun(delta)))).^1.*(-(dH_ddelta-T.*ds_ddelta))/(R.*T.*dphi_fun(delta));
end