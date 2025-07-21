% function dYdt = MFR_Reduction(t,Y,T,p,F_N2_0,pO2_in,n_MO,pO2_fun,delta_fun,delta0) %#ok<*INUSD>
% % This function calculates the reduction extent of an MFR during reduction
% % Input:
% % t -       time [s]
% % Y -       nonstoichiometry extent delta
% % T -       reduction temperature [K]
% % p -       reduction pressure [Pa]
% % F_N2_0 -  N2 sweep gas flow rate [mol/s]
% % pO2_in -  O2 partial pressure in the inlet (in sweep gas) [Pa]
% % n_MO -    metal oxide amount [mol]
% % pO2_fun - O2 partial pressure function handle (T, delta) for the chosen metal oxide
% % delta_fun -   delta as a function of phi
% % delta0 -      Initial nonstoichiometry
% % pO2 = pO2_fun(T,delta_fun(Y));
% pO2 = pO2_fun(T,delta_fun(Y));
% rO2 = (pO2-pO2_in)*F_N2_0/(p-pO2);
% % rO2 = (pO2-pO2_in)*F_N2_0/((p-pO2)*sign(p-pO2));
% % rO2 = (pO2-pO2_in)*F_N2_0/((p-pO2)*Heaviside_step_function(p-pO2,1000));
% dYdt = -(2/n_MO)*rO2;   % NEEDS JUSTIFICATION FOR THE MINUS SIGN
% % dYdt = -(2/n_MO)*rO2;
% % dYdt = (2/n_MO)*rO2*Heaviside_step_function(-rO2,1000);
% end

% ------ Version for including the iron aluminates-type of materials (delta
% decreasing with reduction) - SELECT ONLY ONE OPTION
function dYdt = MFR_Reduction(t,Y,T,p,F_N2_0,pO2_in,n_MO,pO2_fun,red_mode,delta0) %#ok<*INUSD>
% This function calculates the reduction extent of an MFR during reduction
% Input:
% t -       time [s]
% Y -       nonstoichiometry extent delta
% T -       reduction temperature [K]
% p -       reduction pressure [Pa]
% F_N2_0 -  N2 sweep gas flow rate [mol/s]
% pO2_in -  O2 partial pressure in the inlet (in sweep gas) [Pa]
% n_MO -    metal oxide amount [mol]
% pO2_fun - O2 partial pressure function handle (T, delta) for the chosen metal oxide
% red_mode -    Type of material (0-nonstoichiometry increases during reduction, 1-nonstoichiomtery decreases during reduction)
% delta0 -      Initial nonstoichiometry
pO2 = pO2_fun(T,Y);
rO2 = (pO2-pO2_in)*F_N2_0/(p-pO2);
dYdt = (2/n_MO)*rO2*(1-red_mode)+red_mode*(-(2/n_MO)*rO2);
end