% function dYdt = MFR_Oxidation(t,Y,T,p,F_ox_in,K,p_ox_in,p_prod_in,pO2_in,n_MO,pO2_fun,delta_fun,delta0)
% % This function calculates the reduction extent of an MFR during oxidation
% % Input:
% % t -       time [s]
% % Y -       nonstoichiometry extent delta
% % T -       oxidation temperature [K]
% % p -       oxidation pressure [Pa]
% % F_ox_in - oxidizer flow rate [mol/s]
% % K -       thermolysis equilibrium constant
% % p_ox_in - oxidizer partial pressure in the inlet [Pa]
% % p_prod_in -   product partial pressure in the inlet [Pa]
% % pO2_in -  O2 partial pressure in the inlet (in sweep gas) [Pa]
% % n_MO -    metal oxide amount [mol]
% % pO2_fun - O2 partial pressure function handle (T, delta) for the chosen metal oxide
% % delta_fun -   delta as a function of phi
% % delta0 -      Initial nonstoichiometry
% p_ref = 1e5;
% pO2 = pO2_fun(T,delta_fun(Y));
% dzetadt = max((F_ox_in/p)*(K*p_ox_in-sqrt(pO2/p_ref)*p_prod_in)/(sqrt(pO2/p_ref)+K),0);
% rO2 = ((pO2_in-pO2)*F_ox_in+0.5*dzetadt*(p-pO2))/(p-pO2);
% ddeltadt = (2/n_MO)*rO2;
% dYdt = ddeltadt;
% end
% ------ Version for including the iron aluminates-type of materials (delta
% decreasing with reduction) - SELECT ONLY ONE OPTION
function dYdt = MFR_Oxidation(t,Y,T,p,F_ox_in,K,p_ox_in,p_prod_in,pO2_in,n_MO,pO2_fun,red_mode,delta0)
% This function calculates the reduction extent of an MFR during oxidation
% Input:
% t -       time [s]
% Y -       nonstoichiometry extent delta
% T -       oxidation temperature [K]
% p -       oxidation pressure [Pa]
% F_ox_in - oxidizer flow rate [mol/s]
% K -       thermolysis equilibrium constant
% p_ox_in - oxidizer partial pressure in the inlet [Pa]
% p_prod_in -   product partial pressure in the inlet [Pa]
% pO2_in -  O2 partial pressure in the inlet (in sweep gas) [Pa]
% n_MO -    metal oxide amount [mol]
% pO2_fun - O2 partial pressure function handle (T, delta) for the chosen metal oxide
% red_mode -    Type of material (0-nonstoichiometry increases during reduction, 1-nonstoichiomtery decreases during reduction)
% delta0 -      Initial nonstoichiometry
p_ref = 1e5;
pO2 = pO2_fun(T,Y);
dzetadt = max((F_ox_in/p)*(K*p_ox_in-sqrt(pO2/p_ref)*p_prod_in)/(sqrt(pO2/p_ref)+K),0);
rO2 = ((pO2_in-pO2)*F_ox_in+0.5*dzetadt*(p-pO2))/(p-pO2);
ddeltadt = -(2/n_MO)*rO2*(1-red_mode)+red_mode*((2/n_MO)*rO2);
dYdt = ddeltadt;
end