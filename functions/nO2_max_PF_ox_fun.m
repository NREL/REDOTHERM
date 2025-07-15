function nO2_max = nO2_max_PF_ox_fun(T,p,x_r_in,x_p_in,omega,nO2_total,pO2_fun,K,delta_fun,phi_red)
% This function calculates the maximum oxygen exchange for parallel flow
% (PF) oxidation
% Input:
% T -           Temperature [K]
% p -           Pressure [Pa]
% x_r_in -      Equilibrium mole fraction of the reactant in oxidizer inlet stream
% x_p_in -      Equilibrium mole fraction of the reactant in oxidizer inlet stream
% omega -       Ratio of oxidizer molar flow rate to redox material flow rate [mol-ox/s to mol-redox/s]
% nO2_total -   Maximum specific O2 release per mole of redox material [mol-O2/mol-redox]
% pO2_fun -     O2 partial pressure function handle (T, delta) for the chosen metal oxide
% K -           Equilibium constant of CO2/H2O splitting reaction
% delta_fun -   delta(phi) function handle
% phi_red -     Initial O content (end of reduction)
% Output
% nO_2max -     Maximum O2 exchange for parallel flow (PF) oxidation [mol]
nO2 = linspace(0,nO2_total,1000);  % Initialize array of nO2
for I=1:length(nO2)
    n = nO2(I);
    pO2_gas = p*(K*(omega*x_r_in-2*n)/(omega*x_p_in+2*n))^2;    % Equilibrium O2 pressure in gas stream at current conditions [Pa]
    pO2_MO = pO2_fun(T,delta_fun(phi_red+2*n));                 % Equilibrium O2 pressure over metal oxide at current conditions [Pa]
    if (n==nO2_total)               % Check if we reached the maximum O2 exchange
        nO2_max = n;
        return;
    elseif (pO2_gas<pO2_MO)        % Check if O2 pressure in the ceria is larger than the gas
        nO2_max = nO2(max(I-1,1));  % Make sure that if I=1 we choose I=1 and not zero
        return;
    end
end
end