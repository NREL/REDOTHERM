function nO2_max = nO2_max_CF_ox_fun(T,p,x_r_in,x_p_in,omega,nO2_total,pO2_fun,pO2_der_fun,K,delta_fun,phi_red)
% This function calculates the maximum oxygen exchange for countercurrent
% flow (CF) oxidation
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
% Output:
% nO2_max -      Maximum O2 exchange for countercurrent flow (CF) oxidation [mol]
p_ref = 1e5;                            % Reference pressure [Pa]
nO2 = linspace(0,nO2_total,1000);       % Initialize array of nO2
for I=1:length(nO2)
    n = nO2(I);
    % Create empty arrays for pO2 and their derivatives
    pO2_MO_ar = [];
    pO2_gas_ar = [];
    % pO2_MO_der_ar = [];
    % pO2_gas_der_ar = [];
    n_i_total = linspace(0,n,1000); % Initialize the array from 0 exhange coordinate to current exchange coordinate n
    for n_i=n_i_total
        delta = 2*n_i;      % Current reduction extent (twice nO2)
        n_tag = n-n_i;      % Inverse reaction coordinate (for the counter flow)
        % Reacting gas
        pO2_gas = p_ref*(K*(omega*x_r_in-2*n_tag)/(omega*x_p_in+2*n_tag))^2;                                % Equilibrium O2 pressure in gas stream at current conditions [Pa]
        % pO2_gas_der = -4*omega*K^2*p*(x_r_in+x_p_in)*(omega*x_r_in-2*n_tag)/(omega*x_p_in+2*n_tag)^3;   % Derivative of pO2 by reaction coordinate [Pa/mol-O2]
        % Metal oxide
        pO2_MO = pO2_fun(T,delta_fun(phi_red+delta));      % Equlibrium O2 pressure [Pa]
        % pO2_MO_der = pO2_der_fun(T,delta_fun(delta));      % Derivative of MO O2 equilibrium pressure [Pa]
        % Add pO2 values to arrays
        % Metal oxide - add at the end
        pO2_MO_ar = horzcat(pO2_MO_ar,pO2_MO);
        % pO2_MO_der_ar = horzcat(pO2_MO_der_ar,pO2_MO_der);
        % Gas - add at the front (reverse coordinate)
        pO2_gas_ar = horzcat(pO2_gas,pO2_gas_ar);
        % pO2_gas_der_ar = horzcat(pO2_gas_der,pO2_gas_der_ar);
        % Check chemical potential (pressure) conditions to know if
        % maximum exchange reached or we are at the highest possible
        % exchange
        if (n==nO2_total)
            nO2_max = n;
            return;
        elseif (pO2_gas<pO2_MO)||(any(pO2_gas_ar<pO2_MO_ar))
            nO2_max = nO2(max(I-1,1));  % Make sure that if I=1 we choose I=1 and not zero
            return;
        end
    end
end
end