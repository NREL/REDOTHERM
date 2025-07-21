function nO2_max = nO2_max_PF_red_fun(T,p,phi,omega,nO2_total,pO2_fun,delta_fun,phi0)
% This function calculates the maximum oxygen exchange for parallel
% flow (PF) reduction
% Input:
% T -           Temperature [K]
% p -           Pressure [Pa]
% phi -         O2 impurity (mole fraction) in sweep gas
% omega -       Ratio of sweep gas molar flow rate to redox material flow rate [mol-sg/s to mol-redox/s]
% nO2_total -   Maximum specific O2 release per mole of redox material [mol-O2/mol-redox]
% pO2_fun -     O2 partial pressure function handle (T, delta) for the chosen metal oxide
% phi0 -        Initial phi (moles of O missing per moles of MO) [mol-O/mol-MO]
% delta_fun -   delta(phi) function handle
% Output
% nO2max -      Maximum O2 exchange for parallel flow (PF) reduction [mol]
nO2 = linspace(1e-25,nO2_total,1000); % Initialize array for nO2
% Start running for all nO2 from 0 until end or until limit is reached
for I=1:length(nO2)
    n = nO2(I);
    pO2_MO = pO2_fun(T,delta_fun(phi0-n*2));    % Equilibrium O2 partial pressure for metal oxide at current conditions [Pa]
    pO2_sg = (phi*omega+n)/(omega+n)*p;         % Equilibrium O2 partial pressure in the sweep gas at current conditions [Pa]
    % Check if sweep gas O2 pressure is larger than ceria pO2 or we are at end of the array
    if (n==nO2_total)
        nO2_max = n;
        return;
    elseif (pO2_sg>pO2_MO)
        if I==1
            nO2_max = 0;
        else
            nO2_max = nO2(I-1);
        end
        return;
    end
end
end