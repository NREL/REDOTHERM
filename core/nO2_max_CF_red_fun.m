function nO2_max = nO2_max_CF_red_fun(T,p,phi,omega,nO2_total,pO2_fun,pO2_der_fun,delta_fun,phi0)
% This function calculates the maximum oxygen exchange for countercurrent
% flow (CF) reduction
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
% nO2max -      Maximum O2 exchange for countercurrent flow (CF) reduction [mol]
% Note: must not initialize nO2 to 0 but rather use a low value to avoid
% numerical issues.
nO2 = linspace(1e-25,nO2_total,1000);       % Initialize array of nO2
for I=1:length(nO2)
    n = nO2(I);
% for n=nO2
    % Create empty arrays for pO2 and their derivatives - UNCOMMENT IF
    % ARRAY CALCULATIN IS UNCOMMENTED IN THE LOOP
    % pO2_MO_ar = [];
    % pO2_sg_ar = [];
    % pO2_MO_der_ar = [];
    % pO2_sg_der_ar = [];
    n_i_total = linspace(0,n,1000); % Initialize the array from 0 exhange coordinate to current exchange coordinate n
    for n_i=n_i_total
        delta = delta_fun(phi0-2*n_i);              % Current reduction extent (twice nO2)
        n_tag = n-n_i;                              % Inverse reaction coordinate (for the counter flow)
        % Sweep gas
        pO2_sg = (phi*omega+n_tag)/(omega+n_tag)*p; % Equilibrium O2 partial pressure in the sweep gas at current conditions [Pa]
        % pO2_sg_der = -p*omega/(omega+n_tag)^2;      % Derivative of pO2 by reaction coordinate [Pa/mol-O2]
        % Metal oxide
        pO2_MO = pO2_fun(T,delta);                  % Equlibrium O2 pressure [Pa]
        % pO2_MO_der = pO2_der_fun(T,delta);          % Derivative of MO O2 equilibrium pressure [Pa]
        % Add pO2 values to arrays - NOTE: array themselves not used,
        % hence they are commented. Keeping them in case it is coupled with
        % any analysis that would like to know the array
        % Metal oxide - add at the end
        % pO2_MO_ar = horzcat(pO2_MO_ar,pO2_MO);
        % pO2_MO_der_ar = horzcat(pO2_MO_der_ar,pO2_MO_der);
        % Sweep gas - add at the front (reverse coordinate)
        % pO2_sg_ar = horzcat(pO2_sg,pO2_sg_ar);
        % pO2_sg_der_ar = horzcat(pO2_sg_der,pO2_sg_der_ar);
        % Check chemical potential (pressure) conditions to know if
        % maximum exchange reached or we are at the highest possible
        % exchange
        if (n==nO2_total)
            nO2_max = n;
            return;
        % elseif (pO2_sg>=pO2_MO)||((pO2_sg>=pO2_MO)&&(pO2_sg_der>=pO2_MO_der)) % Including the tangent condition
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
end