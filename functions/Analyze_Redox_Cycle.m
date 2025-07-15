function Sol = Analyze_Redox_Cycle(X,S)
% Parameters space X
% X(1) - Reduction temperature (T_red) [K]
% X(2) - Oxidation temperature (T_ox) [K]
% X(3) - Reduction pressure (p_red) [Pa]
% X(4) - Oxidation pressure (p_ox) [Pa]
% X(5) - Sweep gas to oxide molar flow rate ratio (omega_red)
% X(6) - Oxidizer gas to oxide molar flow rate ratio (omega_ox)
% ---Function---
R = 8.3144598;              % Universal Gas constant [J/mol-K]
% Oxidation - pre-calculation of inlet stream compostion
if S.ox_comp_flag==0
    gas = GRI30;
    nsp = nSpecies(gas);
    iCO2 = speciesIndex(gas,'CO2');
    iH2O = speciesIndex(gas,'H2O');
    % Set gas mole fractions of CO2 to 1, and all others to zero for input
    % CO2 inlet equilibrium
    x = zeros(nsp,1);
    x(iCO2) = 1;
    % Calculate inlet state properties during oxidation
    set(gas,'T',X(2),'P',X(4),'X',x);
    equilibrate(gas,'TP');
    xCO_in = moleFraction(gas,'CO');        % Inlet CO mole fraction
    xCO2_in = moleFraction(gas,'CO2');      % Inlet CO2 mole fraction
    % H2O inlet equilibrium
    x = zeros(nsp,1);
    x(iH2O) = 1;
    % Calculate inlet state properties during oxidation
    set(gas,'T',X(2),'P',X(4),'X',x);
    equilibrate(gas,'TP');
    xH2_in = moleFraction(gas,'H2');        % Inlet H2 mole fraction
    xH2O_in = moleFraction(gas,'H2O');      % Inlet H2O mole fraction
else
    xH2_in = 1-S.x_ox_in;     % Inlet H2 mole fraction
    xH2O_in = S.x_ox_in;      % Inlet H2O mole fraction
    xCO_in = 1-S.x_ox_in;     % Inlet CO mole fraction
    xCO2_in = S.x_ox_in;      % Inlet CO2 mole fraction
end
% Select values depending on oxidizer
% g_i (i=H2,H2O,O2_WS,CO,CO2,O2_CDS) - Gibbs energy of species i [J/mol]
% DeltaG (WS/CDS) - Gibbs energy of reaction (WS/CDS) [J/mol]
% K -       equilibrium constant of the overall thermolysis reaction
% (WS/CDS)
% x_r_in -  inlet reactant oxidizer (H2O/CO2) mole fraction
% x_p_in -  inlet product (H2/CO) mole fraction (IN INLET STREAM!)
% ox_str -  string of the oxidizer type (H2O/CO2)
% prod_str - string of the product type (H2/CO)
% p_str -   string of the pressure of product at oxidation feedstock
% ox_title_str - string depicting water/CO2 splitting (WS/CDS)
% HHV -     higher heating value of product [J/mol]
% Delta_H_r - enthalpy of reaction at T_ox [J/mol]
switch S.K_input
    case 1
        [g_H2,g_H2O,g_O2_WS] = Gibbs_energy_WS(X(2));
        DeltaG_WS = g_H2+0.5*g_O2_WS-g_H2O;
        K = exp(-DeltaG_WS./(R.*X(2)));
        x_r_in = xH2O_in;
        x_p_in = xH2_in;
        HHV = 286e3;
        Delta_H_r = Reaction_Enthalpy_WS(X(2));
        S.oxidizer_type = "H2O";
    case 2
        [g_CO,g_CO2,g_O2_CDS] = Gibbs_energy_CDS(X(2));
        DeltaG_CDS = g_CO+0.5*g_O2_CDS-g_CO2;
        K = exp(-DeltaG_CDS./(R.*X(2)));
        x_r_in = xCO2_in;
        x_p_in = xCO_in;
        HHV = 283.4e3;
        Delta_H_r = Reaction_Enthalpy_CDS(X(2));
        S.oxidizer_type = "CO2";
end
% Start iterative solution to find reduction extent
delta_tol = Inf;
phi_red = Inf;
phi_ox = S.phi0;
while abs(delta_tol)>1e-5
    phi_red_old = phi_red;
    phi_ox_old = phi_ox;
    switch S.red_config
        case 1
            nO2_max_red = nO2_max_PF_red_fun(X(1),X(3),S.phi,X(5),S.nO2_total,S.pO2_fun,S.delta_fun,S.phi0);                    % Calculate maximum reaction extent (O2 transfer) for PF reduction
        case 2
            nO2_max_red = nO2_max_CF_red_fun(X(1),X(3),S.phi,X(5),S.nO2_total,S.pO2_fun,S.pO2_der_fun,S.delta_fun,S.phi0);      % Calculate maximum reaction extent (O2 transfer) for CF reduction
    end
    phi_red = S.phi0-nO2_max_red*2;                                                   % Calculate maximum phi reduction
    delta_red = S.delta_fun(phi_red);                                                 % Calculate maximum non-stoichiometry
    switch S.ox_config
        case 1
            nO2_max_ox = nO2_max_PF_ox_fun(X(2),X(4),x_r_in,x_p_in,X(6),nO2_max_red,S.pO2_fun,K,S.delta_fun,phi_red);  % Calculate maximum reaction extent (O2 transfer) for PF oxidation
        case 2
            nO2_max_ox = nO2_max_CF_ox_fun(X(2),X(4),x_r_in,x_p_in,X(6),nO2_max_red,S.pO2_fun,S.pO2_der_fun,K,S.delta_fun,phi_red);  % Calculate maximum reaction extent (O2 transfer) for CF oxidation
    end
    phi_ox = phi_red+2*nO2_max_ox;                                      % Calculate final phi for oxidation
    delta_ox = S.delta_fun(phi_ox);                                     % Calculate final delta for oxidation
    delta_tol = max(phi_ox-phi_ox_old,phi_red-phi_red_old);             % Find tolerance of phi values - difference between iterations
end
delta_delta = (delta_red-delta_ox)*((1-S.red_mode)-S.red_mode);     % Extent of reduction (delta_red-delta_ox) for oxidation
delta_phi = phi_ox-phi_red;                                         % Amount of O atoms exchanged per amount of MO [mol-O/mol-MO]
%% System level analysis
X_conv = min(delta_phi/X(6),1);                                     % Conversion extent
Q_fuel = HHV*delta_phi;                                             % Fuel specific energy [J/mol-MO]
% Oxidizier-product separation and oxidizer sensible heating
switch S.oxidizer_type
    case "H2O"
        h_H2O_ox = CP_PropsSI('HMOLAR','T',X(2),'P',X(4),'H2O');    % H2O specific enthalpy at T_ox [J/mol_H2O]
        h_H2O_satvap = CP_PropsSI('HMOLAR','Q',1,'P',X(4),'H2O');   % H2O specific enthalpy at T_sat_H2O(p_ox) [J/mol_H2O] - saturated vapor
        h_H2_ox = CP_PropsSI('HMOLAR','T',X(2),'P',X(4),'H2');      % H2 specific enthalpy at T_ox [J/mol_H2]
        h_ox = h_H2O_ox*(1-X_conv)+h_H2_ox*X_conv;                  % Effluent specific enthalpy at T_ox [J/mol]
        switch S.prod_sep_flag
            case 1      % Condensing and reboiling
                Q_ox_sep = (S.h_fg+S.Q_liq_heat)*X(6);  % Specific H2O boiling heat and liquid heating to boiling point [J/mol-MO]
                W_ox_sep = 0;                           % Work for effluent stream product separaion [J/mol-MO]
            case 2      % Mechanical vapor recompression cycle separation
                T_water_in = 25+273.15;     % Fresh water inlet temperature [K]
                eta_comp_MVR = 0.87;        % Vapor compressor efficiency
                p_MVR = X(4)*2;           % Compression pressure in the MVR cycle [Pa]
                T_MVR_in = CP_PropsSI('T','P',X(4),'Q',1,'H2O')+15;      % Inlet temperature into the MVR cycle (outlet temperature from the gas-gas HX hot side for the oxidation loop) [K]
                [~,~,~,~,~,~,~,~,~,~,~,~,~,Q_dot_extra,W_dot_comp,~,~,~] = MVR_H2_H2O(X(6),T_MVR_in,X(4),p_MVR,T_water_in,eta_comp_MVR,X_conv);
                W_ox_sep = W_dot_comp;                  % MVR compressor work in [J/mol-MO]
                Q_ox_sep = max(Q_dot_extra,0);          % Specific H2O condensing heat [J/mol-MO]
        end
        Q_ox_h = X(6)*(h_H2O_ox-h_H2O_satvap-(h_ox-h_H2O_satvap)*S.eps_g)+Q_ox_sep;     % Oxidizer heating in [J/mol-MO]
        % Pumping work
        h_ox_in = CP_PropsSI('HMOLAR','T',S.T0,'P',1e5,'H2O');              % Specific enthalpy at p0 [J/mol-H2O]
        s_ox_in = CP_PropsSI('SMOLAR','T',S.T0,'P',1e5,'H2O');              % Specific entropy at p0 [J/mol-K-CO2]
        h_ox_out_s = CP_PropsSI('HMOLAR','SMOLAR',s_ox_in,'P',X(4),'H2O');  % Specific enthalpy at p_ox [J/mol-CO2] - isentropic
    case "CO2"
        h_CO2_ox = CP_PropsSI('HMOLAR','T',X(2),'P',X(4),'CO2');    % CO2 specific enthalpy at T_ox [J/mol_CO2]
        h_CO2_0 = CP_PropsSI('HMOLAR','T',S.T0,'P',X(4),'CO2');        % CO2 specific enthalpy at T0 [J/mol_CO2]
        h_CO_ox = CP_PropsSI('HMOLAR','T',X(2),'P',X(4),'CO');      % CO specific enthalpy at T_ox [J/mol_CO]
        h_ox = h_CO2_ox*(1-X_conv)+h_CO_ox*X_conv;                  % Effluent specific enthalpy at T_ox [J/mol]
        W_ox_sep = -X(6)*R*S.T0*(X_conv*log(X_conv)+(1-X_conv)*log(1-X_conv))/S.eta_CO2_sep;    % Specific CO2 separation power [J/mol-MO]
        Q_ox_h = X(6)*(h_CO2_ox-h_CO2_0-(h_ox-h_CO2_0)*S.eps_g);                                % Oxidizer heating in [J/mol-MO]
        % Pumping work
        h_ox_in = CP_PropsSI('HMOLAR','T',S.T0,'P',1e5,'CO2');              % Specific enthalpy at p0 [J/mol-CO2]
        s_ox_in = CP_PropsSI('SMOLAR','T',S.T0,'P',1e5,'CO2');              % Specific entropy at p0 [J/mol-CO2]
        h_ox_out_s = CP_PropsSI('HMOLAR','SMOLAR',s_ox_in,'P',X(4),'CO2');  % Specific enthalpy at p_ox [J/mol-CO2] - isentropic
end
F_ox_sep = W_ox_sep/Q_fuel;             % Effluent stream product separation energy fraction
F_ox_h = Q_ox_h/Q_fuel;                 % Oxidizer heating energy fraction
% Sweep gas purification
% Starting framework for PSA separation
% xO2_in = (nO2_max_red./(nO2_max_red+X(5)))+S.phi;                   % Inlet O2 mole fraction into PSA
% xO2_out = S.phi;                                                    % Outlet O2 mole fraction
W_inert_sep = S.E_inert*X(5);                   % Inert gas separation energy in [J/mol-MO]
F_inert_sep = W_inert_sep/Q_fuel;               % Inert gas separation energy fraction
% Sensible MO heating
Q_sens_MO = integral(S.cp_s_fun,X(2),X(1),'ArrayValued',true)*S.M_MO*(1-S.eps_HR)*(1+S.f_th_loss);  % Sensible heat in [J/mol-MO]
F_sens_MO = Q_sens_MO/Q_fuel;                                                                       % Sensible heat fraction
% Chemical reduction energy
Q_red = integral(S.dH_fun,delta_ox,delta_red)*(1+S.f_th_loss);  % Reduction energy in [J/mol-MO] - DO WE NEED TO DIVIDE BY delta_delta and multiply by delta_phi?
F_red = Q_red/Q_fuel;                                           % Reduction energy fraction
% Exothermic oxidation heat
Q_ox_exo = (-Q_red+delta_phi*Delta_H_r)*S.eps_HR_ox;    % Exothermic heat recovered in [J/mol-MO]
F_ox_exo = Q_ox_exo/Q_fuel;                     % Exothermic heat fraction
% Sensible heating of inert gas
Q_inert_h = X(5)*(CP_PropsSI('HMOLAR','T',X(1),'P',X(3),'N2')-CP_PropsSI('HMOLAR','T',S.T0,'P',X(3),'N2'))*(1-S.eps_g); % Inert sweep gas heating in [J/mol-MO]
F_inert_h = Q_inert_h/Q_fuel;                                                                                           % Inert sweep gas heating fraction
% Oxidizer pumping work
if X(4)>1e5
    W_ox_pump = X(6).*(h_ox_out_s-h_ox_in)./S.eta_pump;
else
    W_ox_pump = 0;
end
F_ox_pump = W_ox_pump/Q_fuel;   % Oxidizier pumping work fraction
% Vacuum pumping work
if X(3)<1e5
    eta_p_vac = Pump_Efficiency(X(3),1e5);
    W_vac = (R.*S.T_pump/eta_p_vac).*(X(5)+nO2_max_red).*log(1e5./X(3));
else
    W_vac = 0;
end
F_vac = W_vac/Q_fuel;
% Recovered oxidation heat work
W_work_rec = min(Q_ox_h+Q_ox_exo,0)*S.eta_ox_htw;       % Recovered oxidation heat converted to work in [J/mol-MO]
F_work_rec = W_work_rec/Q_fuel;                         % Recovered oxidation heat converted to work energy fraction
% Total efficiency
F_total = (F_ox_sep+F_inert_sep+F_ox_pump+F_vac+F_work_rec)+F_sens_MO+F_red+F_inert_h+max(F_ox_h+F_ox_exo,0);  % Total energy fractions
eta = 1/F_total;            % Cycle efficiency
% Export results
Sol.HHV = HHV;
Sol.eta = eta;
Sol.delta_red = delta_red;
Sol.delta_ox = delta_ox;
Sol.delta_delta = delta_delta;
Sol.phi_red = phi_red;
Sol.phi_ox = phi_ox;
Sol.delta_phi = delta_phi;
Sol.X_conv = X_conv;    
Sol.F_ox_h = F_ox_h;
Sol.F_ox_sep = F_ox_sep;
Sol.F_inert_sep = F_inert_sep;
Sol.F_sens_MO = F_sens_MO;
Sol.F_red = F_red;
Sol.F_ox_exo = F_ox_exo;
Sol.F_inert_h = F_inert_h;
Sol.F_ox_pump = F_ox_pump;
Sol.F_vac = F_vac;
Sol.F_work_rec = F_work_rec;
Sol.Q_fuel = Q_fuel;
Sol.Q_ox_h = Q_ox_h;
Sol.W_ox_sep = W_ox_sep;
Sol.W_inert_sep = W_inert_sep;
Sol.W_ox_pump = W_ox_pump;
Sol.W_vac = W_vac;
Sol.W_work_rec = W_work_rec;
Sol.Q_sens_MO = Q_sens_MO;
Sol.Q_red = Q_red;
Sol.Q_ox_exo = Q_ox_exo;
Sol.Q_inert_h = Q_inert_h;
%% Functions
function result = CP_PropsSI(varargin)
    % Shorthand version of CoolProp for MATLAB
    result = py.CoolProp.CoolProp.PropsSI(varargin{:});
end
end