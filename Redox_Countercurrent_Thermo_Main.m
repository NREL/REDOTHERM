%% ---------------General info---------------
% Created by: Alon Lidor (alon.lidor@nrel.gov)
% National Renewable Energy Laboratory (NREL)
% Date: Nov. 4, 2023
% This code calculates the performance of a counter-current flow (CF) and
% parallel flow (PF) redox reactor using the method developed by Bulfin
% (2019) to avoid violations of the second law of thermodynamics based on
% chemical potential calculations. It extends the original model to include
% both reduction and oxidation, as well as enthalpy and entropy
% dependencies on nonstoichiometry (delta), and the option to analyze more
% oxide materials
%--------------------------------------------------------------------------
%% Main
clearvars;
clc;
close all;
import py.CoolProp.CoolProp.*       % Using CoolProp for thermodynamic properties
addpath(genpath('./functions'));    % Adds the functions in the subfolder 'functions'
addpath(genpath('./materials'));    % Adds the functions in the subfolder 'materials'
addpath(genpath('./core'));         % Adds the functions in the subfolder 'core'
tic                                 % Start time counter
% Plotting defaults
% Settings
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultTextFontSize',14);
fontsize_l = 13;                        % Label size
fontsize_a = 12;                        % Axes text size
% Number of data points (parametric study)
n_x = 2;
n_y = 100;
% Save options - saving or not saving plots
save_input = input("Save figures? Y/N [N]: ","s");
if isempty(save_input)
    save_input = 'N';
end
switch save_input
    case 'N'
        SAVEFLAG = 0;
    case 'Y'
        SAVEFLAG = 1;
end
% MFR flag (0-no MFR, 1-solve for MFR)
MFR_flag = 0;
% Parametric sweeps over omegas - yes/no
para_input_omega = input("Perform parametric runs on omega? Y/N [N]: ","s");
if isempty(para_input_omega)
    para_input_omega = 'N';
end
% If parametric sweeps over omega performed, get user input for omega
% min/max values
if para_input_omega == 'Y'
    omega_red_min = input("Enter minimum flow rate ratio of sweep gas to oxide [omega_red_min=0.01]: ");
    if isempty(omega_red_min)
        omega_red_min = 0.01;
    end
    omega_red_max = input("Enter maximum flow rate ratio of sweep gas to oxide [omega_red_max=100]: ");
    if isempty(omega_red_max)
        omega_red_max = 100;
    end
    omega_ox_min = input("Enter minimum flow rate ratio of oxidizer to oxide [omega_ox_min=0.01]: ");
    if isempty(omega_ox_min)
        omega_ox_min = 0.01;
    end
    omega_ox_max = input("Enter maximum flow rate ratio of oxidizer to oxide [omega_ox_max=100]: ");
    if isempty(omega_ox_max)
        omega_ox_max = 100;
    end
end
% Parametric sweeps over temperatures - yes/no
para_input_T = input("Perform parametric runs on temperatures? Y/N [N]: ","s");
if isempty(para_input_omega)
    para_input_omega = 'N';
end
% If parametric sweeps over temperature performed, get user input for
% temperature min/max values
if para_input_T == 'Y'
    T_red_min = input("Enter minimum reduction temperature in deg. C [Tred_min=1400]: ");
    if isempty(T_red_min)
        T_red_min = 1673.15;
    else
        T_red_min = T_red_min+273.15;
    end
    T_red_max = input("Enter maximum reduction temperature in deg. C [Tred_max=1700]: ");
    if isempty(T_red_max)
        T_red_max = 1973.15;
    else
        T_red_max = T_red_max+273.15;
    end
    T_ox_min = input("Enter minimum oxidation temperature in deg. C [Tox_min=600]: ");
    if isempty(T_ox_min)
        T_ox_min = 873.15;
    else
        T_ox_min = T_ox_min+273.15;
    end
    T_ox_max = input("Enter maximum oxidation temperature in deg. C [Tox_max=1200]: ");
    if isempty(T_ox_max)
        T_ox_max = 1473.15;
    else
        T_ox_max = T_ox_max+273.15;
    end
end
% Constants
R = 8.3144598;              % Universal Gas constant [J/mol-K]
M_CeO2 = 172.1e-3;          % CeO2 molar mass [kg/mol]
M_CeZr = 162.33759e-3;      % Molar mass of Ce0.8Zr0.2O2 [kg/mol]
M_La = 138.9055e-3;         % Molar mass of La [kg/mol]
M_Ca = 40.08e-3;            % Molar mass of Ca [kg/mol]
M_Co = 58.933195e-3;        % Molar mass of Co [kg/mol]
M_Al = 26.98e-3;            % Molar mass of Al [kg/mol]
M_Fe = 55.845e-3;           % Molar mass of Fe [kg/mol]
M_Mn = 54.938045e-3;        % Molar mass of Mn [kg/mol]
M_Sr = 87.62e-3;            % Molar mass of Sr [kg/mol]
M_O = 15.9994e-3;           % Molar mass of O [kg/mol]
% Input
% Reduction pressure [Pa]
p_red = input("Enter reduction pressure (p_red) in Pa [p_red=100000]: ");
if isempty(p_red)
    p_red = 1e5;
end
% Oxidation pressure [Pa]
p_ox = input("Enter oxidation pressure (p_ox) in Pa [p_ox=100000]: ");
if isempty(p_ox)
    p_ox = 1e5;
end
% Reduction temperature [K]
T_red = input("Enter reduction temperature (T_red) in deg. C [T_red=1550]: ");
if isempty(T_red)
    T_red = 1823.15;
else
    T_red = T_red+273.15;
end
% Oxidtion temperature [K]
T_ox = input("Enter oxidation temperature (T_ox) in deg. C [T_ox=900]: ");
if isempty(T_ox)
    T_ox = 1173.15;
else
    T_ox = T_ox+273.15;
end
% Ratio of sweep gas molar flow rate to redox material flow rate [mol-sg/s to mol-redox/s]
omega_red = input("Enter ratio of sweep gas molar flow rate to redox material flow rate in mol-sg/s to mol-redox/s [omega_red=1]: ");
if isempty(omega_red)
    omega_red = 1;
end
% Ratio of oxidizer gas molar flow rate to redox material flow rate [mol-sg/s to mol-redox/s]
omega_ox = input("Enter ratio of oxidizer molar flow rate to redox material flow rate in mol-ox/s to mol-redox/s [omega_ox=1]: ");
if isempty(omega_ox)
    omega_ox = 1;
end
% O2 impurity in sweep gas
phi = input("Enter mole fraction of O2 impurity in sweep gas [phi=1e-5]: ");
if isempty(phi)
    phi = 1e-5;
end
% MFR input data - NOT USED IN THE SOLUTION FOR THE FLOW SYSTEMS
if MFR_flag==1
    n_MO = 1e4;         % Amount of metal oxide (for MFR) [mol]
    t_red = 600;        % Reduction duration (for MFR) [s]
    t_ox = 600;         % Oxidation duration (for MFR) [s]
end
% Oxidizer composition selection
x_ox_in = input('Enter the oxidizer mole fraction in the oxidizer stream (leave empty to use equilibrium composition at T_ox): ');
if isempty(x_ox_in)
    ox_comp_flag = 0;
else
    ox_comp_flag = 1;
end
% Select redox material
redox_material = input("Select redox material: 1-CeO2,2-CeZr20,3-LCMA,4-LSM40,5-Fe33Al67 [1-CeO2]: ");
if isempty(redox_material)
    redox_material = 1;
end
% NOTE: the initial delta for MFR can't be too low due to numerical
% stability issues (solution converges to zero), but it must not be too
% high as well otherwise it will overpredict performance
% NOTE 2: the derviative conditions have been found to be non-applicable,
% they are kept as legacy but not used in the thermodynamic functions
switch redox_material
    case 1
        M_MO = 172.1e-3;                                            % Molar mass of ceria [kg/mol]
        % -------
        % If using MATLAB version of the CeO2 material functions uncomment
        % these:
        % dH_fun = @(delta)Reduction_Enthalpy_CeO2(delta);            % Reduction enthalpy function handle
        % dS_fun = @(delta)Reduction_Entropy_CeO2(delta);             % Reduction entropy function handle
        % -------
        % If using MEX version of the CeO2 material functions uncomment
        % these:
        dH_fun = @(delta)panlener_dhfun_mex(delta);            % Reduction enthalpy function handle
        dS_fun = @(delta)panlener_dsfun_mex(delta);             % Reduction entropy function handle
        % -------
        dH_ddelta_fun = @(delta)Reduction_Enthalpy_Der_CeO2(delta); % Reduction enthalpy derivative function handle
        dS_ddelta_fun = @(delta)Reduction_Entropy_Der_CeO2(delta);  % Reduction entropy derivative function handle
        cp_s_fun = @(T)cp_ceria_only(T);                            % MO specific heat capacity function [J/kg-K]
        delta_0 = 0.01;                                             % Initial nonstoichiometry (for MFR)
        nO2_total = 0.25/2;                                         % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox] (MATERIAL DEPENDENT)
        MO_label = 'CeO_2';             % Material label
        red_mode = 0;                   % Reduction mode (0 - material that exhibits increase in nonstoichiomtery as it is reduced)
        delta0 = 0;                     % Initial delta - should be zero for "type 0" materials
        phi_fun = @(delta)(2-delta);    % Phi function handle
        dphi_fun = @(delta)(-1);        % d(phi)/d(delta) function handle
        delta_fun = @(phi)(2-phi);      % Delta function handle
        ddelta_fun = @(phi)(-1);        % d(delta)/d(phi) function handle
        phi0 = phi_fun(delta0);         % Initial phi
    case 2
        M_MO = 162.33759e-3;                                            % Molar mass of Ce0.8Zr0.2O2 [kg/mol]
        dH_fun = @(delta)Reduction_Enthalpy_CeZr20(delta);              % Reduction enthalpy function handle
        dS_fun = @(delta)Reduction_Entropy_CeZr20(delta);               % Reduction entropy function handle
        dH_ddelta_fun = @(delta)Reduction_Enthalpy_Der_CeZr20(delta);   % Reduction enthalpy derivative function handle
        dS_ddelta_fun = @(delta)Reduction_Entropy_Der_CeZr20(delta);    % Reduction entropy derivative function handle
        cp_s_fun = @(T)cp_ceria_only(T);                                % MO specific heat capacity function [J/kg-K]
        delta_0 = 0.02;                                                 % Initial nonstoichiometry (for MFR)
        nO2_total = 0.25/2;             % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox] (MATERIAL DEPENDENT)
        MO_label = 'CeZr20';            % Material label
        red_mode = 0;                   % Reduction mode (0 - material that exhibits increase in nonstoichiomtery as it is reduced)
        delta0 = 0;                     % Initial delta - should be zero for "type 0" materials
        phi_fun = @(delta)(2-delta);    % Phi function handle
        dphi_fun = @(delta)(-1);        % d(phi)/d(delta) function handle
        delta_fun = @(phi)(2-phi);      % Delta function handle
        ddelta_fun = @(phi)(-1);        % d(delta)/d(phi) function handle
        phi0 = phi_fun(delta0);         % Initial phi
    case 3
        M_MO = 0.6*M_La+0.4*M_Ca+0.6*M_Mn+0.4*M_Al+M_O*3;           % Molar mass of LCMA6464 [kg/mol] - La0.6Ca0.4Mn0.6Al0.4O3
        dH_fun = @(delta)Reduction_Enthalpy_LCMA(delta);            % Reduction enthalpy function handle
        dS_fun = @(delta)Reduction_Entropy_LCMA(delta);             % Reduction entropy function handle
        dH_ddelta_fun = @(delta)Reduction_Enthalpy_Der_LCMA(delta); % Reduction enthalpy derivative function handle
        dS_ddelta_fun = @(delta)Reduction_Entropy_Der_LCMA(delta);  % Reduction entropy derivative function handle
        cp_s_fun = @(T)140/M_MO;                                    % MO specific heat capacity function [J/kg-K]
        delta_0 = 0.03;                                             % Initial nonstoichiometry (for MFR)
        nO2_total = 0.5/2;              % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox] (MATERIAL DEPENDENT)
        MO_label = 'LCMA';              % Material label
        red_mode = 0;                   % Reduction mode (0 - material that exhibits increase in nonstoichiomtery as it is reduced)
        delta0 = 0;                     % Initial delta - should be zero for "type 0" materials
        phi_fun = @(delta)(3-delta);    % Phi function handle
        dphi_fun = @(delta)(-1);        % d(phi)/d(delta) function handle
        delta_fun = @(phi)(3-phi);      % Delta function handle
        ddelta_fun = @(phi)(-1);        % d(delta)/d(phi) function handle
        phi0 = phi_fun(delta0);         % Initial phi
    case 4
        M_MO = 0.6*M_La+0.4*M_Sr+M_O*3;                             % Molar mass of LSM40 [kg/mol] - La0.6Sr0.4O3
        dH_fun = @(delta)Reduction_Enthalpy_LSM40(delta);           % Reduction enthalpy function handle
        dS_fun = @(delta)Reduction_Entropy_LSM40(delta);            % Reduction entropy function handle
        dH_ddelta_fun = @(delta)Reduction_Enthalpy_Der_LSM40(delta);% Reduction enthalpy derivative function handle
        dS_ddelta_fun = @(delta)Reduction_Entropy_Der_LSM40(delta); % Reduction entropy derivative function handle
        cp_s_fun = @(T)140/M_MO;                                    % MO specific heat capacity function [J/kg-K]
        delta_0 = 0.03;                                             % Initial nonstoichiometry (for MFR)
        nO2_total = 0.5/2;              % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox] (MATERIAL DEPENDENT)
        MO_label = 'LSM40';             % Material label
        red_mode = 0;                   % Reduction mode (0 - material that exhibits increase in nonstoichiomtery as it is reduced)
        delta0 = 0;                     % Initial delta - should be zero for "type 0" materials
        phi_fun = @(delta)(3-delta);    % Phi function handle
        dphi_fun = @(delta)(-1);        % d(phi)/d(delta) function handle
        delta_fun = @(phi)(3-phi);      % Delta function handle
        ddelta_fun = @(phi)(-1);        % d(delta)/d(phi) function handle
        phi0 = phi_fun(delta0);         % Initial phi
    case 5
        M_MO = 0.33*M_Fe+0.67*M_Al+(4/(3-1/3))*M_O;                     % Molar mass of Fe33Al67 [kg/mol] - Fe0.33Al0.67.O4
        dH_fun = @(delta)Reduction_Enthalpy_Fe33Al67(delta);            % Reduction enthalpy function handle
        dS_fun = @(delta)Reduction_Entropy_Fe33Al67(delta);             % Reduction entropy function handle
        dH_ddelta_fun = @(delta)Reduction_Enthalpy_Der_Fe33Al67(delta); % Reduction enthalpy derivative function handle
        dS_ddelta_fun = @(delta)Reduction_Entropy_Der_Fe33Al67(delta);  % Reduction entropy derivative function handle
        cp_s_fun = @(T)168/M_MO;                                        % MO specific heat capacity function [J/kg-K]
        delta_0 = 0.26;                                                 % Initial nonstoichiometry (for MFR)
        nO2_total = 0.5/2;                  % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox] (MATERIAL DEPENDENT)
        MO_label = 'Fe33Al67';              % Material label
        red_mode = 1;                       % Reduction mode (1 - material that exhibits decrease in nonstoichiomtery as it is reduced)
        delta0 = 1/3;                       % Initial delta - the oxidized delta for this "type 1" material
        phi_fun = @(delta)(4/(3-delta));    % Phi function handle
        dphi_fun = @(delta)(4/((3-delta)^2));   % d(phi)/d(delta) function handle
        delta_fun = @(phi)((3*phi-4)/phi);  % Delta function handle
        ddelta_fun = @(phi)(4/phi^2);       % d(Delta)/d(phi) function handle
        phi0 = phi_fun(delta0);             % Initial phi
end
% Define pO2 and d(pO2)/d(delta) function handles
pO2_fun = @(T,delta)pO2_calc_fun_generic(T,delta,dH_fun,dS_fun,dphi_fun);
pO2_der_fun = @(T,delta)pO2_der_calc_fun_generic(T,delta,dH_fun,dS_fun,dH_ddelta_fun,dS_ddelta_fun,dphi_fun);
% Equilibrium constants - for oxidation analysis
[g_H2,g_H2O,g_O2_WS] = Gibbs_energy_WS(T_ox);   % Gibbs energy [J/mol];
DeltaG_WS = g_H2+0.5*g_O2_WS-g_H2O;             % Gibbs energy of WS reaction [J/mol]
K_WS = exp(-DeltaG_WS./(R.*T_ox));              % Equilibrium constant (WS)
[g_CO,g_CO2,g_O2_CDS] = Gibbs_energy_CDS(T_ox); % Gibbs energy [J/mol];
DeltaG_CDS = g_CO+0.5*g_O2_CDS-g_CO2;           % Gibbs energy of CDS reaction [J/mol]
K_CDS = exp(-DeltaG_CDS./(R.*T_ox));            % Equilibrium constant (CDS)
% Select CO2 or H2O splitting
K_input = input("Select oxidizer (1-H2O, 2-CO2) [1-H2O]: ");
if isempty(K_input)
    K_input = 1;
end
% Select if saving the data or not - leave empty to skip saving
filename_str = input('Enter data file name (leave empty to skip saving): ',"s");
%% Oxidation - calculation of inlet stream composition
% Inlet oxidier composition - either equilibirium of fixed
gas = GRI30;                    % Create GRI30 Cantera object
nsp = nSpecies(gas);            % Get number of species in object
iCO = speciesIndex(gas,'CO');   % Index for CO
iCO2 = speciesIndex(gas,'CO2'); % Index for CO2
iO2 = speciesIndex(gas,'O2');   % Index for O2
iH2 = speciesIndex(gas,'H2');   % Index for H2
iH2O = speciesIndex(gas,'H2O'); % Index for H2O
% If oxidizer gas composition is calculated from equilibrium
if ox_comp_flag==0
    % Set gas mole fractions of CO2 to 1, and all others to zero for input
    % CO2 inlet equilibrium
    x = zeros(nsp,1);
    x(iCO2) = 1;
    % Calculate inlet state properties during oxidation
    set(gas,'T',T_ox,'P',p_ox,'X',x);
    equilibrate(gas,'TP');
    xCO_in = moleFraction(gas,'CO');        % Inlet CO mole fraction
    xCO2_in = moleFraction(gas,'CO2');      % Inlet CO2 mole fraction
    xO2_CDS_in = moleFraction(gas,'O2');    % Inlet O2 mole fraction
    % H2O inlet equilibrium
    x = zeros(nsp,1);
    x(iH2O) = 1;
    % Calculate inlet state properties during oxidation
    set(gas,'T',T_ox,'P',p_ox,'X',x);
    equilibrate(gas,'TP');
    xH2_in = moleFraction(gas,'H2');        % Inlet H2 mole fraction
    xH2O_in = moleFraction(gas,'H2O');      % Inlet H2O mole fraction
    xO2_WS_in = moleFraction(gas,'O2');     % Inlet O2 mole fraction
% If oxidizer gas composition was entered
else
    xH2_in = 1-x_ox_in;     % Inlet H2 mole fraction
    xH2O_in = x_ox_in;      % Inlet H2O mole fraction
    xO2_WS_in = 0;          % Inlet O2 mole fraction
    xCO_in = 1-x_ox_in;     % Inlet CO mole fraction
    xCO2_in = x_ox_in;      % Inlet CO2 mole fraction
    xO2_CDS_in = 0;         % Inlet O2 mole fraction
end
% Calculate partial pressures - ONLY APPLICABLE TO MFR, NEEDS TO CHECK IF
% MULTIPLICATION BY OMEGA IS NEEDED OR NOT
pO2_WS_in = omega_ox*p_ox*xO2_WS_in;       % Inlet O2 partial pressure [Pa]
pH2_in = omega_ox*p_ox*xH2_in;             % Inlet H2 partial pressure [Pa]
pH2O_in = omega_ox*p_ox*xH2O_in;           % Inlet H2O partial pressure [Pa]
pO2_CDS_in = omega_ox*p_ox*xO2_CDS_in;     % Inlet O2 partial pressure [Pa]
pCO_in = omega_ox*p_ox*xCO_in;             % Inlet CO partial pressure [Pa]
pCO2_in = omega_ox*p_ox*xCO2_in;           % Inlet CO2 partial pressure [Pa]
% Select values depending on oxidizer
% K -       equilibrium constant of the overall thermolysis reaction
% (WS/CDS)
% pO2_in -  inlet O2 pressure during oxidation [Pa]
% xO2_ox_in - inlet O2 mole fraction during oxidation
% x_r_in -  inlet reactant oxidizer (H2O/CO2) mole fraction
% x_p_in -  inlet product (H2/CO) mole fraction (IN INLET STREAM!)
% p_ox_in - inlet reactant oxidizer (H2O/CO2) pressure [Pa]
% p_prod_in - inlet product (H2/CO) pressure [Pa] in oxidizer feedstock
% ox_str -  string of the oxidizer type (H2O/CO2)
% prod_str - string of the product type (H2/CO)
% p_str -   string of the pressure of product at oxidation feedstock
% ox_title_str - string depicting water/CO2 splitting (WS/CDS)
switch K_input
    case 1
        K = K_WS;
        pO2_in = pO2_WS_in;
        xO2_ox_in = xO2_WS_in;
        x_r_in = xH2O_in;
        x_p_in = xH2_in;
        p_ox_in = pH2O_in;
        p_prod_in = pH2_in;
        ox_str = "H2O";
        prod_str = 'H2';
        p_str = 'p_{H2,in}';
        ox_title_str = 'WS';
        HHV = 286e3;
    case 2
        K = K_CDS;
        pO2_in = pO2_CDS_in;
        xO2_ox_in = xO2_CDS_in;
        x_r_in = xCO2_in;
        x_p_in = xCO_in;
        p_ox_in = pCO2_in;
        p_prod_in = pCO_in;
        ox_str = "CO2";
        prod_str = 'CO';
        p_str = 'p_{CO,in}';
        ox_title_str = 'CDS';
        HHV = 283.4e3;
end
%% Overall thermo - iterative
% Initialize values to avoid code breaking and ensure going into the loop
delta_tol = Inf;        % Overall tolerance in 'phi' difference between iterations
phi_red_PF = Inf;
phi_red_CF = Inf;
phi_ox_PF = phi0;
phi_ox_CF = phi0;
if MFR_flag==1
    phi_red_MFR = Inf;
    phi_ox_MFR = phi0;
end
iter_count = 0;     % Counting iterations
while abs(delta_tol)>1e-5
    iter_count = iter_count + 1;        % Update iterations count
    % Use previous phi values as "old"
    phi_red_PF_old = phi_red_PF;
    phi_red_CF_old = phi_red_CF;
    phi_ox_PF_old = phi_ox_PF;
    phi_ox_CF_old = phi_ox_CF;
    if MFR_flag==1
        phi_red_MFR_old = phi_red_MFR;
        phi_ox_MFR_old = phi_ox_MFR;
    end
    %% Reduction
    % Solve for PF
    nO2_max_PF = nO2_max_PF_red_fun(T_red,p_red,phi,omega_red,nO2_total,pO2_fun,delta_fun,phi0);    % Calculate maximum reaction extent (O2 transfer) for PF reduction
    phi_red_PF = phi0-nO2_max_PF*2;                                                                 % Calculate maximum phi for PF reduction
    delta_red_PF = delta_fun(phi_red_PF);                                                           % Calculate reduction delta for PF
    % Solve for CF
    nO2_max_CF = nO2_max_CF_red_fun(T_red,p_red,phi,omega_red,nO2_total,pO2_fun,pO2_der_fun,delta_fun,phi0);    % Calculate maximum reaction extent (O2 transfer) for CF reduction
    phi_red_CF = phi0-nO2_max_CF*2;                                                                             % Calculate maximum phi reduction
    delta_red_CF = delta_fun(phi_red_CF);                                                                       % Calculate reduction delta for CF
    % MFR
    if MFR_flag==1
        F_sg_0 = n_MO*omega_red/t_red;                      % Sweep gas molar flow rate (for MFR) [mol/s]
        n_t_red = t_red*10;                                 % Number of time steps
        timespan_red = linspace(0,t_red,n_t_red)';          % Time vector [s]
        % MFR_red_func = @(t,Y)MFR_Reduction(t,Y,T_red,p_red,F_sg_0,p_red*phi,n_MO,pO2_fun,delta_fun,delta0); % Define ODE function handle - old version
        MFR_red_func = @(t,Y)MFR_Reduction(t,Y,T_red,p_red,F_sg_0,p_red*phi,n_MO,pO2_fun,red_mode,delta0); % Define ODE function handle
        odeoptions = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-10,'InitialStep',1e-3);               % Define ODE solver options
        % Y0 = phi_fun(delta_0);                                              % Initial conditions (initial phi)
        Y0 = delta_0;                   % Initial conditions (initial delta)
        [t_red_sol,Y] = ode45(MFR_red_func,timespan_red,Y0,odeoptions);     % Solve ODE
        % delta_red_MFR = delta_fun(Y(end));                                  % Calculate maximum reduction extent for MFR reduction
        delta_red_MFR = Y(end);                     % Calculate maximum reduction extent for MFR reduction
        phi_red_MFR = phi_fun(delta_red_MFR);       % Calculate phi for MFR reduction
    end
    %% Oxidation
    % Solve for PF
    nO2_max_PF_ox = nO2_max_PF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox,nO2_max_PF,pO2_fun,K,delta_fun,phi_red_PF);  % Calculate maximum reaction extent (O2 transfer) for PF oxidation
    phi_ox_PF = phi_red_PF+2*nO2_max_PF_ox;                                                                         % Final (minimum) phi for PF oxidation
    delta_phi_PF = phi_ox_PF-phi_red_PF;                                                                            % Change in moles of O over moles of MO [mol-O/mol-MO]
    delta_ox_PF = delta_fun(phi_ox_PF);                                                                             % Oxidation delta for PF
    delta_delta_PF = (delta_red_PF-delta_ox_PF)*((1-red_mode)-red_mode);                                            % Extent of reduction (delta_red-delta_ox) for PF oxidation
    X_PF = min(delta_phi_PF/omega_ox,1);                                                                            % Conversion extent (PF)
    prod_red_PF = 0.5*delta_phi_PF*1e6/(M_MO*1e3);                                                                  % Reduction productivity [micro-mole_O2/g_MO)
    prod_ox_PF = 2*prod_red_PF;                                                                                     % Oxidation productivity [micro-mole_prod/g_MO)
    % Solve for CF
    nO2_max_CF_ox = nO2_max_CF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox,nO2_max_CF,pO2_fun,pO2_der_fun,K,delta_fun,phi_red_CF);  % Calculate maximum reaction extent (O2 transfer) for CF oxidation
    phi_ox_CF = phi_red_CF+2*nO2_max_CF_ox;                                                                                     % Final (minimum) nonstoichiometry extent for CF oxidation
    delta_phi_CF = phi_ox_CF-phi_red_CF;                                                                                        % Change in moles of O over moles of MO [mol-O/mol-MO]
    delta_ox_CF = delta_fun(phi_ox_CF);                                                                                         % Calculate phi for CF oxidation
    delta_delta_CF = (delta_red_CF-delta_ox_CF)*((1-red_mode)-red_mode);                                                        % Extent of reduction (delta_red-delta_ox) for CF oxidation
    X_CF = min(delta_phi_CF/omega_ox,1);                                                                                        % Conversion extent (CF)
    prod_red_CF = 0.5*delta_phi_CF*1e6/(M_MO*1e3);                                                                              % Reduction productivity [micro-mole_O2/g_MO)
    prod_ox_CF = 2*prod_red_CF;                                                                                                 % Oxidation productivity [micro-mole_prod/g_MO)
    % Solve for MFR
    if MFR_flag==1
        F_ox_in = n_MO*omega_ox/t_ox;               % Oxidizer (CO2/H2O) molar flow rate [mol/s]
        n_t_ox = t_ox*10;                           % Number of time steps
        timespan_ox = linspace(0,t_ox,n_t_ox)';     % Time vector [s]
        % MFR_ox_func = @(t,Y)MFR_Oxidation(t,Y,T_ox,p_ox,F_ox_in,K,p_ox_in,p_prod_in,pO2_in,n_MO,pO2_fun,delta_fun,delta0);     % Define ODE function handle - old version
        MFR_ox_func = @(t,Y)MFR_Oxidation(t,Y,T_ox,p_ox,F_ox_in,K,p_ox_in,p_prod_in,pO2_in,n_MO,pO2_fun,red_mode,delta0);     % Define ODE function handle
        odeoptions = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-10,'InitialStep',1e-3);           % Set ODE solver options
        % Y0 = phi_fun(delta_red_MFR);           % Define initial conditions (phi)
        Y0 = delta_red_MFR;                    % Define initial conditions (delta)
        [t_ox_sol,Y] = ode45(MFR_ox_func,timespan_ox,Y0,odeoptions);        % Solve ODE
        % [t_ox_sol,Y] = ode15s(MFR_ox_func,timespan_ox,Y0,odeoptions);
        % delta_ox_MFR_min = delta_fun(Y(end));              % Final (minimum) nonstoichiometry extent for MFR oxidation
        delta_ox_MFR = Y(end);              % Final (minimum) nonstoichiometry extent for MFR oxidation
        phi_ox_MFR = phi_fun(delta_ox_MFR); % Final phi for MFR oxidation
        delta_delta_MFR = (delta_red_MFR-delta_ox_MFR)*((1-red_mode)-red_mode);                 % Extent of reduction (delta_red-delta_ox) for MFR oxidation
        delta_phi_MFR = phi_ox_MFR-phi_red_MFR;                                                 % Change in moles of O over moles of MO [mol-O/mol-MO]
        X_MFR = min(delta_phi_MFR*n_MO/(F_ox_in*t_ox),1);                                       % Conversion extent
        prod_red_MFR = 0.5*delta_phi_MFR*1e6/(M_MO*1e3);                                        % Reduction productivity [micro-mole_O2/g_MO)
        prod_ox_MFR = 2*prod_red_MFR;                                                           % Oxidation productivity [micro-mole_prod/g_MO)
    end
    % Find the largest difference between "phi" values between current and
    % previous iterations:
    if MFR_flag==1
        delta_tol = max([phi_red_PF-phi_red_PF_old,phi_red_CF-phi_red_CF_old,phi_red_MFR-phi_red_MFR_old,phi_ox_PF-phi_ox_PF_old,phi_ox_CF-phi_ox_CF_old,phi_ox_MFR-phi_ox_MFR_old]);
    else
        delta_tol = max([phi_red_PF-phi_red_PF_old,phi_red_CF-phi_red_CF_old,phi_ox_PF-phi_ox_PF_old,phi_ox_CF-phi_ox_CF_old]);
    end
end
runtime = toc;      % Stop runtime clock
%% Display results
disp('----- Reduction -----');
disp(['The maximum O2 exchange for parallel flow (PF) reduction is: ',num2str(nO2_max_PF)]);
disp(['The maximum O2 exchange for countercurrent flow (CF) reduction is: ',num2str(nO2_max_CF)]);
disp(['The maximum reduction extent (delta) for parallel flow (PF) reduction is: ',num2str(delta_red_PF)]);
disp(['The maximum reduction extent (delta) for countercurrent flow (CF) reduction is: ',num2str(delta_red_CF)]);
if MFR_flag==1
    disp(['The maximum reduction extent (delta) for mixed flow reactor (MFR) reduction is: ',num2str(delta_red_MFR)]);
end
disp('----- Oxidation -----');
disp(['The product mole fraction in the oxidizer inlet stream is: ',num2str(x_p_in)]);
disp(['The maximum O2 exchange for parallel flow (PF) oxidation is: ',num2str(nO2_max_PF_ox)]);
disp(['The maximum O2 exchange for countercurrent flow (CF) oxidation is: ',num2str(nO2_max_CF_ox)]);
disp(['The minimum reduction extent (delta) for parallel flow (PF) oxidation is: ',num2str(delta_ox_PF)]);
disp(['The minimum reduction extent (delta) for countercurrent flow (CF) oxidation is: ',num2str(delta_ox_CF)]);
disp(['The conversion extent (X) for parallel flow (PF) oxidation is: ',num2str(X_PF)]);
disp(['The conversion extent (X) for countercurrent flow (CF) oxidation is: ',num2str(X_CF)]);
if MFR_flag==1
    disp(['The minimum reduction extent (delta) for mixed flow reactor (MFR) oxidation is: ',num2str(delta_ox_MFR)]);
    disp(['The conversion extent (X) for mixed flow reactor (MFR) oxidation is: ',num2str(X_MFR)]);
end
disp('----- Cycle Performance -----');
disp(['The extent of reduction for parallel flow (PF) redox cycle is: ',num2str(delta_delta_PF)]);
disp(['The extent of reduction for countercurrent flow (CF) redox cycle is: ',num2str(delta_delta_CF)]);
if MFR_flag==1
    disp(['The extent of reduction for mixed flow reactor (MFR) redox cycle is: ',num2str(delta_delta_MFR)]);
end
%% Plotting - Preparation
% Settings default values
set(groot, 'DefaultLineLineWidth',1.5);
set(groot, 'DefaultTextFontSize',12);
% Solve for plotting
kappa = linspace(0,nO2_total,10000);                                % Reaction coordinate kappa (normalized O2 exchange)
kappa_ox_PF = linspace(0,nO2_max_PF,10000);                         % Reaction coordinate kappa (normalized O2 exchange) for PF oxidation
kappa_ox_CF = linspace(0,nO2_max_CF,10000);                         % Reaction coordinate kappa (normalized O2 exchange) for CF oxidation
pO2_MO_red = pO2_fun(T_red,delta_fun(phi0-kappa*2));                % Equilibrium O2 partial pressure for ceria - reduction [Pa]
pO2_sg_PF_red = (phi.*omega_red+kappa)./(omega_red+kappa).*p_red;   % Equilibrium O2 partial pressure in the sweep gas for PF reduction [Pa]
pO2_MO_PF_ox = pO2_fun(T_ox,delta_fun(phi_red_PF+2*kappa_ox_PF));   % Equilibrium O2 partial pressure for ceria - PF oxidation [Pa]
pO2_MO_CF_ox = pO2_fun(T_ox,delta_fun(phi_red_CF+2*kappa_ox_CF));   % Equilibrium O2 partial pressure for ceria - CF oxidation [Pa]
pO2_gas_PF_ox = p_ox.*(K.*(omega_ox.*x_r_in-2.*kappa_ox_PF)./(omega_ox.*x_p_in+2.*kappa_ox_PF)).^2; % Equilibrium O2 partial pressure in oxidizer gas stream for PF oxidation [Pa]
% Find indices of maximum reaction coordinates for all cases
[~,ind_PF_red] = min(abs(kappa-nO2_max_PF));        % kappa_max index for parallel flow reduction
[~,ind_CF_red] = min(abs(kappa-nO2_max_CF));        % kappa_max index for countercurrent flow reduction
[~,ind_PF_ox] = min(abs(kappa-nO2_max_PF_ox));      % kappa_max index for parallel flow oxidation
[~,ind_CF_ox] = min(abs(kappa-nO2_max_CF_ox));      % kappa_max index for countercurrent flow oxidation
kappa_red_max = ceil(max(nO2_max_PF,nO2_max_CF)*100)/100;       % Find maximum kappa for reduction
kappa_ox_max = ceil(max(nO2_max_PF_ox,nO2_max_CF_ox)*100)/100;  % Find maximum kappa for oxidation
% Calculate O2 partial pressures for CF sweep gas during reduction
pO2_sg_CF_red = zeros(1,ind_CF_red);
for I=1:ind_CF_red
    n_tag = nO2_max_CF-kappa(I);                                      % Inverse reaction coordinate (for the countercurrent flow reduction)
    pO2_sg_CF_red(I) = (phi*omega_red+n_tag)/(omega_red+n_tag)*p_red;   % Equilibrium O2 pressure [Pa]
end
% Calculate O2 partial pressures for CF oxidizer gas during oxidation
pO2_gas_CF_ox = zeros(1,ind_CF_ox);
for I=1:ind_CF_ox
    n_tag = nO2_max_CF_ox-kappa_ox_CF(I);         % Inverse reaction coordinate (for the countercurrent flow oxidation)
    pO2_gas_CF_ox(I) = p_ox*(K*(omega_ox*x_r_in-2*n_tag)/(omega_ox*x_p_in+2*n_tag))^2; % Gas O2 pressure [Pa]
end
%% Plotting
% Plot reduction
fig = figure(1);
ax = gca;
plot(kappa,pO2_MO_red,'-k');
hold on;
plot(kappa(1:ind_PF_red),pO2_sg_PF_red(1:ind_PF_red),'--k');
plot(kappa(1:ind_CF_red),pO2_sg_CF_red(1:ind_CF_red),':k');
xlim([0 kappa_red_max]);
ylim([phi*p_red 1e5]);
% ylim([0.01 1e4]);
yscale('log');
xlabel('\kappa');
ylabel('p_{O_2}, Pa');
line(ax,[kappa(ind_PF_red) kappa(ind_PF_red)],[phi*p_red pO2_sg_PF_red(ind_PF_red)],'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',1);
legend({MO_label,'Gas (PF)','Gas (CF)'},'Location','northeast');
ax.FontSize = fontsize_a;
if SAVEFLAG
    fig_name = 'Reduction';
    savefig(fig,fig_name);
    saveas(fig,fig_name,'epsc');
    saveas(fig,fig_name,'emf');
    print(fig,fig_name,'-r1000','-dpng');
end
% Plot oxidation
fig = figure(2);
ax = gca;
plot(kappa_ox_PF,pO2_MO_PF_ox,'-k');
hold on;
plot(kappa(1:ind_PF_ox),pO2_gas_PF_ox(1:ind_PF_ox),'--k');
plot(kappa_ox_CF,pO2_MO_CF_ox,'-.k');
plot(kappa(1:ind_CF_ox),pO2_gas_CF_ox(1:ind_CF_ox),':k');
xlim([0 kappa_ox_max]);
ylim([1e-10 1e5]);
yscale('log');
xlabel('\kappa');
ylabel('p_{O_2}, Pa');
line(ax,[kappa(ind_PF_ox) kappa(ind_PF_ox)],[1e-10 pO2_gas_PF_ox(ind_PF_ox)],'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',1);
line(ax,[kappa(ind_CF_ox) kappa(ind_CF_ox)],[1e-10 pO2_gas_CF_ox(ind_CF_ox)],'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',1);
legend({strcat(MO_label,' (PF)'),'Gas (PF)',strcat(MO_label,' (CF)'),'Gas (CF)'},'Location','best');
ax.FontSize = fontsize_a;
if SAVEFLAG
    fig_name = 'Oxidation';
    savefig(fig,fig_name);
    saveas(fig,fig_name,'epsc');
    saveas(fig,fig_name,'emf');
    print(fig,fig_name,'-r1000','-dpng');
end
%% Parametric study - Omega
if para_input_omega=='Y'
    omega_red_par = logspace(log10(omega_red_min),log10(omega_red_max),n_x);    % Create array of omega_red for parametric sweeps
    % Add default omega_red to array, sort, remove duplicates (in case
    % omega_red was part of the original omega_red_apr
    omega_red_par = [omega_red_par omega_red];
    omega_red_par = unique(sort(omega_red_par));
    omega_ox_par = logspace(log10(omega_ox_min),log10(omega_ox_max),n_y);       % Create array of omega_ox for parametric sweeps
    % Add default omega_ox to array, sort, remove duplicates (in case
    % omega_ox was part of the original omega_ox_apr
    omega_ox_par = [omega_ox_par omega_ox];
    omega_ox_par = unique(sort(omega_ox_par));
    nO2_red_PF_par_omega = zeros(length(omega_red_par),1);
    nO2_red_CF_par_omega = zeros(length(omega_red_par),1);
    phi_red_PF_par_omega = zeros(length(omega_red_par),1);
    phi_red_CF_par_omega = zeros(length(omega_red_par),1);
    delta_red_PF_par_omega = zeros(length(omega_red_par),1);
    delta_red_CF_par_omega = zeros(length(omega_red_par),1);
    nO2_ox_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    nO2_ox_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    phi_ox_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    phi_ox_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    delta_ox_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    delta_ox_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    if MFR_flag==1
        phi_red_MFR_par_omega = zeros(length(omega_red_par),1);
        delta_red_MFR_par_omega = zeros(length(omega_red_par),1);
        phi_ox_MFR_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
        delta_ox_MFR_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
        prod_ox_MFR_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
    end
    for I=1:length(omega_red_par)
        % Solve for PF
        nO2_red_PF_par_omega(I) = nO2_max_PF_red_fun(T_red,p_red,phi,omega_red_par(I),nO2_total,pO2_fun,delta_fun,phi0);    % Calculate maximum reaction extent (O2 transfer) for PF reduction
        phi_red_PF_par_omega(I) = phi0-2*nO2_red_PF_par_omega(I);           % Calculate maximum phi for PF
        delta_red_PF_par_omega(I) = delta_fun(phi_red_PF_par_omega(I));     % Calculate maximum reduction extent for PF reduction
        % Solve for CF
        nO2_red_CF_par_omega(I) = nO2_max_CF_red_fun(T_red,p_red,phi,omega_red_par(I),nO2_total,pO2_fun,pO2_der_fun,delta_fun,phi0);   % Calculate maximum reaction extent (O2 transfer) for CF reduction
        phi_red_CF_par_omega(I) = phi0-2*nO2_red_CF_par_omega(I);           % Calculate maximum phi for CF
        delta_red_CF_par_omega(I) = delta_fun(phi_red_CF_par_omega(I));     % Calculate maximum reduction extent for CF reduction
        % Solve for MFR
        if MFR_flag==1
            F_sg_0_par = n_MO*omega_red_par(I)/t_red;         % Sweep gas molar flow rate (for MFR) [mol/s]
            % MFR_red_func =
            % @(t,Y)MFR_Reduction(t,Y,T_red,p_red,F_sg_0_par,p_red*phi,n_MO,pO2_fun,delta_fun,delta0);    % Define ODE function handle - old version
            MFR_red_func = @(t,Y)MFR_Reduction(t,Y,T_red,p_red,F_sg_0_par,p_red*phi,n_MO,pO2_fun,red_mode,delta0);  % Define ODE function handle
            % Y0 = delta_0;       % Define initial conditions - TOO LOW VALUE CAUSES NUMERICAL INSTABILITY
            Y0 = nO2_red_CF_par_omega(I);       % Define initial conditions - TOO LOW VALUE CAUSES NUMERICAL INSTABILITY
            [t_red_sol,Y] = ode45(MFR_red_func,timespan_red,Y0,odeoptions);
            % [t_red_sol,Y] = ode15s(MFR_red_func,timespan_red,Y0,odeoptions);
            delta_red_MFR_par_omega(I) = Y(end);         % Calculate maximum reduction extent for MFR reduction
            phi_red_MFR_par_omega(I) = phi_fun(delta_red_MFR_par_omega(I));   % Calculate phi during reduction for MFR
        end
        % Solve for oxidation, covering all oxidation molar flow rates for
        % the specific reduction molar flow rate
        for J=1:length(omega_ox_par)
            % PF:
            nO2_ox_PF_par_omega(I,J) = nO2_max_PF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox_par(J),nO2_red_PF_par_omega(I),pO2_fun,K,delta_fun,phi_red_PF_par_omega(I));   % Calculate maximum reaction extent (O2 transfer) for PF oxidation
            phi_ox_PF_par_omega(I,J) = phi_red_PF_par_omega(I)+2*nO2_ox_PF_par_omega(I,J);      % Final phi for PF oxidation
            delta_ox_PF_par_omega(I,J) = delta_fun(phi_ox_PF_par_omega(I,J));            % Final (minimum) nonstoichiometry extent for PF oxidation
            % CF:
            nO2_ox_CF_par_omega(I,J) = nO2_max_CF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox_par(J),nO2_red_CF_par_omega(I),pO2_fun,pO2_der_fun,K,delta_fun,phi_red_CF_par_omega(I));   % Calculate maximum reaction extent (O2 transfer) for CF oxidation
            phi_ox_CF_par_omega(I,J) = phi_red_CF_par_omega(I)+2*nO2_ox_CF_par_omega(I,J);      % Final phi for CF oxidation
            delta_ox_CF_par_omega(I,J) = delta_fun(phi_ox_CF_par_omega(I,J));            % Final (minimum) nonstoichiometry extent for CF oxidation
            % Solve for MFR
            if MFR_flag==1
                Y0 = delta_red_MFR_par_omega(I);           % Define initial conditions
                F_ox_in_par = n_MO*omega_ox_par(J)/t_ox;       % Oxidizer (CO2/H2O) molar flow rate [mol/s]
                MFR_ox_func = @(t,Y)MFR_Oxidation(t,Y,T_ox,p_ox,F_ox_in_par,K,p_ox_in,p_prod_in,pO2_in,n_MO,pO2_fun,red_mode,delta0);     % Define ODE function handle
                [t_ox_sol,Y] = ode45(MFR_ox_func,timespan_ox,Y0,odeoptions);
                % [t_ox_sol,Y] = ode15s(MFR_ox_func,timespan_ox,Y0,odeoptions);
                delta_ox_MFR_par_omega(I,J) = Y(end);                                    % Final (minimum) nonstoichiometry extent for MFR oxidation
                phi_ox_MFR_par_omega(I,J) = phi_fun(delta_ox_MFR_par_omega(I,J));      % Final phi for MFR oxidation
            end
        end
    end
    delta_delta_PF_par_omega = 2*nO2_ox_PF_par_omega;                                   % Extent of reduction (delta_red-delta_ox) for PF oxidation
    delta_delta_CF_par_omega = 2*nO2_ox_CF_par_omega;                                   % Extent of reduction (delta_red-delta_ox) for CF oxidation
    delta_phi_PF_par_omega = phi_ox_PF_par_omega-phi_red_PF_par_omega;                  % Extent of change in atoms of O over moles of MO - PF
    delta_phi_CF_par_omega = phi_ox_CF_par_omega-phi_red_CF_par_omega;                  % Extent of change in atoms of O over moles of MO - CF
    prod_ox_PF_par_omega = delta_phi_PF_par_omega.*1e6./(M_MO.*1e3);                    % Oxidation productivity [micro-mole_prod/g_MO)
    prod_ox_CF_par_omega = delta_phi_CF_par_omega.*1e6./(M_MO.*1e3);                    % Oxidation productivity [micro-mole_prod/g_MO)
    X_PF_par_omega = min(delta_phi_PF_par_omega./omega_ox_par,1);                       % Conversion extent - PF
    X_CF_par_omega = min(delta_phi_CF_par_omega./omega_ox_par,1);                       % Conversion extent - CF
    if MFR_flag==1
        delta_delta_MFR_par_omega = delta_red_MFR_par_omega-delta_ox_MFR_par_omega;         % Extent of reduction (delta_red-delta_ox) for MFR oxidation
        delta_phi_MFR_par_omega = phi_ox_MFR_par_omega-phi_red_MFR_par_omega;               % Extent of change in atoms of O over moles of MO - MFR
        prod_ox_MFR_par_omega = delta_phi_MFR_par_omega.*1e6./(M_MO.*1e3);                  % Oxidation productivity [micro-mole_prod/g_MO)
        X_MFR_par_omega = min(delta_phi_MFR_par_omega./omega_ox_par,1);                     % Conversion extent - MFR
    end
    % Clean unphysical results (delta_delta<0 and X<0)
    [row,col] = find(delta_delta_PF_par_omega<0);
    delta_delta_PF_par_omega(row,col) = NaN;
    [row,col] = find(delta_delta_CF_par_omega<0);
    delta_delta_CF_par_omega(row,col) = NaN;
    [row,col] = find(X_PF_par_omega<0);
    X_PF_par_omega(row,col) = NaN;
    [row,col] = find(X_CF_par_omega<0);
    X_CF_par_omega(row,col) = NaN;
    if MFR_flag==1
        [row,col] = find(delta_delta_MFR_par_omega<0);
        delta_delta_MFR_par_omega(row,col) = NaN;
        [row,col] = find(X_MFR_par_omega<0);
        X_MFR_par_omega(row,col) = NaN;
    end
end
%% Parametric study - Temperatures
if para_input_T=='Y'
    T_red_par = linspace(T_red_min,T_red_max,n_x);
    T_red_par = [T_red_par T_red];
    T_red_par = unique(sort(T_red_par));
    T_ox_par = linspace(T_ox_min,T_ox_max,n_y);
    T_ox_par = [T_ox_par T_ox];
    T_ox_par = unique(sort(T_ox_par));
    nO2_red_PF_par_T = zeros(length(T_red_par),1);
    nO2_red_CF_par_T = zeros(length(T_red_par),1);
    delta_red_PF_par_T = zeros(length(T_red_par),1);
    delta_red_CF_par_T = zeros(length(T_red_par),1);
    phi_red_PF_par_T = zeros(length(T_red_par),1);
    phi_red_CF_par_T = zeros(length(T_red_par),1);
    nO2_ox_PF_par_T = zeros(length(T_red_par),length(T_ox_par));
    nO2_ox_CF_par_T = zeros(length(T_red_par),length(T_ox_par));
    phi_ox_PF_par_T = zeros(length(T_red_par),length(T_ox_par));
    phi_ox_CF_par_T = zeros(length(T_red_par),length(T_ox_par));
    delta_ox_PF_par_T = zeros(length(T_red_par),length(T_ox_par));
    delta_ox_CF_par_T = zeros(length(T_red_par),length(T_ox_par));
    if MFR_flag==1
        delta_red_MFR_par_T = zeros(length(T_red_par),1);
        phi_red_MFR_par_T = zeros(length(T_red_par),1);
        phi_ox_MFR_par_T = zeros(length(T_red_par),length(T_ox_par));
        delta_ox_MFR_par_T = zeros(length(T_red_par),length(T_ox_par));
    end
    K_par = zeros(1,length(T_ox_par));
    x_p_in_par = zeros(1,length(T_ox_par));
    x_r_in_par = zeros(1,length(T_ox_par));
    xO2_ox_in_par = zeros(1,length(T_ox_par));
    pO2_in_par = zeros(1,length(T_ox_par));
    p_prod_in_par = zeros(1,length(T_ox_par));
    p_ox_in_par = zeros(1,length(T_ox_par));
    for I=1:length(T_red_par)
        % Solve for PF
        nO2_red_PF_par_T(I) = nO2_max_PF_red_fun(T_red_par(I),p_red,phi,omega_red,nO2_total,pO2_fun,delta_fun,phi0);             % Calculate maximum reaction extent (O2 transfer) for PF reduction
        phi_red_PF_par_T(I) = phi0-nO2_red_PF_par_T(I)*2;           % Calculate maximum phi for PF
        delta_red_PF_par_T(I) = delta_fun(phi_red_PF_par_T(I));     % Calculate maximum reduction extent for PF reduction
        % Solve for CF
        nO2_red_CF_par_T(I) = nO2_max_CF_red_fun(T_red_par(I),p_red,phi,omega_red,nO2_total,pO2_fun,pO2_der_fun,delta_fun,phi0); % Calculate maximum reaction extent (O2 transfer) for CF reduction
        phi_red_CF_par_T(I) = phi0-nO2_red_CF_par_T(I)*2;           % Calculate maximum phi for CF
        delta_red_CF_par_T(I) = delta_fun(phi_red_CF_par_T(I));     % Calculate maximum reduction extent for CF reduction
        % Solve for MFR
        if MFR_flag==1
            MFR_red_func = @(t,Y)MFR_Reduction(t,Y,T_red_par(I),p_red,F_sg_0,p_red*phi,n_MO,pO2_fun,red_mode,delta0);  % Define ODE function handle
            % Y0 = delta_0;       % Define initial conditions - TOO LOW VALUE CAUSES NUMERICAL INSTABILITY
            Y0 = nO2_red_CF_par_T(I);       % Define initial conditions - TOO LOW VALUE CAUSES NUMERICAL INSTABILITY
            [t_red_sol,Y] = ode15s(MFR_red_func,timespan_red,Y0,odeoptions);
            delta_red_MFR_par_T(I) = Y(end);                        % Calculate maximum reduction extent for MFR reduction
            phi_red_MFR_par_T(I) = phi_fun(delta_red_MFR_par_T(I)); % Calculate phi during reduction for MFR
        end
        for J=1:length(T_ox_par)
            % Calculate new K            
            % Select values depending on oxidizer
            switch K_input
                case 1
                    [g_H2_par,g_H2O_par,g_O2_WS_par] = Gibbs_energy_WS(T_ox_par(J));    % Gibbs energy [J/mol];
                    DeltaG_CDS_par = g_H2_par+0.5*g_O2_WS_par-g_H2O_par;                % Gibbs energy of WS reaction [J/mol]
                    K_par(J) = exp(-DeltaG_CDS_par./(R.*T_ox_par(J)));                  % Equilibrium constant (WS)
                    % CO2 inlet equilibrium
                    x = zeros(nsp,1);
                    x(iH2O) = 1;
                    % Calculate inlet state properties during oxidation
                    set(gas,'T',T_ox_par(J),'P',p_ox,'X',x);
                    equilibrate(gas,'TP');
                    x_p_in_par(J) = moleFraction(gas,'H2');             % Inlet H2 mole fraction
                    x_r_in_par(J) = moleFraction(gas,'H2O');            % Inlet H2O mole fraction
                    xO2_ox_in_par(J) = moleFraction(gas,'O2');          % Inlet O2 mole fraction
                case 2
                    [g_CO_par,g_CO2_par,g_O2_CDS_par] = Gibbs_energy_CDS(T_ox_par(J));  % Gibbs energy [J/mol];
                    DeltaG_CDS_par = g_CO_par+0.5*g_O2_CDS_par-g_CO2_par;               % Gibbs energy of CDS reaction [J/mol]
                    K_par(J) = exp(-DeltaG_CDS_par./(R.*T_ox_par(J)));                  % Equilibrium constant (CDS)
                    % CO2 inlet equilibrium
                    x = zeros(nsp,1);
                    x(iCO2) = 1;
                    % Calculate inlet state properties during oxidation
                    set(gas,'T',T_ox_par(J),'P',p_ox,'X',x);
                    equilibrate(gas,'TP');
                    x_p_in_par(J) = moleFraction(gas,'CO');             % Inlet CO mole fraction
                    x_r_in_par(J) = moleFraction(gas,'CO2');            % Inlet CO2 mole fraction
                    xO2_ox_in_par(J) = moleFraction(gas,'O2');          % Inlet O2 mole fraction
            end
            if ox_comp_flag==1      % If oxidizer gas composition is fixed and not solved by equilibrium
                xO2_ox_in_par(J) = 0;
                x_p_in_par(J) = 1-x_ox_in;
                x_r_in_par(J) = x_ox_in;
            end
            pO2_in_par(J) = omega_ox*p_ox*xO2_ox_in_par(J);     % Inlet O2 partial pressure [Pa]
            p_prod_in_par(J) = omega_ox*p_ox*x_p_in_par(J);     % Inlet CO partial pressure [Pa]
            p_ox_in_par(J) = omega_ox*p_ox*x_r_in_par(J);       % Inlet CO2 partial pressure [Pa]
            % Solve for PF
            nO2_ox_PF_par_T(I,J) = nO2_max_PF_ox_fun(T_ox_par(J),p_ox,x_r_in_par(J),x_p_in_par(J),omega_ox,nO2_red_PF_par_T(I),pO2_fun,K_par(J),delta_fun,phi_red_PF_par_T(I));   % Calculate maximum reaction extent (O2 transfer) for PF oxidation
            phi_ox_PF_par_T(I,J) = phi_red_PF_par_T(I)+2*nO2_ox_PF_par_T(I,J);              % Final phi for PF oxidation
            delta_ox_PF_par_T(I,J) = delta_fun(phi_ox_PF_par_T(I,J));                       % Final (minimum) nonstoichiometry extent for PF oxidation
            % Solve for CF
            nO2_ox_CF_par_T(I,J) = nO2_max_CF_ox_fun(T_ox_par(J),p_ox,x_r_in_par(J),x_p_in_par(J),omega_ox,nO2_red_CF_par_T(I),pO2_fun,pO2_der_fun,K_par(J),delta_fun,phi_red_CF_par_T(I));   % Calculate maximum reaction extent (O2 transfer) for CF oxidation
            phi_ox_CF_par_T(I,J) = phi_red_CF_par_T(I)+2*nO2_ox_CF_par_T(I,J);              % Final phi for CF oxidation
            delta_ox_CF_par_T(I,J) = delta_fun(phi_ox_CF_par_T(I,J));                       % Final (minimum) nonstoichiometry extent for CF oxidation
            % Solve for MFR
            if MFR_flag==1
                MFR_ox_func = @(t,Y)MFR_Oxidation(t,Y,T_ox_par(J),p_ox,F_ox_in,K_par(J),p_ox_in_par(J),p_prod_in_par(J),pO2_in_par(J),n_MO,pO2_fun,red_mode,delta0);     % Define ODE function handle
                Y0 = delta_red_MFR_par_T(I);           % Define initial conditions
                [t_ox_sol,Y] = ode45(MFR_ox_func,timespan_ox,Y0,odeoptions);
                % [t_ox_sol,Y] = ode15s(MFR_ox_func,timespan_ox,Y0,odeoptions);
                delta_ox_MFR_par_T(I,J) = Y(end);              % Final (minimum) nonstoichiometry extent for MFR oxidation
                phi_ox_MFR_par_T(I,J) = phi_fun(delta_ox_MFR_par_T(I,J));   % Phi at oxidation for MFR
            end
        end
    end
    delta_delta_PF_par_T = 2*nO2_ox_PF_par_T;                           % Extent of reduction (delta_red-delta_ox) for PF oxidation
    delta_delta_CF_par_T = 2*nO2_ox_CF_par_T;                           % Extent of reduction (delta_red-delta_ox) for CF oxidation
    delta_phi_PF_par_T = phi_ox_PF_par_T-phi_red_PF_par_T;              % Extent of O atoms change in solid
    delta_phi_CF_par_T = phi_ox_CF_par_T-phi_red_CF_par_T;              % Extent of O atoms change in solid
    prod_ox_PF_par_T = delta_phi_PF_par_T.*1e6./(M_MO.*1e3);            % Oxidation productivity [micro-mole_prod/g_MO)
    prod_ox_CF_par_T = delta_phi_CF_par_T.*1e6./(M_MO.*1e3);            % Oxidation productivity [micro-mole_prod/g_MO)
    X_PF_par_T = min(delta_phi_PF_par_T./omega_ox,1);                   % Conversion extent - PF
    X_CF_par_T = min(delta_phi_CF_par_T./omega_ox,1);                   % Conversion extent - CF
    if MFR_flag==1
        delta_delta_MFR_par_T = delta_red_MFR_par_T-delta_ox_MFR_par_T;     % Extent of reduction (delta_red-delta_ox) for MFR oxidation
        delta_phi_MFR_par_T = phi_ox_MFR_par_T-phi_red_MFR_par_T;           % Extent of O atoms change in solid
        prod_ox_MFR_par_T = delta_phi_MFR_par_T.*1e6./(M_MO.*1e3);          % Oxidation productivity [micro-mole_prod/g_MO)
        X_MFR_par_T = min(delta_phi_MFR_par_T./omega_ox,1);                 % Conversion extent - MFR
    end
    % Clean unphysical results (delta_delta<0 and X<0)
    [row,col] = find(delta_delta_PF_par_T<0);
    delta_delta_PF_par_T(row,col) = NaN;
    [row,col] = find(delta_delta_CF_par_T<0);
    delta_delta_CF_par_T(row,col) = NaN;
    [row,col] = find(X_PF_par_T<0);
    X_PF_par_T(row,col) = NaN;
    [row,col] = find(X_CF_par_T<0);
    X_CF_par_T(row,col) = NaN;
    if MFR_flag==1
        [row,col] = find(delta_delta_MFR_par_T<0);
        delta_delta_MFR_par_T(row,col) = NaN;
        [row,col] = find(X_MFR_par_T<0);
        X_MFR_par_T(row,col) = NaN;
    end
end
%% Save specific run data
SAVE_SPEC_RUN = 0;
spec_run_filename = 'Prod_comp.xlsx';
if SAVE_SPEC_RUN == 1
    if MFR_flag==1
        Tab = table(T_red,p_red,omega_red,phi,T_ox,p_ox,omega_ox,x_p_in,delta_red_PF,delta_red_CF,delta_red_MFR,...
            delta_ox_PF,delta_ox_CF,delta_ox_MFR,delta_delta_PF,delta_delta_CF,delta_delta_MFR,delta_phi_PF,...
            delta_phi_CF,delta_phi_MFR,X_PF,X_CF,X_MFR,prod_red_PF,prod_red_CF,prod_red_MFR,prod_ox_PF,prod_ox_CF,...
            prod_ox_MFR,ox_str,string(MO_label));
    else
        Tab = table(T_red,p_red,omega_red,phi,T_ox,p_ox,omega_ox,x_p_in,delta_red_PF,delta_red_CF,...
            delta_ox_PF,delta_ox_CF,delta_delta_PF,delta_delta_CF,delta_phi_PF,...
            delta_phi_CF,X_PF,X_CF,prod_red_PF,prod_red_CF,prod_ox_PF,prod_ox_CF,...
            ox_str,string(MO_label));
    end
    if isfile(spec_run_filename)
        writetable(Tab,spec_run_filename,'WriteRowNames',false,'WriteVariableNames',false,'WriteMode','append');
    else
        writetable(Tab,spec_run_filename,'WriteRowNames',true);
    end
end
%% Save data
if ~isempty(filename_str)
    % Get a list of all variables
    allvars = whos;
    % Identify the variables that ARE NOT graphics handles. This uses a regular
    % expression on the class of each variable to check if it's a graphics object
    tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
    % Pass these variable names to save
    save(filename_str, allvars(tosave).name);
end
%% Comparison with experiments of Venstrom et al. (2014+2015) and Davenport et al. (2016)
% Venstrom
% v_N2_Ven_max = 0.6e-3/60*(1773/300);
% v_N2_Ven_min = 0.05e-3/60*(1773/300);
% v_CO2_Ven_max = 0.6e-3/60*(1773/300);
% v_CO2_Ven_min = 0.05e-3/60*(1773/300);
% m_MO_Ven = 1.0245e-3;
% V_MO_Ven = 0.8*10e-3*0.25*pi*(9.5e-3)^2;
% tau_min = V_MO_Ven/v_N2_Ven_max;
% tau_max = V_MO_Ven/v_N2_Ven_min;
% 600e-6/60*1e5/(R*300)
% Davenport