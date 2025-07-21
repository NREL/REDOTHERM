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
% oxide materials. It utilizes an optimization scheme
%--------------------------------------------------------------------------
%% Main
clearvars;
clc;
close all;
import py.CoolProp.CoolProp.*       % Load Python CoolProp package
% cd functions\;
addpath(genpath('./functions'));    % Adds the functions in the subfolder 'functions'
addpath(genpath('./materials'));    % Adds the functions in the subfolder 'materials'
addpath(genpath('./core'));         % Adds the functions in the subfolder 'core'
% Plotting defaults
% Settings
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultTextFontSize',14);
fontsize_l = 13;
fontsize_a = 12;
% Save options
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
omega_red_min = input("Enter minimum flow rate ratio of sweep gas to oxide [omega_red_min=0.001]: ");
if isempty(omega_red_min)
    omega_red_min = 0.001;
end
omega_red_max = input("Enter maximum flow rate ratio of sweep gas to oxide [omega_red_max=1000]: ");
if isempty(omega_red_max)
    omega_red_max = 1000;
end
omega_ox_min = input("Enter minimum flow rate ratio of oxidizer to oxide [omega_ox_min=0.001]: ");
if isempty(omega_ox_min)
    omega_ox_min = 0.001;
end
omega_ox_max = input("Enter maximum flow rate ratio of oxidizer to oxide [omega_ox_max=1000]: ");
if isempty(omega_ox_max)
    omega_ox_max = 1000;
end
T_red_min = input("Enter minimum reduction temperature in deg. C [Tred_min=1300]: ");
if isempty(T_red_min)
    T_red_min = 1573.15;
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
% Fixed input
eta_CO2_sep = 0.1;          % CO-CO2 separation efficiency
eta_pump = 0.9;             % Efficiency of increasing the pressure of the oxidizier (cold) - pump/compressor
E_inert = 15e3;             % Cryogenic N2 separation energy [J/mol-N2]
eps_HR = 0.5;               % Solid heat recovery effectiveness
eps_HR_ox = 0.8;            % Exothermic oxidation heat recovery effectiveness
eta_ox_htw = 0.4;           % Efficiency of converting heat to work from excess exothermic heat
eps_g = 0.8;                % Gas-gas heat recovery effectiveness
T0 = 298.15;                % Ambient temperature [K]
red_config = 2;             % Reduction reactor configuration (1-PF, 2-CF)
ox_config = 2;              % Oxidation reactor configuration (1-PF, 2-CF)
T_pump = 300;               % Vacuum pump temperature [K]
f_th_loss = 0.0;            % Thermal losses (reduction reactor)
% Data from validation case in Bulfin (2019) (taken from Scheffe et al.
% (2014)
% omega_red_BB = (1.487E-05/39.95e-3)/(61e-6/172.1e-3);
% Input
p_red = input("Enter reduction pressure (p_red) in Pa [p_red=100000]: ");
if isempty(p_red)
    p_red = 1e5;                % Reduction pressure [Pa]
end
p_ox = input("Enter oxidation pressure (p_ox) in Pa [p_ox=100000]: ");
if isempty(p_ox)
    p_ox = 1e5;                 % Oxidation pressure [Pa]
end
% Reduction temperature [K]
T_red = input("Enter reduction temperature (T_red) in deg. C [T_red=1600]: ");
if isempty(T_red)
    T_red = 1873.15;
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
        dH_fun = @(delta)panlener_dhfun_mex(delta);                 % Reduction enthalpy function handle
        dS_fun = @(delta)panlener_dsfun_mex(delta);                 % Reduction entropy function handle
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
        M_MO = 0.33*M_Fe+0.67*M_Al+4*M_O;                               % Molar mass of Fe33Al67 [kg/mol] - Fe0.33Al0.67.O4
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
% Equilibrium constants
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
% Select H2-H2O separation technology
prod_sep_flag = input('Select H2-H2O separation technology (1-condensation,2-mechanical vapor recompression) [1]: ');
if isempty(prod_sep_flag)
    prod_sep_flag = 1;
end
% Select if saving the data or not
filename_str = input('Enter data file name (leave empty to skip saving): ',"s");
%% Reduction
disp('----- Reduction -----');
% Solve for PF
nO2_max_PF = nO2_max_PF_red_fun(T_red,p_red,phi,omega_red,nO2_total,pO2_fun,delta_fun,phi0);    % Calculate maximum reaction extent (O2 transfer) for PF reduction
phi_red_PF = phi0-nO2_max_PF*2;                                                                 % Calculate maximum phi for PF reduction
delta_red_PF = delta_fun(phi_red_PF);                                                           % Calculate reduction delta for PF
% Solve for CF
nO2_max_CF = nO2_max_CF_red_fun(T_red,p_red,phi,omega_red,nO2_total,pO2_fun,pO2_der_fun,delta_fun,phi0);    % Calculate maximum reaction extent (O2 transfer) for CF reduction
phi_red_CF = phi0-nO2_max_CF*2;                                                                             % Calculate maximum phi reduction
delta_red_CF = delta_fun(phi_red_CF);                                                                       % Calculate reduction delta for CF
disp(['The maximum O2 exchange for parallel flow (PF) reduction is: ',num2str(nO2_max_PF)]);
disp(['The maximum O2 exchange for countercurrent flow (CF) reduction is: ',num2str(nO2_max_CF)]);
disp(['The maximum reduction extent (delta) for parallel flow (PF) reduction is: ',num2str(delta_red_PF)]);
disp(['The maximum reduction extent (delta) for countercurrent flow (CF) reduction is: ',num2str(delta_red_CF)]);
%% Oxidation
disp('----- Oxidation -----');
% Inlet oxidier composition - either equilibirium of fixed
gas = GRI30;
nsp = nSpecies(gas);
iCO = speciesIndex(gas,'CO');
iCO2 = speciesIndex(gas,'CO2');
iO2 = speciesIndex(gas,'O2');
iH2 = speciesIndex(gas,'H2');
iH2O = speciesIndex(gas,'H2O');
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
else
    xH2_in = 1-x_ox_in;     % Inlet H2 mole fraction
    xH2O_in = x_ox_in;      % Inlet H2O mole fraction
    xO2_WS_in = 0;          % Inlet O2 mole fraction
    xCO_in = 1-x_ox_in;     % Inlet CO mole fraction
    xCO2_in = x_ox_in;      % Inlet CO2 mole fraction
    xO2_CDS_in = 0;         % Inlet O2 mole fraction
end
% Calculate partial pressures
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
% HHV -     higher heating value of product [J/mol]
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
% Solve for PF
nO2_max_PF_ox = nO2_max_PF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox,nO2_max_PF,pO2_fun,K,delta_fun,phi_red_PF);    % Calculate maximum reaction extent (O2 transfer) for PF oxidation
phi_ox_PF = phi_red_PF+2*nO2_max_PF_ox;                                                                     % Final (minimum) phi for PF oxidation
delta_phi_PF = phi_ox_PF-phi_red_PF;                                                                        % Change in moles of O over moles of MO [mol-O/mol-MO]
delta_ox_PF = delta_fun(phi_ox_PF);                                                                         % Oxidation delta for PF
delta_delta_PF = (delta_red_PF-delta_ox_PF)*((1-red_mode)-red_mode);                                        % Extent of reduction (delta_red-delta_ox) for PF oxidation
X_PF = min(delta_phi_PF/omega_ox,1);                                                                        % Conversion extent (PF)
% Solve for CF
nO2_max_CF_ox = nO2_max_CF_ox_fun(T_ox,p_ox,x_r_in,x_p_in,omega_ox,nO2_max_CF,pO2_fun,pO2_der_fun,K,delta_fun,phi_red_CF);   % Calculate maximum reaction extent (O2 transfer) for CF oxidation
phi_ox_CF = phi_red_CF+2*nO2_max_CF_ox;                                                                                 % Final (minimum) nonstoichiometry extent for CF oxidation
delta_phi_CF = phi_ox_CF-phi_red_CF;                                                                                    % Change in moles of O over moles of MO [mol-O/mol-MO]
delta_ox_CF = delta_fun(phi_ox_CF);                                                                                     % Calculate phi for CF oxidation
delta_delta_CF = (delta_red_CF-delta_ox_CF)*((1-red_mode)-red_mode);                                                    % Extent of reduction (delta_red-delta_ox) for CF oxidation
X_CF = min(delta_phi_CF/omega_ox,1);                                                                                    % Conversion extent (CF)
disp(['The product mole fraction in the oxidizer inlet stream is: ',num2str(x_p_in)]);
disp(['The maximum O2 exchange for parallel flow (PF) oxidation is: ',num2str(nO2_max_PF_ox)]);
disp(['The maximum O2 exchange for countercurrent flow (CF) oxidation is: ',num2str(nO2_max_CF_ox)]);
disp(['The minimum reduction extent (delta) for parallel flow (PF) oxidation is: ',num2str(delta_ox_PF)]);
disp(['The minimum reduction extent (delta) for countercurrent flow (CF) oxidation is: ',num2str(delta_ox_CF)]);
disp(['The conversion extent (X) for parallel flow (PF) oxidation is: ',num2str(X_PF)]);
disp(['The conversion extent (X) for countercurrent flow (CF) oxidation is: ',num2str(X_CF)]);
disp('----- Cycle Performance -----');
disp(['The extent of reduction for parallel flow (PF) redox cycle is: ',num2str(delta_delta_PF)]);
disp(['The extent of reduction for countercurrent flow (CF) redox cycle is: ',num2str(delta_delta_CF)]);
%% Optimization
tic
% Select optimization variable target (1-delta_delta,2-conversion,3-efficiency)
opt_target = 3;
% Create structure with auxiliary parameters
S.phi = phi;                    % O2 mole fraction in sweep gas
S.nO2_total = nO2_total;        % Maximum specific O2 release per mole of redox material [mol-O2/mol-redox]
S.pO2_fun = pO2_fun;            % MO partial pressure function handle [Pa]
S.pO2_der_fun = pO2_der_fun;    % MO partial pressure derivative function handle [Pa]
S.red_mode = red_mode;          % Material type reduction mode
S.delta_fun = delta_fun;        % Function handle of delta(phi)
S.delta0 = delta0;              % Initial nonstoichiometry extent at fully oxidized state
S.phi0 = phi0;                  % Initial moles of O atoms in MO per moles of MO in fully oxidized state
S.ox_comp_flag = ox_comp_flag;  % Oxidizer composition flag (0-equilibrium,1-fixed)
S.x_ox_in = x_ox_in;            % Mole fraction of oxidizer gas in feed during oxidation (i.e. xH2O_in or xCO2_in)
S.K_input = K_input;            % Type of oxidation reaction (1-WS,2-CDS)
S.red_config = red_config;      % Reduction reactor configuration (1-PF,2-CF)
S.ox_config = ox_config;        % Oxidation reactor configuration (1-PF,2-CF)
S.h_fg = CP_PropsSI('HMOLAR','P',p_ox,'Q',1,'H2O')-py.CoolProp.CoolProp.PropsSI('HMOLAR','P',p_ox,'Q',0,'H2O');     % Enthalpy of vaporization at p_ox [J/mol-H2O]
S.Q_liq_heat = max(CP_PropsSI('HMOLAR','Q',0,'P',p_ox,'H2O')-CP_PropsSI('HMOLAR','T',T0,'P',p_ox,'H2O'),0);         % Specific heat to heat liquid at T0 to Tsat [J/mol_H2O]
S.eta_pump = eta_pump;                          % Oxidizier pumping/compression efficiency
S.eta_CO2_sep = eta_CO2_sep;                    % CO2 separation energy [J/mol-CO2]
S.E_inert = E_inert;                            % Inert gas separation energy [J/mol-N2]
S.cp_s_fun = cp_s_fun;                          % Specific heat capacity of MO function handle [J/mol-MO]
S.M_MO = M_MO;                                  % MO molar mass [kg/mol]
S.dH_fun = dH_fun;                              % MO reduction enthalpy function handle [J/mol-O]
S.oxidizer_type = ox_str;                       % Oxidizer name
S.eps_HR = eps_HR;                              % MO sensible heat recovery effectiveness
S.eps_HR_ox = eps_HR_ox;                        % Exothermic oxidation heat recovery effectiveness
S.eta_ox_htw = eta_ox_htw;                      % Efficiency of converting excess exothermic heat to work
S.eps_g = eps_g;                                % Gas-gas heat recovery effectiveness
S.p_red = p_red;                                % Reduction pressure [Pa]
S.p_ox = p_ox;                                  % Oxidation pressure [Pa]
S.T_red_min = T_red_min;                        % Minimum reduction temperature [K]
S.T_red_max = T_red_max;                        % Maximum reduction temperature [K]
S.T_ox_min = T_ox_min;                          % Minimum oxidation temperature [K]
S.T_ox_max = T_ox_max;                          % Maximum oxidation temperature [K]
S.T0 = T0;                                      % Ambient temperature [K]
S.T_pump = T_pump;                              % Vacuum pump temperature [K]
S.f_th_loss = f_th_loss;                        % Thermal losses fraction (reduction reactor)
S.omega_red_min = omega_red_min;                % Minimum omega_red [mol-N2/s / mol-MO/s]
S.omega_red_max = omega_red_max;                % Maximum omega_red [mol-N2/s / mol-MO/s]
S.omega_ox_min = omega_ox_min;                  % Minimum omega_ox [mol-ox/s / mol-MO/s]
S.omega_ox_max = omega_ox_max;                  % Maximum omega_ox [mol-ox/s / mol-MO/s]
S.prod_sep_flag = prod_sep_flag;                % Product separation technology flag
opt_flag = 3;                                   % Choose optimizer
% Optimization
X0 = [T_red,T_ox,p_red,p_ox,omega_red,omega_ox];                % Set initial values
[X_sol,fvalue,exitflag] = OptimizeRedoxCycle(X0,S,opt_flag);    % Solve optimization problem
Sol = Analyze_Redox_Cycle(X_sol,S);                             % Get full solution structure for optimal solution
runtime = toc;
% Display results
disp('---- Optimal Solution ----');
disp(['Total runtime: ',num2str(runtime),' s']);
disp(['System efficiency: ',num2str(Sol.eta)]);
disp(['Feedstock conversion extent: ',num2str(Sol.X_conv)]);
disp(['Reduction temperature: ',num2str(X_sol(1)),' K']);
disp(['Oxidation temperature: ',num2str(X_sol(2)),' K']);
disp(['Molar flow rate of sweep gas to redox material: ',num2str(X_sol(5)),' mol-N2/mol_MO']);
disp(['Molar flow rate of oxidizer gas to redox material: ',num2str(X_sol(6)),' mol-ox/mol_MO']);
% Get back to main directory
% cd ..
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
%% Functions
function result = CP_PropsSI(varargin)
    % Shorthand version of CoolProp for MATLAB
    result = py.CoolProp.CoolProp.PropsSI(varargin{:});
end