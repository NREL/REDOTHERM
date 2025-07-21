%% ---------------General info---------------
% Created by: Alon Lidor (alon.lidor@nrel.gov)
% National Renewable Energy Laboratory (NREL)
% Date: Dec. 11, 2024
% This code calculates the performance of a counter-current flow (CF) and
% parallel flow (PF) redox reactor using the method developed by Bulfin
% (2019) to avoid violations of the second law of thermodynamics based on
% chemical potential calculations. It extends the original model to include
% both reduction and oxidation, as well as enthalpy and entropy
% dependencies on nonstoichiometry (delta), and the option to analyze more
% oxide materials. It utilizes an optimization scheme. This is a
% non-interactable code that can perform multiple optimizaton runs,
% covering a wide range of cases
%--------------------------------------------------------------------------
%% Main
clearvars;
clc;
close all;
import py.CoolProp.CoolProp.*       % Load Python CoolProp package
% cd functions\;
addpath(genpath('./functions'));    % Adds the functions in the subfolder 'functions'
addpath(genpath('./materials'));    % Adds the functions in the subfolder 'materials'
tic
% Set optimization boundaries
omega_red_min = 0.01;   % Minimum molar flow rate ratio between sweep gas and MO [mol-N2/mol-MO]
omega_red_max = 100;    % Maximum molar flow rate ratio between sweep gas and MO [mol-N2/mol-MO]
omega_ox_min = 0.01;    % Minimum molar flow rate ratio between oxidizer gas and MO [mol-ox/mol-MO]
omega_ox_max = 100;     % Maximum molar flow rate ratio between oxidizer gas and MO [mol-ox/mol-MO]
T_red_min = 1673.15;    % Minimum reduction temperature [K]
T_red_max = 1973.15;    % Maximum reduction temperature [K]
T_ox_min = 873.15;     % Minimum oxidation temperature [K]
T_ox_max = 1473.15;     % Maximum oxidation temperature [K]
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
phi = 1e-5;                 % O2 impurity in sweep gas
T0 = 298.15;                % Ambient temperature [K]
red_config = 2;             % Reduction reactor configuration (1-PF, 2-CF)
ox_config = 2;              % 0.Oxidation reactor configuration (1-PF, 2-CF)
T_pump = 300;               % Vacuum pump temperature [K]
f_th_loss = 0.0;            % Thermal losses (reduction reactor)
% Initial Guess Input
p_red = [1e5];                  % Reduction pressure [Pa]
p_ox = [1e5];                   % Oxidation pressure [Pa]
T_red = T_red_max;              % Reduction temperature [K]
T_ox = 0.5*(T_ox_min+T_ox_max); % Oxidation temperature [K]
omega_red = 1;                  % Ratio of sweep gas molar flow rate to redox material flow rate [mol-sg/s to mol-redox/s]
omega_ox = 1;                   % Ratio of oxidizer gas molar flow rate to redox material flow rate [mol-sg/s to mol-redox/s]
% Oxidizer composition selection
x_ox_in = input('Enter the oxidizer mole fraction in the oxidizer stream (leave empty to use equilibrium composition at T_ox): ');
if isempty(x_ox_in)
    ox_comp_flag = 0;
else
    ox_comp_flag = 1;
end
% --- Input to sweep over ---
% Heat recovery
eps_HR = [0 0.25 0.5 0.75 1];       % Solid heat recovery effectiveness
% eps_HR = 0.5;
eps_HR_ox = [0.4 0.8];              % Exothermic oxidation heat recovery effectiveness
% eps_HR_ox = 0.8;              % Exothermic oxidation heat recovery effectiveness
eta_ox_htw = 0.4;                   % Efficiency of converting heat to work from excess exothermic heat
eps_g = 0.8;                        % Gas-gas heat recovery effectiveness
% Select redox material 1-CeO2,2-CeZr20,3-LCMA,4-LSM40,5-Fe33Al67
redox_material = [ 3 4];
% redox_material = [1];
% Select CO2 or H2O splitting (1-H2O, 2-CO2)
K_input = 1;
% Select H2-H2O separation technology (1-condensation,2-mechanical vapor recompression)
prod_sep_flag = [1];
% Choose optimizer
% 1 - fmincon
% 2 - Surrogate Optimization (surrogateopt) - NOT WORKING
% 3 - Pattern Search (patternsearch)
% 4 - Genetic Algorithm + fmincon (Hybrid)
% 5 - Global Search
% 6 - Global Search + Multi Start
% 7 - Pattern Search with parallel computing (patternsearch) - NO BENEFIT WITHOUT PARFOR
opt_flag = 3;
% Strings for filenames
MO_name = {'_CeO2','_CeZr20','_LCM','_LSM40','_Fe33Al67'};
Ox_name = {'_H2O','_CO2'};
Prod_sep_name = {'_cond','_MVR'};
num_runs = length(redox_material)*length(prod_sep_flag)*length(eps_HR)*length(eps_HR_ox);   % Number of total runs
%% Start looping over all parameters
for I=1:length(redox_material)
    for J=1:length(prod_sep_flag)
        for K=1:length(eps_HR)
            for L=1:length(eps_HR_ox)
                for M=1:length(p_red)
                    for N=1:length(p_ox)
                        % disp(['---Solving run'])
                        switch redox_material(I)
                            case 1
                                M_MO = 172.1e-3;                                            % Molar mass of ceria [kg/mol]
                                rho_MO = 7220;                                              % Solid density [kg/m^3]
                                % dH_fun = @(delta)Reduction_Enthalpy_CeO2(delta);            % Reduction enthalpy function handle
                                % dS_fun = @(delta)Reduction_Entropy_CeO2(delta);             % Reduction entropy function handle
                                dH_fun = @(delta)panlener_dhfun_mex(delta);            % Reduction enthalpy function handle
                                dS_fun = @(delta)panlener_dsfun_mex(delta);             % Reduction entropy function handle
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
                                rho_MO = 7220;                                                  % Solid density [kg/m^3]
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
                                rho_MO = 7000;                                              % Solid density [kg/m^3] - PLACEHOLDER VALUE FOR THIS MATERIAL
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
                                rho_MO = 7000;                                              % Solid density [kg/m^3] - PLACEHOLDER VALUE FOR THIS MATERIAL
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
                                M_MO = 0.33*M_Fe+0.67*M_Al+4*M_O/(3-1/3);                       % Molar mass of Fe33Al67 [kg/mol] - Fe_0.33-Al_0.67-O_4/(3-1/3)
                                rho_MO = 7000;                                                  % Solid density [kg/m^3] - PLACEHOLDER VALUE FOR THIS MATERIAL
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
                        pO2_fun = @(T,delta)pO2_calc_fun_generic(T,delta,dH_fun,dS_fun,dphi_fun);
                        pO2_der_fun = @(T,delta)pO2_der_calc_fun_generic(T,delta,dH_fun,dS_fun,dH_ddelta_fun,dS_ddelta_fun,dphi_fun);
                        %% Optimization
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
                        S.h_fg = CP_PropsSI('HMOLAR','P',p_ox(N),'Q',1,'H2O')-CP_PropsSI('HMOLAR','P',p_ox(N),'Q',0,'H2O');               % Enthalpy of vaporization at p_ox [J/mol-H2O]
                        S.Q_liq_heat = max(CP_PropsSI('HMOLAR','Q',0,'P',p_ox(N),'H2O')-CP_PropsSI('HMOLAR','T',T0,'P',p_ox(N),'H2O'),0); % Specific heat to heat liquid at T0 to Tsat [J/mol_H2O]
                        S.eta_CO2_sep = eta_CO2_sep;                    % CO2 separation energy [J/mol-CO2]
                        S.E_inert = E_inert;                            % Inert gas separation energy [J/mol-N2]
                        S.cp_s_fun = cp_s_fun;                          % Specific heat capacity of MO function handle [J/mol-MO]
                        S.M_MO = M_MO;                                  % MO molar mass [kg/mol]
                        S.dH_fun = dH_fun;                              % MO reduction enthalpy function handle [J/mol-O]
                        S.eps_HR = eps_HR(K);                           % MO sensible heat recovery effectiveness
                        S.eps_HR_ox = eps_HR_ox(L);                     % Exothermic oxidation heat recovery effectiveness
                        S.eta_ox_htw = eta_ox_htw;                      % Efficiency of converting excess exothermic heat to work
                        S.eta_pump = eta_pump;                          % Efficiency of pumping oxidizer to higher pressure
                        S.eps_g = eps_g;                                % Gas-gas heat recovery effectiveness
                        S.p_red = p_red(M);                             % Reduction pressure [Pa]
                        S.p_ox = p_ox(N);                               % Oxidation pressure [Pa]
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
                        S.prod_sep_flag = prod_sep_flag(J);             % Product separation technology flag
                        % Optimization
                        X0 = [T_red,T_ox,p_red(M),p_ox(N),omega_red,omega_ox];                % Set initial values
                        [X_sol,fvalue,exitflag] = OptimizeRedoxCycle(X0,S,opt_flag);    % Solve optimization problem
                        Sol = Analyze_Redox_Cycle(X_sol,S);                             % Get full solution structure for optimal solution
                        % Save data
                        filename_str = strcat(string(datetime('today','Format','uuuuMMdd')),MO_name(redox_material(I)),Ox_name(K_input),Prod_sep_name(prod_sep_flag(J)),'_epsHRs',num2str(eps_HR(K)*100),'_epsHRox',num2str(eps_HR_ox(L)*100),'_p_red',num2str(p_red(M)/1e5),'_p_ox',num2str(p_ox(N)/1e5));
                        save(filename_str)
                    end
                end
            end
        end
    end
end
%% Save data
% if ~isempty(filename_str)
%     % Get a list of all variables
%     allvars = whos;
%     % Identify the variables that ARE NOT graphics handles. This uses a regular
%     % expression on the class of each variable to check if it's a graphics object
%     tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
%     % Pass these variable names to save
%     save(filename_str, allvars(tosave).name);
% end
%% Functions
function result = CP_PropsSI(varargin)
    % Shorthand version of CoolProp for MATLAB
    result = py.CoolProp.CoolProp.PropsSI(varargin{:});
end