%% ---------------General info---------------
% Created by: Alon Lidor (alon.lidor@nrel.gov)
% National Renewable Energy Laboratory (NREL)
% Date: June 1, 2025
% This example script demonstrate results analysis from multiple
% optimization runs using the REDOTHERM model. This example code is
% structured to load all results file that are in a single folder, and
% analyze energy terms as well as performance indicators such as efficiency
% and conversion, and compare results for different materials and H2
% separation technologies. This code can easily be adapted to perform
% different analysis
%--------------------------------------------------------------------------
%% Initialization
clearvars
close all
clc
% Get folder with all results file
selpath = uigetdir;
cd(selpath);
d = dir('*.mat');
fn = {d.name};
[indx,tf] = listdlg('PromptString','Select file to include in plotting:','SelectionMode','multiple','ListString',fn,'ListSize',[350 500]);
len = length(indx);
Tab = table;
eta = zeros(len,1);
X_conv = zeros(len,1);
% Load selected results files
for Ind=1:len
    LoadS = load(d(indx(Ind)).name);
    Temp_tab = struct2table(LoadS.Sol);
    MO_name = regexprep(string(LoadS.MO_name(LoadS.redox_material(LoadS.I))),'_','');
    Ox_name = regexprep(string(LoadS.Ox_name(LoadS.K_input)),'_','');
    ProdSepName = regexprep(string(LoadS.Prod_sep_name(LoadS.prod_sep_flag(LoadS.J))),'_','');
    Temp_tab2 = array2table(LoadS.X_sol,'VariableNames',{'T_red','T_ox','p_red','p_ox','omega_red','omega_ox'});
    Temp_tab2 = addvars(Temp_tab2,MO_name,Ox_name,ProdSepName,LoadS.T0,LoadS.S.eps_HR,LoadS.S.eps_HR_ox,LoadS.S.eps_g,LoadS.eta_ox_htw,'Before',"T_red",'NewVariableNames',{'Redox_material','Oxidizer','Prod_sep','T0','eps_HR','eps_HR_ox','eps_g','eta_ox_htw'});
    Tab = [Tab ; [Temp_tab2 Temp_tab]]; %#ok<AGROW>
end
cd ..;
% Calculate F_total (total energy fraction)
F_total = (Tab.F_ox_sep+Tab.F_inert_sep+Tab.F_ox_pump+Tab.F_vac+Tab.F_work_rec)+Tab.F_sens_MO+Tab.F_red+Tab.F_inert_h+max(Tab.F_ox_h+Tab.F_ox_exo,0);
Tab = addvars(Tab,F_total,'After',"F_work_rec",'NewVariableNames','F_total');
% Calculate specific energy/heat terms in energy per mole of product
w_s_ox_sep = Tab.F_ox_sep.*Tab.HHV;          % Effluent stream product separation energy per mole of product [J/mol-H2/CO]
q_s_ox_h = Tab.F_ox_h.*Tab.HHV;              % Effluent stream product separation energy per mole of product [J/mol-H2/CO]
w_s_inert = Tab.F_inert_sep.*Tab.HHV;        % Inert gas separation energy per mole of product [J/mol-H2/CO]
q_s_sens_MO = Tab.F_sens_MO.*Tab.HHV;        % Sensible heat per mole of product [J/mol-H2/CO]
q_s_red = Tab.F_red.*Tab.HHV;                % Reduction energy per mole of product [J/mol-H2/CO]
q_s_ox_exo = Tab.F_ox_exo.*Tab.HHV;          % Exothermic heat per mole of product [J/mol-H2/CO]
q_s_inert_h = Tab.F_inert_h.*Tab.HHV;        % Inert sweep gas heating per mole of product [J/mol-H2/CO]
w_s_work_rec = Tab.F_work_rec.*Tab.HHV;      % Work recovered from oxidation heat per mole of product [J/mol-H2/CO]
q_t_s = Tab.F_total.*Tab.HHV;                % Total energy per mole of product [J/mol-H2/CO]
Tab = addvars(Tab,w_s_ox_sep,q_s_ox_h,w_s_inert,q_s_sens_MO,q_s_red,q_s_ox_exo,q_s_inert_h,w_s_work_rec,q_t_s,...
        'After',"Q_inert_h",'NewVariableNames',{'w_s_ox_sep','q_s_ox_h','w_s_inert','q_s_sens_MO','q_s_red','q_s_ox_exo','q_s_inert_h','w_s_work_rec','q_t_s'});
% Tab = [Tab ; [Temp_tab2  Temp_tab]];
%% Post processing
Tab = sortrows(Tab,{'Redox_material','Oxidizer','Prod_sep','eps_HR'},"ascend");                     % Reorder results table
indx_baseline_cond = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));              % Get indices of baseline case using condenser+boiler
indx_baseline_MVR = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));                % Get indices of baseline case using MVR separation
indx_no_solid_HR_cond = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));             % Get indices of no solid HR case using condenser+boiler
indx_no_solid_HR_MVR = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));               % Get indices of no solid HR case using MVR separation
indx_full_solid_HR_cond = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));           % Get indices of full solid HR case using condenser+boiler
indx_full_solid_HR_MVR = (Tab.p_red==1e5)&(Tab.p_ox==1e5)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));             % Get indices of full solid HR case using MVR separation
indx_baseline_cond_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));     % Get indices of baseline case using condenser+boiler
indx_baseline_MVR_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));       % Get indices of baseline case using MVR separation
indx_no_solid_HR_cond_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));    % Get indices of no solid HR case using condenser+boiler
indx_no_solid_HR_MVR_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));      % Get indices of no solid HR case using MVR separation
indx_full_solid_HR_cond_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));  % Get indices of full solid HR case using condenser+boiler
indx_full_solid_HR_MVR_int_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==5e5)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));    % Get indices of full solid HR case using MVR separation
indx_baseline_cond_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));    % Get indices of baseline case using condenser+boiler
indx_baseline_MVR_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==0.5)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));      % Get indices of baseline case using MVR separation
indx_no_solid_HR_cond_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond'));   % Get indices of no solid HR case using condenser+boiler
indx_no_solid_HR_MVR_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==0)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));     % Get indices of no solid HR case using MVR separation
indx_full_solid_HR_cond_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'cond')); % Get indices of full solid HR case using condenser+boiler
indx_full_solid_HR_MVR_high_p_ox = (Tab.p_red==1e5)&(Tab.p_ox==1e6)&(Tab.eps_HR==1)&(Tab.eps_HR_ox==0.8)&(strcmp(string(Tab.Prod_sep),'MVR'));   % Get indices of full solid HR case using MVR separation
indx_all_cases = indx_baseline_cond|indx_baseline_MVR|indx_no_solid_HR_cond|indx_no_solid_HR_MVR|indx_full_solid_HR_cond|indx_full_solid_HR_MVR;    % Get indices of all cases together
%% Plotting
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
% Get baseline values - condenser case
eta_baseline_cond_array = Tab.eta(indx_baseline_cond,:);                            % System efficiency
X_conv_baseline_cond_array = Tab.X_conv(indx_baseline_cond,:);                      % Conversion extent
T_red_baseline_cond_array = Tab.T_red(indx_baseline_cond,:);                        % Reduction temperature [K]
T_ox_baseline_cond_array = Tab.T_ox(indx_baseline_cond,:);                          % Oxidation temperature [K]
Delta_T_baseline_cond_array = T_red_baseline_cond_array-T_ox_baseline_cond_array;   % Temperature swing [K]
omega_red_baseline_cond_array = Tab.omega_red(indx_baseline_cond,:);                % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_baseline_cond_array = Tab.omega_ox(indx_baseline_cond,:);                  % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_baseline_cond_array = Tab.delta_red(indx_baseline_cond,:);                % Reduction extent - end of reduction
delta_ox_baseline_cond_array = Tab.delta_ox(indx_baseline_cond,:);                  % Reduction extent - end of oxidation
delta_phi_baseline_cond_array = Tab.delta_phi(indx_baseline_cond,:);                % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get no solid HR values - condenser case
eta_no_solid_HR_cond_array = Tab.eta(indx_no_solid_HR_cond,:);                              % System efficiency
X_conv_no_solid_HR_cond_array = Tab.X_conv(indx_no_solid_HR_cond,:);                        % Conversion extent
T_red_no_solid_HR_cond_array = Tab.T_red(indx_no_solid_HR_cond,:);                          % Reduction temperature [K]
T_ox_no_solid_HR_cond_array = Tab.T_ox(indx_no_solid_HR_cond,:);                            % Oxidation temperature [K]
Delta_T_no_solid_HR_cond_array = T_red_no_solid_HR_cond_array-T_ox_no_solid_HR_cond_array;  % Temperature swing [K]
omega_red_no_solid_HR_cond_array = Tab.omega_red(indx_no_solid_HR_cond,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_no_solid_HR_cond_array = Tab.omega_ox(indx_no_solid_HR_cond,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_no_solid_HR_cond_array = Tab.delta_red(indx_no_solid_HR_cond,:);                  % Reduction extent - end of reduction
delta_ox_no_solid_HR_cond_array = Tab.delta_ox(indx_no_solid_HR_cond,:);                    % Reduction extent - end of oxidation
delta_phi_no_solid_HR_cond_array = Tab.delta_phi(indx_no_solid_HR_cond,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get full solid HR values - condenser case
eta_full_solid_HR_cond_array = Tab.eta(indx_full_solid_HR_cond,:);                              % System efficiency
X_conv_full_solid_HR_cond_array = Tab.X_conv(indx_full_solid_HR_cond,:);                        % Conversion extent
T_red_full_solid_HR_cond_array = Tab.T_red(indx_full_solid_HR_cond,:);                          % Reduction temperature [K]
T_ox_full_solid_HR_cond_array = Tab.T_ox(indx_full_solid_HR_cond,:);                            % Oxidation temperature [K]
Delta_T_full_solid_HR_cond_array = T_red_full_solid_HR_cond_array-T_ox_full_solid_HR_cond_array;% Temperature swing [K]
omega_red_full_solid_HR_cond_array = Tab.omega_red(indx_full_solid_HR_cond,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_full_solid_HR_cond_array = Tab.omega_ox(indx_full_solid_HR_cond,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_full_solid_HR_cond_array = Tab.delta_red(indx_full_solid_HR_cond,:);                  % Reduction extent - end of reduction
delta_ox_full_solid_HR_cond_array = Tab.delta_ox(indx_full_solid_HR_cond,:);                    % Reduction extent - end of oxidation
delta_phi_full_solid_HR_cond_array = Tab.delta_phi(indx_full_solid_HR_cond,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Plot 1 - "baseline" conditions efficiency and production rate ratio
fig = figure; 
hold on
nil = zeros(length(eta_baseline_cond_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_cond,:),[eta_baseline_cond_array nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_cond,:),[nil X_conv_baseline_cond_array]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta','X_{ox}');
if SAVEFLAG
    filename = 'Baseline_cond_eta_X';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 2 - "baseline" conditions energy terms (absolute)
fig = figure; 
w_s_ox_sep_array = Tab.w_s_ox_sep(indx_baseline_cond,:).*1e-3;      % Product separation specific energy [kJ/mol-MO]
q_s_ox_h_array = Tab.q_s_ox_h(indx_baseline_cond,:).*1e-3;          % Oxidizer gas heating specific energy [kJ/mol-MO]
w_s_inert_array = Tab.w_s_inert(indx_baseline_cond,:).*1e-3;        % Inert sweep gas separation specific energy [kJ/mol-MO]
q_s_ox_exo_array = Tab.q_s_ox_exo(indx_baseline_cond,:).*1e-3;      % Exothermic heat of oxidation specific energy [kJ/mol-MO]
q_s_red_array = Tab.q_s_red(indx_baseline_cond,:).*1e-3;            % Endothermic reduction specific energy [kJ/mol-MO]
q_s_inert_h_array = Tab.q_s_inert_h(indx_baseline_cond,:).*1e-3;    % Inert sweep gas heating specific energy [kJ/mol-MO]
w_s_work_rec_array = Tab.w_s_work_rec(indx_baseline_cond,:).*1e-3;  % Recovered heat converetd to work specific energy [kJ/mol-MO]
q_s_sens_MO_array =Tab.q_s_sens_MO(indx_baseline_cond,:).*1e-3;     % MO sensible heating required specific energy [kJ/mol-MO]
q_t_s_array = Tab.q_t_s(indx_baseline_cond,:).*1e-3;                % Total required specific energy [kJ/mol-MO]
bar(Tab.Redox_material(indx_baseline_cond,:),[w_s_ox_sep_array q_s_ox_h_array w_s_inert_array q_s_inert_h_array q_s_sens_MO_array q_s_red_array -q_s_ox_exo_array]);
ylabel('Energy per produced H_2, kJ/mol_{H_2}');
legend({'w_{ox,sep}','q_{ox,h}','w_{sg,sep}','q_{sg,h}','q_{MO,h}','q_{red}','q_{ox}'},'NumColumns',2);
yscale('log');
if SAVEFLAG
    filename = 'Baseline_cond_energy_absolute';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 3 - "baseline" conditions energy terms (relative)
fig = figure;   
% Get array of selected cases
F_s_ox_sep_array = Tab.F_ox_sep(indx_baseline_cond,:);      % Product separation fraction
F_s_ox_h_array = Tab.F_ox_h(indx_baseline_cond,:);          % Oxidizer gas heating fraction
F_s_inert_array = Tab.F_inert_sep(indx_baseline_cond,:);    % Inert sweep gas separation fraction
F_s_ox_exo_array = Tab.F_ox_exo(indx_baseline_cond,:);      % Exothermic heat of oxidation fraction
F_s_red_array = Tab.F_red(indx_baseline_cond,:);            % Endothermic reduction fraction
F_s_inert_h_array = Tab.F_inert_h(indx_baseline_cond,:);    % Inert sweep gas heating fraction
F_s_work_rec_array = Tab.F_work_rec(indx_baseline_cond,:);  % Recovered heat converetd to work fraction
F_s_sens_MO_array =Tab.F_sens_MO(indx_baseline_cond,:);     % MO sensible heating required fraction
F_t_s_array = Tab.F_total(indx_baseline_cond,:);            % Total required fraction
b = bar(Tab.Redox_material(indx_baseline_cond,:),[F_s_ox_sep_array F_s_ox_h_array F_s_inert_array F_s_inert_h_array F_s_sens_MO_array F_s_red_array F_s_ox_exo_array]);
ylabel('Energy term ratio to fuel energy');
legend({'w_{ox,sep}','q_{ox,h}','w_{sg,sep}','q_{sg,h}','q_{MO,h}','q_{red}','q_{ox}'});
if SAVEFLAG
    filename = 'Baseline_cond_energy_relative';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Get baseline values - MVR case
eta_baseline_MVR_array = Tab.eta(indx_baseline_MVR,:);
X_conv_baseline_MVR_array = Tab.X_conv(indx_baseline_MVR,:);
T_red_baseline_MVR_array = Tab.T_red(indx_baseline_MVR,:);
T_ox_baseline_MVR_array = Tab.T_ox(indx_baseline_MVR,:);
Delta_T_baseline_MVR_array = T_red_baseline_MVR_array-T_ox_baseline_MVR_array;
omega_red_baseline_MVR_array = Tab.omega_red(indx_baseline_MVR,:);
omega_ox_baseline_MVR_array = Tab.omega_ox(indx_baseline_MVR,:);
delta_red_baseline_MVR_array = Tab.delta_red(indx_baseline_MVR,:);
delta_ox_baseline_MVR_array = Tab.delta_ox(indx_baseline_MVR,:);
delta_phi_baseline_MVR_array = Tab.delta_phi(indx_baseline_MVR,:);
% Get no solid HR values - MVR case
eta_no_solid_HR_MVR_array = Tab.eta(indx_no_solid_HR_MVR,:);
X_conv_no_solid_HR_MVR_array = Tab.X_conv(indx_no_solid_HR_MVR,:);
T_red_no_solid_HR_MVR_array = Tab.T_red(indx_no_solid_HR_MVR,:);
T_ox_no_solid_HR_MVR_array = Tab.T_ox(indx_no_solid_HR_MVR,:);
Delta_T_no_solid_HR_MVR_array = T_red_no_solid_HR_MVR_array-T_ox_no_solid_HR_MVR_array;
omega_red_no_solid_HR_MVR_array = Tab.omega_red(indx_no_solid_HR_MVR,:);
omega_ox_no_solid_HR_MVR_array = Tab.omega_ox(indx_no_solid_HR_MVR,:);
delta_red_no_solid_HR_MVR_array = Tab.delta_red(indx_no_solid_HR_MVR,:);
delta_ox_no_solid_HR_MVR_array = Tab.delta_ox(indx_no_solid_HR_MVR,:);
delta_phi_no_solid_HR_MVR_array = Tab.delta_phi(indx_no_solid_HR_MVR,:);
% Get full solid HR values - MVR case
eta_full_solid_HR_MVR_array = Tab.eta(indx_full_solid_HR_MVR,:);
X_conv_full_solid_HR_MVR_array = Tab.X_conv(indx_full_solid_HR_MVR,:);
T_red_full_solid_HR_MVR_array = Tab.T_red(indx_full_solid_HR_MVR,:);
T_ox_full_solid_HR_MVR_array = Tab.T_ox(indx_full_solid_HR_MVR,:);
Delta_T_full_solid_HR_MVR_array = T_red_full_solid_HR_MVR_array-T_ox_full_solid_HR_MVR_array;
omega_red_full_solid_HR_MVR_array = Tab.omega_red(indx_full_solid_HR_MVR,:);
omega_ox_full_solid_HR_MVR_array = Tab.omega_ox(indx_full_solid_HR_MVR,:);
delta_red_full_solid_HR_MVR_array = Tab.delta_red(indx_full_solid_HR_MVR,:);
delta_ox_full_solid_HR_MVR_array = Tab.delta_ox(indx_full_solid_HR_MVR,:);
delta_phi_full_solid_HR_MVR_array = Tab.delta_phi(indx_full_solid_HR_MVR,:);
% Plot 4 - "baseline" conditions efficiency and production rate ratio
fig = figure;   
hold on
nil = zeros(length(eta_baseline_MVR_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[eta_baseline_MVR_array nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[nil X_conv_baseline_MVR_array]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta','X_{ox}');
if SAVEFLAG
    filename = 'Baseline_MVR_eta_X';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 5 - "baseline" conditions energy terms (absolute)
fig = figure;   
w_s_ox_sep_array = Tab.w_s_ox_sep(indx_baseline_MVR,:).*1e-3;       % Product separation specific energy [kJ/mol-MO]
q_s_ox_h_array = Tab.q_s_ox_h(indx_baseline_MVR,:).*1e-3;           % Oxidizer gas heating specific energy [kJ/mol-MO]
w_s_inert_array = Tab.w_s_inert(indx_baseline_MVR,:).*1e-3;         % Inert sweep gas separation specific energy [kJ/mol-MO]
q_s_ox_exo_array = Tab.q_s_ox_exo(indx_baseline_MVR,:).*1e-3;       % Exothermic heat of oxidation specific energy [kJ/mol-MO]
q_s_red_array = Tab.q_s_red(indx_baseline_MVR,:).*1e-3;             % Endothermic reduction specific energy [kJ/mol-MO]
q_s_inert_h_array = Tab.q_s_inert_h(indx_baseline_MVR,:).*1e-3;     % Inert sweep gas heating specific energy [kJ/mol-MO]
w_s_work_rec_array = Tab.w_s_work_rec(indx_baseline_MVR,:).*1e-3;   % Recovered heat converetd to work specific energy [kJ/mol-MO]
q_s_sens_MO_array =Tab.q_s_sens_MO(indx_baseline_MVR,:).*1e-3;      % MO sensible heating required specific energy [kJ/mol-MO]
q_t_s_array = Tab.q_t_s(indx_baseline_MVR,:).*1e-3;                 % Total required specific energy [kJ/mol-MO]
bar(Tab.Redox_material(indx_baseline_MVR,:),[w_s_ox_sep_array q_s_ox_h_array w_s_inert_array q_s_inert_h_array q_s_sens_MO_array q_s_red_array -q_s_ox_exo_array]);
ylabel('Energy per produced H_2, kJ/mol_{H_2}');
legend({'w_{ox,sep}','q_{ox,h}','w_{sg,sep}','q_{sg,h}','q_{MO,h}','q_{red}','q_{ox}'},'NumColumns',2);
yscale('log');
if SAVEFLAG
    filename = 'Baseline_MVR_energy_absolute';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 6 - "baseline" conditions energy terms (relative)
fig = figure;   
F_s_ox_sep_array = Tab.F_ox_sep(indx_baseline_MVR,:);
F_s_ox_h_array = Tab.F_ox_h(indx_baseline_MVR,:);
F_s_inert_array = Tab.F_inert_sep(indx_baseline_MVR,:);
F_s_ox_exo_array = Tab.F_ox_exo(indx_baseline_MVR,:);
F_s_red_array = Tab.F_red(indx_baseline_MVR,:);
F_s_inert_h_array = Tab.F_inert_h(indx_baseline_MVR,:);
F_s_work_rec_array = Tab.F_work_rec(indx_baseline_MVR,:);
F_s_sens_MO_array =Tab.F_sens_MO(indx_baseline_MVR,:);
F_t_s_array = Tab.F_total(indx_baseline_MVR,:);
b = bar(Tab.Redox_material(indx_baseline_MVR,:),[F_s_ox_sep_array F_s_ox_h_array F_s_inert_array F_s_inert_h_array F_s_sens_MO_array F_s_red_array F_s_ox_exo_array]);
ylabel('Energy term ratio to fuel energy');
legend({'w_{ox,sep}','q_{ox,h}','w_{sg,sep}','q_{sg,h}','q_{MO,h}','q_{red}','q_{ox}'});
if SAVEFLAG
    filename = 'Baseline_MVR_energy_relative';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 7 - "baseline" conditions efficiency and conversion extent for
% cond and MVR
fig = figure;   
hold on
nil = zeros(length(eta_baseline_MVR_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[eta_baseline_cond_array eta_baseline_MVR_array nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[nil nil X_conv_baseline_cond_array X_conv_baseline_MVR_array]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_eta_X';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 8 - no solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;   
hold on
nil = zeros(length(eta_no_solid_HR_MVR_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR,:),[eta_no_solid_HR_cond_array eta_no_solid_HR_MVR_array nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR,:),[nil nil X_conv_no_solid_HR_cond_array X_conv_no_solid_HR_MVR_array]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_eta_X';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 9 - full solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;   
hold on
nil = zeros(length(eta_full_solid_HR_MVR_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR,:),[eta_full_solid_HR_cond_array eta_full_solid_HR_MVR_array nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR,:),[nil nil X_conv_full_solid_HR_cond_array X_conv_full_solid_HR_MVR_array]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_eta_X';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 10 - temperature swing for different materials at all case - cond
T_red_all_case_array = Tab.T_red(indx_all_cases,:);
T_ox_all_case_array = Tab.T_ox(indx_all_cases,:);
DeltaT_all_case_array = T_red_all_case_array-T_ox_all_case_array;
fig = figure;   
ax = gca;
bar(Tab.Redox_material(indx_baseline_cond,:),[Delta_T_baseline_cond_array Delta_T_no_solid_HR_cond_array Delta_T_full_solid_HR_cond_array]);
ylabel('Temperature swing \DeltaT, K');
ax.YLim = [0 max(ax.YLim)+25];
legend(['\epsilon=',num2str(unique(Tab.eps_HR(indx_baseline_cond)))],['\epsilon=',num2str(unique(Tab.eps_HR(indx_no_solid_HR_cond)))],['\epsilon=',num2str(unique(Tab.eps_HR(indx_full_solid_HR_cond)))]);
if SAVEFLAG
    filename = 'DeltaT_cond';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 11 - temperature swing for different materials at all case - MVR
fig = figure;
ax = gca;
bar(Tab.Redox_material(indx_baseline_MVR,:),[Delta_T_baseline_MVR_array Delta_T_no_solid_HR_MVR_array Delta_T_full_solid_HR_MVR_array]);
ylabel('Temperature swing \DeltaT, K');
ax.YLim = [0 max(ax.YLim)+25];
legend(['\epsilon=',num2str(unique(Tab.eps_HR(indx_baseline_MVR)))],['\epsilon=',num2str(unique(Tab.eps_HR(indx_no_solid_HR_MVR)))],['\epsilon=',num2str(unique(Tab.eps_HR(indx_full_solid_HR_MVR)))]);
if SAVEFLAG
    filename = 'DeltaT_MVR';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 12 - "baseline" conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_baseline_cond_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[Delta_T_baseline_cond_array Delta_T_baseline_MVR_array nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[nil nil omega_ox_baseline_cond_array omega_ox_baseline_MVR_array]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_DelT_omega_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 13 - no solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_no_solid_HR_cond_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR,:),[Delta_T_no_solid_HR_cond_array Delta_T_no_solid_HR_MVR_array nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR,:),[nil nil omega_ox_no_solid_HR_cond_array omega_ox_no_solid_HR_MVR_array]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_DelT_omega_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 14 - full solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_full_solid_HR_cond_array),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR,:),[Delta_T_full_solid_HR_cond_array Delta_T_full_solid_HR_MVR_array nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR,:),[nil nil omega_ox_full_solid_HR_cond_array omega_ox_full_solid_HR_MVR_array]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_DelT_omega_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
%% Intermediate pressure
% Get baseline values - condensing case
eta_baseline_cond_array_int_p_ox = Tab.eta(indx_baseline_cond_int_p_ox,:);                            % System efficiency
X_conv_baseline_cond_array_int_p_ox = Tab.X_conv(indx_baseline_cond_int_p_ox,:);                      % Conversion extent
T_red_baseline_cond_array_int_p_ox = Tab.T_red(indx_baseline_cond_int_p_ox,:);                        % Reduction temperature [K]
T_ox_baseline_cond_array_int_p_ox = Tab.T_ox(indx_baseline_cond_int_p_ox,:);                          % Oxidation temperature [K]
Delta_T_baseline_cond_array_int_p_ox = T_red_baseline_cond_array_int_p_ox-T_ox_baseline_cond_array_int_p_ox;   % Temperature swing [K]
omega_red_baseline_cond_array_int_p_ox = Tab.omega_red(indx_baseline_cond_int_p_ox,:);                % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_baseline_cond_array_int_p_ox = Tab.omega_ox(indx_baseline_cond_int_p_ox,:);                  % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_baseline_cond_array_int_p_ox = Tab.delta_red(indx_baseline_cond_int_p_ox,:);                % Reduction extent - end of reduction
delta_ox_baseline_cond_array_int_p_ox = Tab.delta_ox(indx_baseline_cond_int_p_ox,:);                  % Reduction extent - end of oxidation
delta_phi_baseline_cond_array_int_p_ox = Tab.delta_phi(indx_baseline_cond_int_p_ox,:);                % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get no solid HR values - condensing case
eta_no_solid_HR_cond_array_int_p_ox = Tab.eta(indx_no_solid_HR_cond_int_p_ox,:);                              % System efficiency
X_conv_no_solid_HR_cond_array_int_p_ox = Tab.X_conv(indx_no_solid_HR_cond_int_p_ox,:);                        % Conversion extent
T_red_no_solid_HR_cond_array_int_p_ox = Tab.T_red(indx_no_solid_HR_cond_int_p_ox,:);                          % Reduction temperature [K]
T_ox_no_solid_HR_cond_array_int_p_ox = Tab.T_ox(indx_no_solid_HR_cond_int_p_ox,:);                            % Oxidation temperature [K]
Delta_T_no_solid_HR_cond_array_int_p_ox = T_red_no_solid_HR_cond_array_int_p_ox-T_ox_no_solid_HR_cond_array_int_p_ox;  % Temperature swing [K]
omega_red_no_solid_HR_cond_array_int_p_ox = Tab.omega_red(indx_no_solid_HR_cond_int_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_no_solid_HR_cond_array_int_p_ox = Tab.omega_ox(indx_no_solid_HR_cond_int_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_no_solid_HR_cond_array_int_p_ox = Tab.delta_red(indx_no_solid_HR_cond_int_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_no_solid_HR_cond_array_int_p_ox = Tab.delta_ox(indx_no_solid_HR_cond_int_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_no_solid_HR_cond_array_int_p_ox = Tab.delta_phi(indx_no_solid_HR_cond_int_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get full solid HR values - condensing case
eta_full_solid_HR_cond_array_int_p_ox = Tab.eta(indx_full_solid_HR_cond_int_p_ox,:);                              % System efficiency
X_conv_full_solid_HR_cond_array_int_p_ox = Tab.X_conv(indx_full_solid_HR_cond_int_p_ox,:);                        % Conversion extent
T_red_full_solid_HR_cond_array_int_p_ox = Tab.T_red(indx_full_solid_HR_cond_int_p_ox,:);                          % Reduction temperature [K]
T_ox_full_solid_HR_cond_array_int_p_ox = Tab.T_ox(indx_full_solid_HR_cond_int_p_ox,:);                            % Oxidation temperature [K]
Delta_T_full_solid_HR_cond_array_int_p_ox = T_red_full_solid_HR_cond_array_int_p_ox-T_ox_full_solid_HR_cond_array_int_p_ox;% Temperature swing [K]
omega_red_full_solid_HR_cond_array_int_p_ox = Tab.omega_red(indx_full_solid_HR_cond_int_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_full_solid_HR_cond_array_int_p_ox = Tab.omega_ox(indx_full_solid_HR_cond_int_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_full_solid_HR_cond_array_int_p_ox = Tab.delta_red(indx_full_solid_HR_cond_int_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_full_solid_HR_cond_array_int_p_ox = Tab.delta_ox(indx_full_solid_HR_cond_int_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_full_solid_HR_cond_array_int_p_ox = Tab.delta_phi(indx_full_solid_HR_cond_int_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get baseline values - MVR case
eta_baseline_MVR_array_int_p_ox = Tab.eta(indx_baseline_MVR_int_p_ox,:);                            % System efficiency
X_conv_baseline_MVR_array_int_p_ox = Tab.X_conv(indx_baseline_MVR_int_p_ox,:);                      % Conversion extent
T_red_baseline_MVR_array_int_p_ox = Tab.T_red(indx_baseline_MVR_int_p_ox,:);                        % Reduction temperature [K]
T_ox_baseline_MVR_array_int_p_ox = Tab.T_ox(indx_baseline_MVR_int_p_ox,:);                          % Oxidation temperature [K]
Delta_T_baseline_MVR_array_int_p_ox = T_red_baseline_MVR_array_int_p_ox-T_ox_baseline_MVR_array_int_p_ox;   % Temperature swing [K]
omega_red_baseline_MVR_array_int_p_ox = Tab.omega_red(indx_baseline_MVR_int_p_ox,:);                % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_baseline_MVR_array_int_p_ox = Tab.omega_ox(indx_baseline_MVR_int_p_ox,:);                  % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_baseline_MVR_array_int_p_ox = Tab.delta_red(indx_baseline_MVR_int_p_ox,:);                % Reduction extent - end of reduction
delta_ox_baseline_MVR_array_int_p_ox = Tab.delta_ox(indx_baseline_MVR_int_p_ox,:);                  % Reduction extent - end of oxidation
delta_phi_baseline_MVR_array_int_p_ox = Tab.delta_phi(indx_baseline_MVR_int_p_ox,:);                % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get no solid HR values - MVR case
eta_no_solid_HR_MVR_array_int_p_ox = Tab.eta(indx_no_solid_HR_MVR_int_p_ox,:);                              % System efficiency
X_conv_no_solid_HR_MVR_array_int_p_ox = Tab.X_conv(indx_no_solid_HR_MVR_int_p_ox,:);                        % Conversion extent
T_red_no_solid_HR_MVR_array_int_p_ox = Tab.T_red(indx_no_solid_HR_MVR_int_p_ox,:);                          % Reduction temperature [K]
T_ox_no_solid_HR_MVR_array_int_p_ox = Tab.T_ox(indx_no_solid_HR_MVR_int_p_ox,:);                            % Oxidation temperature [K]
Delta_T_no_solid_HR_MVR_array_int_p_ox = T_red_no_solid_HR_MVR_array_int_p_ox-T_ox_no_solid_HR_MVR_array_int_p_ox;  % Temperature swing [K]
omega_red_no_solid_HR_MVR_array_int_p_ox = Tab.omega_red(indx_no_solid_HR_MVR_int_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_no_solid_HR_MVR_array_int_p_ox = Tab.omega_ox(indx_no_solid_HR_MVR_int_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_no_solid_HR_MVR_array_int_p_ox = Tab.delta_red(indx_no_solid_HR_MVR_int_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_no_solid_HR_MVR_array_int_p_ox = Tab.delta_ox(indx_no_solid_HR_MVR_int_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_no_solid_HR_MVR_array_int_p_ox = Tab.delta_phi(indx_no_solid_HR_MVR_int_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get full solid HR values - MVR case
eta_full_solid_HR_MVR_array_int_p_ox = Tab.eta(indx_full_solid_HR_MVR_int_p_ox,:);                              % System efficiency
X_conv_full_solid_HR_MVR_array_int_p_ox = Tab.X_conv(indx_full_solid_HR_MVR_int_p_ox,:);                        % Conversion extent
T_red_full_solid_HR_MVR_array_int_p_ox = Tab.T_red(indx_full_solid_HR_MVR_int_p_ox,:);                          % Reduction temperature [K]
T_ox_full_solid_HR_MVR_array_int_p_ox = Tab.T_ox(indx_full_solid_HR_MVR_int_p_ox,:);                            % Oxidation temperature [K]
Delta_T_full_solid_HR_MVR_array_int_p_ox = T_red_full_solid_HR_MVR_array_int_p_ox-T_ox_full_solid_HR_MVR_array_int_p_ox;% Temperature swing [K]
omega_red_full_solid_HR_MVR_array_int_p_ox = Tab.omega_red(indx_full_solid_HR_MVR_int_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_full_solid_HR_MVR_array_int_p_ox = Tab.omega_ox(indx_full_solid_HR_MVR_int_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_full_solid_HR_MVR_array_int_p_ox = Tab.delta_red(indx_full_solid_HR_MVR_int_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_full_solid_HR_MVR_array_int_p_ox = Tab.delta_ox(indx_full_solid_HR_MVR_int_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_full_solid_HR_MVR_array_int_p_ox = Tab.delta_phi(indx_full_solid_HR_MVR_int_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Plot 15 - "baseline" conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_baseline_MVR_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_int_p_ox,:),[eta_baseline_cond_array_int_p_ox eta_baseline_MVR_array_int_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_int_p_ox,:),[nil nil X_conv_baseline_cond_array_int_p_ox X_conv_baseline_MVR_array_int_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_eta_X_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 16 - no solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_no_solid_HR_MVR_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_int_p_ox,:),[eta_no_solid_HR_cond_array_int_p_ox eta_no_solid_HR_MVR_array_int_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_int_p_ox,:),[nil nil X_conv_no_solid_HR_cond_array_int_p_ox X_conv_no_solid_HR_MVR_array_int_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_eta_X_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 17 - full solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_full_solid_HR_MVR_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_int_p_ox,:),[eta_full_solid_HR_cond_array_int_p_ox eta_full_solid_HR_MVR_array_int_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_int_p_ox,:),[nil nil X_conv_full_solid_HR_cond_array_int_p_ox X_conv_full_solid_HR_MVR_array_int_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_eta_X_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 18 - "baseline" conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_baseline_cond_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_int_p_ox,:),[Delta_T_baseline_cond_array_int_p_ox Delta_T_baseline_MVR_array_int_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_int_p_ox,:),[nil nil omega_ox_baseline_cond_array_int_p_ox omega_ox_baseline_MVR_array_int_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_DelT_omega_ox_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 19 - no solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_no_solid_HR_cond_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_int_p_ox,:),[Delta_T_no_solid_HR_cond_array_int_p_ox Delta_T_no_solid_HR_MVR_array_int_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_int_p_ox,:),[nil nil omega_ox_no_solid_HR_cond_array_int_p_ox omega_ox_no_solid_HR_MVR_array_int_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_DelT_omega_ox_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 20 - full solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_full_solid_HR_cond_array_int_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_int_p_ox,:),[Delta_T_full_solid_HR_cond_array_int_p_ox Delta_T_full_solid_HR_MVR_array_int_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_int_p_ox,:),[nil nil omega_ox_full_solid_HR_cond_array_int_p_ox omega_ox_full_solid_HR_MVR_array_int_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_DelT_omega_ox_int_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
%% High pressure
% Get baseline values - condensing case
eta_baseline_cond_array_high_p_ox = Tab.eta(indx_baseline_cond_high_p_ox,:);                            % System efficiency
X_conv_baseline_cond_array_high_p_ox = Tab.X_conv(indx_baseline_cond_high_p_ox,:);                      % Conversion extent
T_red_baseline_cond_array_high_p_ox = Tab.T_red(indx_baseline_cond_high_p_ox,:);                        % Reduction temperature [K]
T_ox_baseline_cond_array_high_p_ox = Tab.T_ox(indx_baseline_cond_high_p_ox,:);                          % Oxidation temperature [K]
Delta_T_baseline_cond_array_high_p_ox = T_red_baseline_cond_array_high_p_ox-T_ox_baseline_cond_array_high_p_ox;   % Temperature swing [K]
omega_red_baseline_cond_array_high_p_ox = Tab.omega_red(indx_baseline_cond_high_p_ox,:);                % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_baseline_cond_array_high_p_ox = Tab.omega_ox(indx_baseline_cond_high_p_ox,:);                  % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_baseline_cond_array_high_p_ox = Tab.delta_red(indx_baseline_cond_high_p_ox,:);                % Reduction extent - end of reduction
delta_ox_baseline_cond_array_high_p_ox = Tab.delta_ox(indx_baseline_cond_high_p_ox,:);                  % Reduction extent - end of oxidation
delta_phi_baseline_cond_array_high_p_ox = Tab.delta_phi(indx_baseline_cond_high_p_ox,:);                % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get no solid HR values - condensing case
eta_no_solid_HR_cond_array_high_p_ox = Tab.eta(indx_no_solid_HR_cond_high_p_ox,:);                              % System efficiency
X_conv_no_solid_HR_cond_array_high_p_ox = Tab.X_conv(indx_no_solid_HR_cond_high_p_ox,:);                        % Conversion extent
T_red_no_solid_HR_cond_array_high_p_ox = Tab.T_red(indx_no_solid_HR_cond_high_p_ox,:);                          % Reduction temperature [K]
T_ox_no_solid_HR_cond_array_high_p_ox = Tab.T_ox(indx_no_solid_HR_cond_high_p_ox,:);                            % Oxidation temperature [K]
Delta_T_no_solid_HR_cond_array_high_p_ox = T_red_no_solid_HR_cond_array_high_p_ox-T_ox_no_solid_HR_cond_array_high_p_ox;  % Temperature swing [K]
omega_red_no_solid_HR_cond_array_high_p_ox = Tab.omega_red(indx_no_solid_HR_cond_high_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_no_solid_HR_cond_array_high_p_ox = Tab.omega_ox(indx_no_solid_HR_cond_high_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_no_solid_HR_cond_array_high_p_ox = Tab.delta_red(indx_no_solid_HR_cond_high_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_no_solid_HR_cond_array_high_p_ox = Tab.delta_ox(indx_no_solid_HR_cond_high_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_no_solid_HR_cond_array_high_p_ox = Tab.delta_phi(indx_no_solid_HR_cond_high_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get full solid HR values - condensing case
eta_full_solid_HR_cond_array_high_p_ox = Tab.eta(indx_full_solid_HR_cond_high_p_ox,:);                              % System efficiency
X_conv_full_solid_HR_cond_array_high_p_ox = Tab.X_conv(indx_full_solid_HR_cond_high_p_ox,:);                        % Conversion extent
T_red_full_solid_HR_cond_array_high_p_ox = Tab.T_red(indx_full_solid_HR_cond_high_p_ox,:);                          % Reduction temperature [K]
T_ox_full_solid_HR_cond_array_high_p_ox = Tab.T_ox(indx_full_solid_HR_cond_high_p_ox,:);                            % Oxidation temperature [K]
Delta_T_full_solid_HR_cond_array_high_p_ox = T_red_full_solid_HR_cond_array_high_p_ox-T_ox_full_solid_HR_cond_array_high_p_ox;% Temperature swing [K]
omega_red_full_solid_HR_cond_array_high_p_ox = Tab.omega_red(indx_full_solid_HR_cond_high_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_full_solid_HR_cond_array_high_p_ox = Tab.omega_ox(indx_full_solid_HR_cond_high_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_full_solid_HR_cond_array_high_p_ox = Tab.delta_red(indx_full_solid_HR_cond_high_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_full_solid_HR_cond_array_high_p_ox = Tab.delta_ox(indx_full_solid_HR_cond_high_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_full_solid_HR_cond_array_high_p_ox = Tab.delta_phi(indx_full_solid_HR_cond_high_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get baseline values - MVR case
eta_baseline_MVR_array_high_p_ox = Tab.eta(indx_baseline_MVR_high_p_ox,:);                            % System efficiency
X_conv_baseline_MVR_array_high_p_ox = Tab.X_conv(indx_baseline_MVR_high_p_ox,:);                      % Conversion extent
T_red_baseline_MVR_array_high_p_ox = Tab.T_red(indx_baseline_MVR_high_p_ox,:);                        % Reduction temperature [K]
T_ox_baseline_MVR_array_high_p_ox = Tab.T_ox(indx_baseline_MVR_high_p_ox,:);                          % Oxidation temperature [K]
Delta_T_baseline_MVR_array_high_p_ox = T_red_baseline_MVR_array_high_p_ox-T_ox_baseline_MVR_array_high_p_ox;   % Temperature swing [K]
omega_red_baseline_MVR_array_high_p_ox = Tab.omega_red(indx_baseline_MVR_high_p_ox,:);                % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_baseline_MVR_array_high_p_ox = Tab.omega_ox(indx_baseline_MVR_high_p_ox,:);                  % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_baseline_MVR_array_high_p_ox = Tab.delta_red(indx_baseline_MVR_high_p_ox,:);                % Reduction extent - end of reduction
delta_ox_baseline_MVR_array_high_p_ox = Tab.delta_ox(indx_baseline_MVR_high_p_ox,:);                  % Reduction extent - end of oxidation
delta_phi_baseline_MVR_array_high_p_ox = Tab.delta_phi(indx_baseline_MVR_high_p_ox,:);                % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get no solid HR values - MVR case
eta_no_solid_HR_MVR_array_high_p_ox = Tab.eta(indx_no_solid_HR_MVR_high_p_ox,:);                              % System efficiency
X_conv_no_solid_HR_MVR_array_high_p_ox = Tab.X_conv(indx_no_solid_HR_MVR_high_p_ox,:);                        % Conversion extent
T_red_no_solid_HR_MVR_array_high_p_ox = Tab.T_red(indx_no_solid_HR_MVR_high_p_ox,:);                          % Reduction temperature [K]
T_ox_no_solid_HR_MVR_array_high_p_ox = Tab.T_ox(indx_no_solid_HR_MVR_high_p_ox,:);                            % Oxidation temperature [K]
Delta_T_no_solid_HR_MVR_array_high_p_ox = T_red_no_solid_HR_MVR_array_high_p_ox-T_ox_no_solid_HR_MVR_array_high_p_ox;  % Temperature swing [K]
omega_red_no_solid_HR_MVR_array_high_p_ox = Tab.omega_red(indx_no_solid_HR_MVR_high_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_no_solid_HR_MVR_array_high_p_ox = Tab.omega_ox(indx_no_solid_HR_MVR_high_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_no_solid_HR_MVR_array_high_p_ox = Tab.delta_red(indx_no_solid_HR_MVR_high_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_no_solid_HR_MVR_array_high_p_ox = Tab.delta_ox(indx_no_solid_HR_MVR_high_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_no_solid_HR_MVR_array_high_p_ox = Tab.delta_phi(indx_no_solid_HR_MVR_high_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Get full solid HR values - MVR case
eta_full_solid_HR_MVR_array_high_p_ox = Tab.eta(indx_full_solid_HR_MVR_high_p_ox,:);                              % System efficiency
X_conv_full_solid_HR_MVR_array_high_p_ox = Tab.X_conv(indx_full_solid_HR_MVR_high_p_ox,:);                        % Conversion extent
T_red_full_solid_HR_MVR_array_high_p_ox = Tab.T_red(indx_full_solid_HR_MVR_high_p_ox,:);                          % Reduction temperature [K]
T_ox_full_solid_HR_MVR_array_high_p_ox = Tab.T_ox(indx_full_solid_HR_MVR_high_p_ox,:);                            % Oxidation temperature [K]
Delta_T_full_solid_HR_MVR_array_high_p_ox = T_red_full_solid_HR_MVR_array_high_p_ox-T_ox_full_solid_HR_MVR_array_high_p_ox;% Temperature swing [K]
omega_red_full_solid_HR_MVR_array_high_p_ox = Tab.omega_red(indx_full_solid_HR_MVR_high_p_ox,:);                  % Omega reduction [mol_N2*s^-1 / mol_MO*s^-1]
omega_ox_full_solid_HR_MVR_array_high_p_ox = Tab.omega_ox(indx_full_solid_HR_MVR_high_p_ox,:);                    % Omega oxidation [mol_ox*s^-1 / mol_MO*s^-1]
delta_red_full_solid_HR_MVR_array_high_p_ox = Tab.delta_red(indx_full_solid_HR_MVR_high_p_ox,:);                  % Reduction extent - end of reduction
delta_ox_full_solid_HR_MVR_array_high_p_ox = Tab.delta_ox(indx_full_solid_HR_MVR_high_p_ox,:);                    % Reduction extent - end of oxidation
delta_phi_full_solid_HR_MVR_array_high_p_ox = Tab.delta_phi(indx_full_solid_HR_MVR_high_p_ox,:);                  % Oxygen released/product generated over amount of MO [mol-O/mol-MO] or [mol-prod/mol-MO]
% Plot 21 - "baseline" conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_baseline_MVR_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_high_p_ox,:),[eta_baseline_cond_array_high_p_ox eta_baseline_MVR_array_high_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_high_p_ox,:),[nil nil X_conv_baseline_cond_array_high_p_ox X_conv_baseline_MVR_array_high_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_eta_X_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 22 - no solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_no_solid_HR_MVR_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_high_p_ox,:),[eta_no_solid_HR_cond_array_high_p_ox eta_no_solid_HR_MVR_array_high_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_high_p_ox,:),[nil nil X_conv_no_solid_HR_cond_array_high_p_ox X_conv_no_solid_HR_MVR_array_high_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_eta_X_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 23 - full solid HR conditions efficiency and production rate ratio for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(eta_full_solid_HR_MVR_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_high_p_ox,:),[eta_full_solid_HR_cond_array_high_p_ox eta_full_solid_HR_MVR_array_high_p_ox nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_high_p_ox,:),[nil nil X_conv_full_solid_HR_cond_array_high_p_ox X_conv_full_solid_HR_MVR_array_high_p_ox]);
ylabel('Oxidizer conversion X_{ox}');
legend('\eta (cond)','\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_eta_X_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 24 - "baseline" conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_baseline_cond_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_high_p_ox,:),[Delta_T_baseline_cond_array_high_p_ox Delta_T_baseline_MVR_array_high_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR_high_p_ox,:),[nil nil omega_ox_baseline_cond_array_high_p_ox omega_ox_baseline_MVR_array_high_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_DelT_omega_ox_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 25 - no solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_no_solid_HR_cond_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_high_p_ox,:),[Delta_T_no_solid_HR_cond_array_high_p_ox Delta_T_no_solid_HR_MVR_array_high_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_no_solid_HR_MVR_high_p_ox,:),[nil nil omega_ox_no_solid_HR_cond_array_high_p_ox omega_ox_no_solid_HR_MVR_array_high_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'No_solid_HR_DelT_omega_ox_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
% Plot 26 - full solid HR conditions temperature swing and omega _ox for
% cond and MVR
fig = figure;
hold on
nil = zeros(length(Delta_T_full_solid_HR_cond_array_high_p_ox),1);
yyaxis left
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_high_p_ox,:),[Delta_T_full_solid_HR_cond_array_high_p_ox Delta_T_full_solid_HR_MVR_array_high_p_ox nil nil]);
ylabel('Temperature swing \Delta_T, K');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_full_solid_HR_MVR_high_p_ox,:),[nil nil omega_ox_full_solid_HR_cond_array_high_p_ox omega_ox_full_solid_HR_MVR_array_high_p_ox]);
ylabel('Oxidizer to oxide flow rate ratio \omega_{ox}');
legend('\DeltaT (cond)','\DeltaT (MVR)','\omega_{ox} (cond)','\omega_{ox} (MVR)');
if SAVEFLAG
    filename = 'Full_solid_HR_DelT_omega_ox_high_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
%% Pressure sweeps
p_ox_array = Tab.p_ox(indx_baseline_cond);
p_ox_int_array = Tab.p_ox(indx_baseline_cond_int_p_ox);
p_ox_high_array = Tab.p_ox(indx_baseline_cond_high_p_ox);
% Plot 27 - system efficiency and conversion as a function of oxidation
% pressure - baseline case, condensing
fig = figure;
hold on
nil = zeros(length(eta_baseline_cond_array),1);
yyaxis left
ax = gca;
ax.YColor = 'k';
colororder('default')
bar(Tab.Redox_material(indx_baseline_MVR,:),[eta_baseline_cond_array eta_baseline_cond_array_int_p_ox nil nil]);
% bar(Tab.Redox_material(indx_baseline_MVR_high_p_ox,:),[eta_baseline_cond_array eta_baseline_cond_array_int_p_ox eta_baseline_cond_array_high_p_ox nil nil nil]);
ylabel('System efficiency \eta');
yyaxis right
ax = gca;
colororder('default')
ax.YColor = 'k';
bar(Tab.Redox_material(indx_baseline_MVR,:),[nil nil X_conv_baseline_cond_array X_conv_baseline_cond_array_int_p_ox ]);
% bar(Tab.Redox_materialindx_baseline_MVR_high_p_ox,:),[nil nil nil X_conv_baseline_cond_array X_conv_baseline_cond_array_int_p_ox X_conv_baseline_cond_array_high_p_ox ]);
ylabel('Oxidizer conversion X_{ox}');
legend(['\eta, (p_{ox}=',num2str(p_ox_array(1)),' Pa)'],'\eta (MVR)','X_{ox} (cond)','X_{ox} (MVR)');
if SAVEFLAG
    filename = 'Baseline_eta_X_cond_vs_p_ox';
    savefig(fig,filename);
    saveas(fig,filename,'epsc');
    saveas(fig,filename,'emf');
    print(fig,filename,'-r1000','-dpng');
end
%% Functions
function result = CP_PropsSI(varargin) %#ok<DEFNU>
    % Shorthand version of CoolProp for MATLAB
    result = py.CoolProp.CoolProp.PropsSI(varargin{:});
end