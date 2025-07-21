% ------------------ General info -----------------
% Created by: Alon Lidor (alon.lidor@nrel.gov)
% National Renewable Energy Laboratory (NREL)
% Date: June 3, 2024
% This script is used to plot and analyze results from the 0D
% counter-current flow (CF) redox thermodynamic model:
% Redox_Countercurrent_Thermo_Main.m
% ------------------------------------------------
clearvars;
close all;
import py.CoolProp.CoolProp.*
% Load Results File
% uiopen('load');
[file_str,path_str] = uigetfile;
load(fullfile(path_str,file_str));
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
        % Check if folder exists
        temp_str = regexp(file_str,'.mat','split');
        folder_str = temp_str{1};
        if exist(folder_str,"dir")==7
            folder_flag = 1;
        else
            folder_flag = 0;
        end
        count = 1;
        while folder_flag==1
            folder_str = strcat(folder_str,'_',num2str(count));
            if exist(folder_str,"dir")==7
                count = count+1;
            else
                folder_flag = 0;
            end
        end
        mkdir(folder_str);
        addpath(folder_str);
end
% Selection if full process efficiency is calculated or not (1-calculated)
extra_calc_flag = 1;
% Fixed input
eta_CO2_sep = 0.1;          % CO-CO2 separation efficiency
eta_pump = 0.9;             % Efficiency of increasing the pressure of the oxidizier (cold) - pump/compressor
E_inert = 15e3;             % Cryogenic N2 separation energy [J/mol-N2]
eps_HR = 0.5;               % Solid heat recovery effectiveness
eps_HR_ox = 0.8;            % Exothermic oxidation heat recovery effectiveness
eta_ox_htw = 0.4;           % Efficiency of converting heat to work from excess exothermic heat
eps_g = 0.8;                % Gas-gas heat recovery effectiveness
T0 = 298.15;                % Ambient temperature [K]
T_pump = 300;               % Vacuum pump temperature [K]
f_th_loss = 0.1;            % Thermal losses (reduction reactor)
% Select H2-H2O separation technology
prod_sep_flag = input('Select H2-H2O separation technology (1-condensation,2-mechanical vapor recompression) [1]: ');
if isempty(prod_sep_flag)
    prod_sep_flag = 1;
end
% MFR flag - Note: MFR not fully implemented yet - keep 0 if
% 'extra_calc_flag' is equal to 1
MFR_flag  = 0;
%% Display general data
disp(['Metal oxide: ',MO_label]);
disp(['Target product: ',prod_str]);
disp('---Reduction default values---');
disp(['Reduction temperature (T_red): ',num2str(T_red-273.15), ' deg. C']);
disp(['Reduction pressure (p_red): ',num2str(p_red), ' Pa']);
disp(['Molar flow rate ratio of sweep gas to metal oxide (omega_red): ',num2str(omega_red)]);
disp(['O2 mole fraction in sweep gas (phi_O2): ',num2str(phi)]);
disp('---Oxidation default values---');
disp(['Oxidation temperature (T_ox): ',num2str(T_ox-273.15), ' deg. C']);
disp(['Oxidation pressure (p_ox): ',num2str(p_ox), ' Pa']);
disp(['Molar flow rate ratio of oxidizer gas to metal oxide (omega_ox): ',num2str(omega_ox)]);
disp(['Product inlet partial pressure (',p_str,'): ',num2str(p_prod_in),' Pa']);
%% Plotting
% Plot reduction
fig = figure;
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
    savefig(fig,fullfile(folder_str,fig_name));
    saveas(fig,fullfile(folder_str,fig_name),'epsc');
    saveas(fig,fullfile(folder_str,fig_name),'emf');
    print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
end
% Plot oxidation
fig = figure;
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
    savefig(fig,fullfile(folder_str,fig_name));
    saveas(fig,fullfile(folder_str,fig_name),'epsc');
    saveas(fig,fullfile(folder_str,fig_name),'emf');
    print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
end
%% Plot parametric sweep results over 'omega_red' and 'omega_ox'
if para_input_omega=='Y'
    % Plot reduction
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    p1 = plot(omega_red_par,delta_red_PF_par_omega,'-k');
    hold on;
    p2 = plot(omega_red_par,delta_red_CF_par_omega,'--k');
    if MFR_flag==1
        p3 = plot(omega_red_par,delta_red_MFR_par_omega,'-.k');
    end
    % xlim([0 kappa_ox_max]);
    % ylim([1e-10 1e5]);
    % yscale('log');
    xscale(ax,'log');
    xlabel('\omega_{red}','FontSize',fontsize_l);
    ylabel('\delta_{max}','FontSize',fontsize_l);
    % line(ax,[kappa(ind_PF_ox) kappa(ind_PF_ox)],[1e-10 pO2_gas_PF_ox(ind_PF_ox)],'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',1);
    % line(ax,[kappa(ind_CF_ox) kappa(ind_CF_ox)],[1e-10 pO2_gas_CF_ox(ind_CF_ox)],'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',1);
    if MFR_flag==1
        legend([p1 p2 p3],{'PF','CF','MFR'},'Location','best');
    else
        legend([p1 p2],{'PF','CF'},'Location','best');
    end
    title(['T_{red}=',num2str(T_red),' K, \phi=',num2str(phi)]);
    if SAVEFLAG
        fig_name = 'Reduction_Par_Omega';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot oxidation
    fig = figure;
    yyaxis left
    ax = gca;
    ax.YColor = 'k';
    ax.FontSize = fontsize_a;
    p1 = plot(omega_ox_par,delta_delta_PF_par_omega(omega_red_par==omega_red,:),'-k');
    hold on;
    p2 = plot(omega_ox_par,delta_delta_CF_par_omega(omega_red_par==omega_red,:),'-r');
    if MFR_flag==1
        p3 = plot(omega_ox_par,delta_delta_MFR_par_omega(omega_red_par==omega_red,:),'-b');
    end
    % xlim([0 kappa_ox_max]);
    % ylim([1e-10 1e5]);
    % yscale('log');
    xscale(ax,'log');
    xlabel('\omega_{ox}');
    ylabel('\Delta\delta_{max}');
    line(ax,[omega_ox_min omega_ox_max],[delta_red_PF delta_red_PF],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    line(ax,[omega_ox_min omega_ox_max],[delta_red_CF delta_red_CF],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    if MFR_flag==1
        line(ax,[omega_ox_min omega_ox_max],[delta_red_MFR delta_red_MFR],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    end
    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    yscale(ax,'log');
    plot(omega_ox_par,X_PF_par_omega(omega_red_par==omega_red,:),'--k');
    hold on;
    plot(omega_ox_par,X_CF_par_omega(omega_red_par==omega_red,:),'--r');
    if MFR_flag==1
        plot(omega_ox_par,X_MFR_par_omega(omega_red_par==omega_red,:),'--b');
    end
    ylabel('X');
    if MFR_flag==1
        legend([p1 p2 p3],{'PF','CF','MFR'},'Location','best');
    else
        legend([p1 p2],{'PF','CF'},'Location','best');
    end
    % title(['T_{red}=',num2str(T_red),' K, \phi=',num2str(phi),', T_{ox}=',num2str(T_ox),' K, ',p_str,'=',num2str(pCO_in),' Pa']);
    title(['T_{red}=',num2str(T_red),' K, \phi=',num2str(phi),', T_{ox}=',num2str(T_ox),' K, ',ox_title_str]);
    ax.FontSize = fontsize_a;
    if SAVEFLAG
        fig_name = 'Oxidation_Par_Omega';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (PF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,delta_delta_PF_par_omega'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,delta_delta_PF_par_omega','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_Delta_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (CF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,delta_delta_CF_par_omega'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,delta_delta_CF_par_omega','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_Delta_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (MFR)
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,delta_delta_MFR_par_omega'); 
        hold on;
        shading interp;
        colormap(jet);
        cb = colorbar;
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(omega_red_par,omega_ox_par,delta_delta_MFR_par_omega','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_Delta_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    % Plot reduction and oxidation - conversion (PF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,X_PF_par_omega'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,X_PF_par_omega','LineColor','k','ShowText','on','LabelSpacing',150); 
    ylabel(cb,'X','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_X_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - conversion (CF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,X_CF_par_omega'); %#ok<*NASGU>
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [~,~] = contour(omega_red_par,omega_ox_par,X_CF_par_omega','LineColor','k','ShowText','on','LabelSpacing',150);
    ylabel(cb,'X','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_X_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - conversion (MFR)
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,X_MFR_par_omega');
        hold on;
        shading interp;
        colormap(jet);
        cb = colorbar;
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(omega_red_par,omega_ox_par,X_MFR_par_omega','LineColor','k','ShowText','on','LabelSpacing',150);
        ylabel(cb,'X','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        % text(min(omega_red_par)+0.01,max(omega_ox_par)/5,'MFR','Color','w');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_X_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
end
%% Plot parametric sweep results over 'T_red' and 'T_ox'
if para_input_T=='Y'
    % Plot reduction
    fig = figure;
    ax = gca;
    plot(T_red_par,delta_red_PF_par_T,'-k');
    hold on;
    plot(T_red_par,delta_red_CF_par_T,'-.k');
    if MFR_flag==1
        plot(T_red_par,delta_red_MFR_par_T,'--k');
    end
    % xlim([0 kappa_ox_max]);
    % ylim([1e-10 1e5]);
    % yscale('log');
    xscale(ax,'log');
    xlabel('T_{red}, K');
    ylabel('\delta_{max}');
    if MFR_flag==1
        legend({'PF','CF','MFR'},'Location','best');
    else
        legend({'PF','CF'},'Location','best');
    end
    title(['\omega_{red}=',num2str(omega_red),', \phi=',num2str(phi)]);
    ax.FontSize = fontsize_a;
    if SAVEFLAG
        fig_name = 'Reduction_Par_T';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot oxidation
    fig = figure;
    yyaxis left
    ax = gca;
    ax.YColor = 'k';
    p1 = plot(T_ox_par,delta_delta_PF_par_T(T_red_par==T_red,:),'-k');
    hold on;
    p2 = plot(T_ox_par,delta_delta_CF_par_T(T_red_par==T_red,:),'-r');
    if MFR_flag==1
        p3 = plot(T_ox_par,delta_delta_MFR_par_T(T_red_par==T_red,:),'-b');
    end
    l1 = line(ax,[T_ox_min T_ox_max],[delta_red_PF delta_red_PF],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    l2 = line(ax,[T_ox_min T_ox_max],[delta_red_CF delta_red_CF],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    if MFR_flag==1
        l3 = line(ax,[T_ox_min T_ox_max],[delta_red_MFR delta_red_MFR],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',1,'MarkerEdgeColor','none');
    end
    xlim([T_ox_min T_ox_max]);
    % ylim([1e-10 1e5]);
    % yscale('log');
    xlabel('T_{ox}, K');
    ylabel('\Delta\delta_{max}');
    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    yscale(ax,'log');
    plot(T_ox_par,X_PF_par_T(T_red_par==T_red,:),'--k');
    hold on;
    plot(T_ox_par,X_CF_par_T(T_red_par==T_red,:),'--r');
    if MFR_flag==1
        plot(T_ox_par,X_MFR_par_T(T_red_par==T_red,:),'--b');
    end
    ylabel('X');
    if MFR_flag==1
        legend([p1 p2 p3],{'PF','CF','MFR'},'Location','best');
    else
        legend([p1 p2],{'PF','CF'},'Location','best');
    end
    title(['T_{red}=',num2str(T_red),' K, ','\omega_{red}=',num2str(omega_red),', \phi=',num2str(phi),' \omega_{ox}=',num2str(omega_ox),', ',p_str,'=',num2str(pCO_in),' Pa']);
    ax.FontSize = fontsize_a;
    if SAVEFLAG
        fig_name = 'Oxidation_Par_T';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (PF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,delta_delta_PF_par_T'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,delta_delta_PF_par_T','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_Delta_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (CF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,delta_delta_CF_par_T'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,delta_delta_CF_par_T','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_Delta_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - Delta_delta (MFR)
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,delta_delta_MFR_par_T'); 
        hold on;
        shading interp;
        colormap(jet);
        cb = colorbar;
        xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(T_red_par,T_ox_par,delta_delta_MFR_par_T','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'\Delta\delta','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_Delta_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    % Plot reduction and oxidation - conversion (PF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,X_PF_par_T'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,X_PF_par_T','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'X','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_X_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - conversion (CF)
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,X_CF_par_T');
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,X_CF_par_T','LineColor','k','ShowText','on','LabelSpacing',150); 
    ylabel(cb,'X','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_X_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot reduction and oxidation - conversion (MFR)
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,X_MFR_par_T');
        hold on;
        shading interp;
        colormap(jet);
        cb = colorbar;
        xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        % levels = 0:0.001:max(delta_delta_PF_par_omega);
        % [Mn,cn] = contour(omega_red_par,omega_ox_par,delta_delta_PF_par_omega,levels,'LineColor','k','ShowText','on','LabelSpacing',150);
        [Mn,cn] = contour(T_red_par,T_ox_par,X_MFR_par_T','LineColor','k','ShowText','on','LabelSpacing',150);
        ylabel(cb,'X','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_X_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
end
%% Calculations of effects - efficiency
if extra_calc_flag == 1
    % Fuel specific energy
    % PF
    Q_fuel_PF_par_omega = HHV.*delta_phi_PF_par_omega;        % Specific fuel energy [J/mol-MO]
    Q_fuel_PF_par_T = HHV.*delta_phi_PF_par_T;                % Specific fuel energy [J/mol-MO]
    % CF
    Q_fuel_CF_par_omega = HHV.*delta_phi_CF_par_omega;        % Specific fuel energy [J/mol-MO]
    Q_fuel_CF_par_T = HHV.*delta_phi_CF_par_T;                % Specific fuel energy [J/mol-MO]
    % MFR
    if MFR_flag==1
        Q_fuel_MFR_par_omega = HHV.*delta_phi_MFR_par_omega;  % Specific fuel energy [J/mol-MO]
        Q_fuel_MFR_par_T = HHV.*delta_phi_MFR_par_T;          % Specific fuel energy [J/mol-MO]
    end
    % Oxidizier-product separation and oxidizer sensible heating
    switch K_input      % Select whether it's H2O or CO2 splitting
        case 1          % H2-H2O separation
            h_fg = py.CoolProp.CoolProp.PropsSI('HMOLAR','P',p_ox,'Q',1,'H2O')-py.CoolProp.CoolProp.PropsSI('HMOLAR','P',p_ox,'Q',0,'H2O');     % Enthalpy of vaporization at 1 atm [J/mol-H2O]
            Q_liq_heat = max(CP_PropsSI('HMOLAR','Q',0,'P',p_ox,'H2O')-CP_PropsSI('HMOLAR','T',T0,'P',p_ox,'H2O'),0);                           % Specific heat to heat liquid at T0 to Tsat [J/mol_H2O]
            h_H2O_ox = CP_PropsSI('HMOLAR','T',T_ox,'P',p_ox,'H2O');    % H2O specific enthalpy at T_ox [J/mol_H2O]
            h_H2O_satvap = CP_PropsSI('HMOLAR','Q',1,'P',p_ox,'H2O');   % H2O specific enthalpy at T_sat_H2O(p_ox) [J/mol_H2O] - saturated vapor
            h_H2_ox = CP_PropsSI('HMOLAR','T',T_ox,'P',p_ox,'H2');      % H2 specific enthalpy at T_ox [J/mol_H2]
            h_H2O_ox_par_T = zeros(length(T_red_par),length(T_ox_par));
            h_H2_ox_par_T = zeros(length(T_red_par),length(T_ox_par));
            for J=1:length(T_ox_par)
                h_H2O_ox_par_T(:,J) = CP_PropsSI('HMOLAR','T',T_ox_par(J),'P',p_ox,'H2O');    % H2O specific enthalpy at T_ox [J/mol_H2O]
                h_H2_ox_par_T(:,J) = CP_PropsSI('HMOLAR','T',T_ox_par(J),'P',p_ox,'H2');      % H2 specific enthalpy at T_ox [J/mol_H2]
            end
            h_ox_PF_par_omega = h_H2O_ox.*(1-X_PF_par_omega)+h_H2_ox.*X_PF_par_omega;       % Effluent specific enthalpy at T_ox [J/mol] - PF, omega_par
            h_ox_CF_par_omega = h_H2O_ox.*(1-X_CF_par_omega)+h_H2_ox.*X_CF_par_omega;       % Effluent specific enthalpy at T_ox [J/mol] - CF, omega_par
            h_ox_PF_par_T = h_H2O_ox_par_T.*(1-X_PF_par_T)+h_H2_ox_par_T.*X_PF_par_T;       % Effluent specific enthalpy at T_ox [J/mol] - PF, omega_T
            h_ox_CF_par_T = h_H2O_ox_par_T.*(1-X_CF_par_T)+h_H2_ox_par_T.*X_CF_par_T;       % Effluent specific enthalpy at T_ox [J/mol] - CF, omega_T
            switch prod_sep_flag
                case 1      % Condensing and reboiling
                    F_ox_h_max = 10;                        % Maximum ratio of oxidizer sensible heating energy to fuel energy (for plotting)
                    levels_ox_h = 0:F_ox_h_max;             % Contour levels for colormap plotting
                    Q_ox_sep_PF_par_omega = (h_fg+Q_liq_heat).*omega_ox_par;    % Specific H2O boiling heat and liquid heating to boiling point [J/mol-MO] - PF, omega_par
                    Q_ox_sep_CF_par_omega = (h_fg+Q_liq_heat).*omega_ox_par;    % Specific H2O boiling heat and liquid heating to boiling point [J/mol-MO] - CF, omega_par
                    Q_ox_sep_PF_par_T = (h_fg+Q_liq_heat).*omega_ox;            % Specific H2O boiling heat and liquid heating to boiling point [J/mol-MO] - PF, omega_T
                    Q_ox_sep_CF_par_T = (h_fg+Q_liq_heat).*omega_ox;            % Specific H2O boiling heat and liquid heating to boiling point [J/mol-MO] - CF, omega_T
                    W_ox_sep_PF_par_omega = 0;                                  % Work for effluent stream product separaion [J/mol-MO] - PF, omega_par
                    W_ox_sep_CF_par_omega = 0;                                  % Work for effluent stream product separaion [J/mol-MO] - CF, omega_par
                    W_ox_sep_PF_par_T = 0;                                      % Work for effluent stream product separaion [J/mol-MO] - PF, T_par
                    W_ox_sep_CF_par_T = 0;                                      % Work for effluent stream product separaion [J/mol-MO] - CF, T_par
                case 2      % Mechanical vapor recompression cycle separation
                    F_ox_h_max = 2;                             % Maximum ratio of oxidizer sensible heating energy to fuel energy (for plotting)
                    F_prod_sep_max = 2;                         % Maximum ratio of oxidizer sensible heating energy to fuel energy (for plotting)
                    levels_ox_h = 0:0.2:F_ox_h_max;             % Contour levels for colormap plotting
                    levels_prod_sep = 0:0.2:F_prod_sep_max;     % Contour levels for colormap plotting
                    T_water_in = 25+273.15;     % Fresh water inlet temperature [K]
                    eta_comp_MVR = 0.87;        % Vapor compressor efficiency
                    p_MVR = p_ox*2;             % Compression pressure in the MVR cycle [Pa]
                    T_MVR_in = CP_PropsSI('T','P',p_ox,'Q',1,'H2O')+15;      % Inlet temperature into the MVR cycle (outlet temperature from the gas-gas HX hot side for the oxidation loop) [K]
                    Q_ox_sep_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
                    Q_ox_sep_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
                    Q_ox_sep_PF_par_T = zeros(length(T_red_par),length(T_ox_par));
                    Q_ox_sep_CF_par_T = zeros(length(T_red_par),length(T_ox_par));
                    W_ox_sep_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
                    W_ox_sep_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));
                    W_ox_sep_PF_par_T = zeros(length(T_red_par),length(T_ox_par));
                    W_ox_sep_CF_par_T = zeros(length(T_red_par),length(T_ox_par));
                    for I=1:length(omega_red_par)
                        for J=1:length(omega_ox_par)
                            % PF, omega_par
                            [~,~,~,~,~,~,~,~,~,~,~,~,~,Q_dot_extra_PF_par_omega,W_dot_comp_PF_par_omega,~,~,~] = MVR_H2_H2O(omega_ox_par(J),T_MVR_in,p_ox,p_MVR,T_water_in,eta_comp_MVR,X_PF_par_omega(I,J));
                            W_ox_sep_PF_par_omega(I,J) = W_dot_comp_PF_par_omega;                  % MVR compressor work in [J/mol-MO]
                            Q_ox_sep_PF_par_omega(I,J) = max(Q_dot_extra_PF_par_omega,0);          % Specific H2O condensing heat [J/mol-MO]
                            % CF, omega_par
                            [~,~,~,~,~,~,~,~,~,~,~,~,~,Q_dot_extra_CF_par_omega,W_dot_comp_CF_par_omega,~,~,~] = MVR_H2_H2O(omega_ox_par(J),T_MVR_in,p_ox,p_MVR,T_water_in,eta_comp_MVR,X_CF_par_omega(I,J));
                            W_ox_sep_CF_par_omega(I,J) = W_dot_comp_CF_par_omega;                  % MVR compressor work in [J/mol-MO]
                            Q_ox_sep_CF_par_omega(I,J) = max(Q_dot_extra_CF_par_omega,0);          % Specific H2O condensing heat [J/mol-MO]
                        end
                    end
                    for I=1:length(T_red_par)
                        for J=1:length(T_ox_par)
                            % PF, T_par
                            [~,~,~,~,~,~,~,~,~,~,~,~,~,Q_dot_extra_PF_par_T,W_dot_comp_PF_par_T,~,~,~] = MVR_H2_H2O(omega_ox,T_MVR_in,p_ox,p_MVR,T_water_in,eta_comp_MVR,X_PF_par_T(I,J));
                            W_ox_sep_PF_par_T(I,J) = W_dot_comp_PF_par_T;                  % MVR compressor work in [J/mol-MO]
                            Q_ox_sep_PF_par_T(I,J) = max(Q_dot_extra_PF_par_T,0);          % Specific H2O condensing heat [J/mol-MO]
                            % CF, T_par
                            [~,~,~,~,~,~,~,~,~,~,~,~,~,Q_dot_extra_CF_par_T,W_dot_comp_CF_par_T,~,~,~] = MVR_H2_H2O(omega_ox,T_MVR_in,p_ox,p_MVR,T_water_in,eta_comp_MVR,X_CF_par_T(I,J));
                            W_ox_sep_CF_par_T(I,J) = W_dot_comp_CF_par_T;                  % MVR compressor work in [J/mol-MO]
                            Q_ox_sep_CF_par_T(I,J) = max(Q_dot_extra_CF_par_T,0);          % Specific H2O condensing heat [J/mol-MO]
                        end
                    end  
            end
            % Oxidizer heating in [J/mol-MO]
            Q_ox_h_PF_par_omega = omega_ox_par.*(h_ox_PF_par_omega-h_H2O_satvap-(h_ox_PF_par_omega-h_H2O_satvap).*eps_g)+Q_ox_sep_PF_par_omega;
            Q_ox_h_CF_par_omega = omega_ox_par.*(h_ox_CF_par_omega-h_H2O_satvap-(h_ox_CF_par_omega-h_H2O_satvap).*eps_g)+Q_ox_sep_CF_par_omega;
            Q_ox_h_PF_par_T = omega_ox.*(h_ox_PF_par_T-h_H2O_satvap-(h_ox_PF_par_T-h_H2O_satvap).*eps_g)+Q_ox_sep_PF_par_T;
            Q_ox_h_CF_par_T = omega_ox.*(h_ox_CF_par_T-h_H2O_satvap-(h_ox_CF_par_T-h_H2O_satvap).*eps_g)+Q_ox_sep_CF_par_T;
            F_prod_sep_str = 'F_{prod,sep}';
            % Pumping work
            h_ox_in = CP_PropsSI('HMOLAR','T',T0,'P',1e5,'H2O');                % Specific enthalpy at p0 [J/mol-H2O]
            s_ox_in = CP_PropsSI('SMOLAR','T',T0,'P',1e5,'H2O');                % Specific entropy at p0 [J/mol-K-CO2]
            h_ox_out_s = CP_PropsSI('HMOLAR','SMOLAR',s_ox_in,'P',T_ox,'H2O');  % Specific enthalpy at p_ox [J/mol-CO2] - isentropic
            h_ox_out_s_par_T = zeros(length(T_red_par),length(T_ox_par));             % Specific enthalpy at p_ox [J/mol-CO2] - isentropic
            for I=1:length(T_ox_par)
                h_ox_out_s_par_T = CP_PropsSI('HMOLAR','SMOLAR',s_ox_in,'P',T_ox_par(I),'H2O');  % Specific enthalpy at p_ox [J/mol-CO2] - isentropic
            end
        case 2      % CO-CO2 separation
            F_ox_h_max = 1;                             % Maximum ratio of CO-CO2 separation energy to fuel energy
            F_prod_sep_max = 1;                         % Maximum ratio of oxidizer sensible heating energy to fuel energy (for plotting)
            levels_ox_h = 0:0.2:F_ox_h_max;             % Contour levels for colormap plotting
            levels_prod_sep = 0:0.2:F_prod_sep_max;     % Contour levels for colormap plotting
            Q_ox_h_PF_par_omega = -omega_ox_par.*R.*T0.*(X_PF_par_omega.*log(X_PF_par_omega)+(1-X_PF_par_omega).*log(1-X_PF_par_omega))/eta_CO2_sep;        % Specific CO2 separation power [J/mol-MO]
            Q_ox_h_PF_par_T = -omega_ox.*R.*T0.*(X_PF_par_T.*log(X_PF_par_T)+(1-X_PF_par_T).*log(1-X_PF_par_T))/eta_CO2_sep;                    % Specific CO2 separation power [J/mol-MO]
            Q_ox_h_CF_par_omega = -omega_ox_par.*R.*T0.*(X_CF_par_omega.*log(X_CF_par_omega)+(1-X_CF_par_omega).*log(1-X_CF_par_omega))/eta_CO2_sep;        % Specific CO2 separation power [J/mol-MO]
            Q_ox_h_CF_par_T = -omega_ox.*R.*T0.*(X_CF_par_T.*log(X_CF_par_T)+(1-X_CF_par_T).*log(1-X_CF_par_T))/eta_CO2_sep;                    % Specific CO2 separation power [J/mol-MO]
            if MFR_flag==1
                Q_ox_h_MFR_par_omega = -omega_ox_par.*R.*T0.*(X_CF_par_omega.*log(X_CF_par_omega)+(1-X_CF_par_omega).*log(1-X_CF_par_omega))/eta_CO2_sep;       % Specific CO2 separation power [J/mol-MO]
                Q_ox_h_MFR_par_T = -omega_ox.*R.*T0.*(X_MFR_par_T.*log(X_MFR_par_T)+(1-X_MFR_par_T).*log(1-X_MFR_par_T))/eta_CO2_sep;               % Specific CO2 separation power [J/mol-MO]
            end
            F_prod_sep_str = 'F_{CO-CO_2,sep}';
    end
    if p_ox>1e5
        W_ox_pump_par_omega = omega_ox_par.*(h_ox_out_s-h_ox_in)./eta_pump;
        W_ox_pump_par_T = omega_ox.*(h_ox_out_s_par_T-h_ox_in)./eta_pump;
    else
        W_ox_pump_par_omega = 0;
        W_ox_pump_par_T = 0;
    end
    % Ratio of oxidizier pumping work to fuel energy
    F_ox_pump_PF_par_omega = W_ox_pump_par_omega./Q_fuel_PF_par_omega;
    F_ox_pump_CF_par_omega = W_ox_pump_par_omega./Q_fuel_CF_par_omega;
    F_ox_pump_PF_par_T = W_ox_pump_par_T./Q_fuel_PF_par_T;
    F_ox_pump_CF_par_T = W_ox_pump_par_T./Q_fuel_CF_par_T;
    if MFR_flag==1
        F_ox_pump_MFR_par_omega = W_ox_pump_par_omega./Q_fuel_MFR_par_omega;
        F_ox_pump_MFR_par_T = W_ox_pump_par_T./Q_fuel_MFR_par_T;
    end
    % Vacuum pumping work - IN WORK
    % if p_red<1e5
    %     eta_p_vac = Pump_Efficiency(p_red,1e5);
    %     W_vac_PF_par_omega = (R.*T_pump/eta_p_vac).*(omega_red_par+nO2_max_PF).*log(1e5./p_red);
    % else
    %     W_vac = 0;
    % end
    % F_vac = W_vac/Q_fuel;
    %% Plot F_prod_sep for all cases
    % Ratio of product separation energy to fuel energy
    F_prod_sep_PF_par_omega = W_ox_sep_PF_par_omega./Q_fuel_PF_par_omega;
    F_prod_sep_CF_par_omega = W_ox_sep_CF_par_omega./Q_fuel_CF_par_omega;
    F_prod_sep_PF_par_T = W_ox_sep_PF_par_T./Q_fuel_PF_par_T;
    F_prod_sep_CF_par_T = W_ox_sep_CF_par_T./Q_fuel_CF_par_T;
    if (prod_sep_flag==2)||(K_input==2)
        % --- NO MFR OPTION HERE YET ---
        % Plot F_prod_sep_PF_par_omega
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,F_prod_sep_PF_par_omega'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_prod_sep_max]);
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(omega_red_par,omega_ox_par,F_prod_sep_PF_par_omega',levels_prod_sep,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,F_prod_sep_str,'FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_F_prod_sep_PF';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
        % Plot F_prod_sep_CF_par_omega
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,F_prod_sep_CF_par_omega'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_prod_sep_max]);
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(omega_red_par,omega_ox_par,F_prod_sep_CF_par_omega',levels_prod_sep,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,F_prod_sep_str,'FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_F_prod_sep_CF';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
        % Plot F_prod_sep_PF_par_T
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,F_prod_sep_PF_par_T'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_prod_sep_max]);
        xlabel(ax,'T_{red}','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(T_red_par,T_ox_par,F_prod_sep_PF_par_T',levels_prod_sep,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,F_prod_sep_str,'FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_F_prod_sep_PF';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
        % Plot F_prod_sep_CF_par_T
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,F_prod_sep_CF_par_T'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_prod_sep_max]);
        xlabel(ax,'T_{red}','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        [Mn,cn] = contour(T_red_par,T_ox_par,F_prod_sep_CF_par_T',levels_prod_sep,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,F_prod_sep_str,'FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_F_prod_sep_CF';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    %% Plot F_ox_h for all cases
    % Ratio of oxidizer sensible heating energy to fuel energy
    F_ox_h_PF_par_omega = Q_ox_h_PF_par_omega./Q_fuel_PF_par_omega;
    F_ox_h_PF_par_T = Q_ox_h_PF_par_T./Q_fuel_PF_par_T;
    F_ox_h_CF_par_omega = Q_ox_h_CF_par_omega./Q_fuel_CF_par_omega;
    F_ox_h_CF_par_T = Q_ox_h_CF_par_T./Q_fuel_CF_par_T;
    if MFR_flag==1
        F_ox_h_MFR_par_omega = Q_ox_h_MFR_par_omega./Q_fuel_MFR_par_omega;
        F_ox_h_MFR_par_T = Q_ox_h_MFR_par_T./Q_fuel_MFR_par_T;
    end
    % if ~any([F_cond_PF_par_omega F_cond_CF_par_omega])<F_cond_max
    %     F_cond_max = ceil(max([F_cond_PF_par_omega F_cond_CF_par_omega],[], ...
    %         "all"));
    % end
    % Plot F_ox_h_PF_par_omega
    % F_cond_PF_par_omega(F_cond_PF_par_omega>7) = NaN;
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,F_ox_h_PF_par_omega'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_ox_h_max]);
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,F_ox_h_PF_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_F_ox_h_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_ox_h_CF_par_omega
    % F_cond_CF_par_omega(F_cond_CF_par_omega>F_cond_max) = F_cond_max;
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,F_ox_h_CF_par_omega');
    hold on;
    shading interp;
    cm = colormap(jet);
    cb = colorbar;
    clim(ax,[0 F_ox_h_max]);
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,F_ox_h_CF_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_F_ox_h_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_ox_h_MFR
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,F_ox_h_MFR_par_omega');
        hold on;
        shading interp;
        cm = colormap(jet);
        cb = colorbar;
        clim(ax,[0 F_ox_h_max]);
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        levels_ox_h = 0:F_ox_h_max;
        [Mn,cn] = contour(omega_red_par,omega_ox_par,F_ox_h_MFR_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_F_ox_h_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    % Plot F_ox_h_PF_par_T
    % F_cond_PF_par_T(F_cond_PF_par_T>7) = NaN;
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,F_ox_h_PF_par_T'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_ox_h_max]);
    xlabel(ax,'T_{red}','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,F_ox_h_PF_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_F_ox_h_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_ox_h_CF_par_T
    % F_cond_CF_par_T(F_cond_CF_par_T>F_cond_max) = F_cond_max;
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,F_ox_h_CF_par_T');
    hold on;
    shading interp;
    cm = colormap(jet);
    cb = colorbar;
    clim(ax,[0 F_ox_h_max]);
    xlabel(ax,'T_{red}','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,F_ox_h_CF_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_F_ox_h_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_ox_h_MFR
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,F_ox_h_MFR_par_T');
        hold on;
        shading interp;
        cm = colormap(jet);
        cb = colorbar;
        clim(ax,[0 F_ox_h_max]);
        xlabel(ax,'T_{red}','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        levels_ox_h = 0:F_ox_h_max;
        [Mn,cn] = contour(T_red_par,T_ox_par,F_ox_h_MFR_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'F_{ox,h}','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_F_ox_h_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    %% Calculation of sweep gas separation energy
    F_inert_max = 5;     % Maximum ratio of inert gas separation energy to fuel energy
    % Cryogenic separation
    E_inert_min = 0.15*3600/py.CoolProp.CoolProp.PropsSI('DMOLAR','T',273.15,'P',1e5,'N2');     % Minimum N2 separation energy [J/mol]
    E_inert_max = 0.25*3600/py.CoolProp.CoolProp.PropsSI('DMOLAR','T',273.15,'P',1e5,'N2');     % Maximum N2 separation energy [J/mol]
    E_inert = 15e3;     % Selected N2 separation energy [J/mol-N2]
    % Pressure Swing Absorption
    % Input
    xO2_in_PF_para_T_red = (phi*omega_red+nO2_red_PF_par_T./(nO2_red_PF_par_T+omega_red));                          % Inlet O2 mole fraction (PF)
    xO2_in_CF_para_T_red = (phi*omega_red+nO2_red_CF_par_T./(nO2_red_CF_par_T+omega_red));                          % Inlet O2 mole fraction (CF)
    xO2_in_PF_para_omega_red = (phi*omega_red_par'+nO2_red_PF_par_omega./(nO2_red_PF_par_omega+omega_red_par'));    % Inlet O2 mole fraction (PF)
    xO2_in_CF_para_omega_red = (phi*omega_red_par'+nO2_red_CF_par_omega./(nO2_red_CF_par_omega+omega_red_par'));    % Inlet O2 mole fraction (CF)
    xO2_out = phi;              % Outlet O2 mole fraction
    eta_comp = 0.85;            % Compressor efficiency
    eta_PSA = 0.9655*xO2_out^0.1853;    % PSA efficiency
    p_in = p_red;               % Inlet PSA pressure [Pa]
    p_out = 7e5;                % Outlet PSA pressure [Pa]
    T_in = 300;                 % Inlet temperature [K]
    x = zeros(nsp,1);
    gas = GRI30;
    iN2 = speciesIndex(gas,'N2');
    % Calculation for T_red_par
    w_PSA_PF_par_T_red = zeros(length(T_red_par),1);
    w_PSA_CF_par_T_red = zeros(length(T_red_par),1);
    for I=1:length(T_red_par)
        % Solve for PF
        x(iN2) = 1-xO2_in_PF_para_T_red(I);    % Set N2 mole fraction
        x(iO2) = xO2_in_PF_para_T_red(I);      % Set O2 mole fraction
        % Calculate inlet state properties
        set(gas,'T',T_in,'P',p_in,'X',x);
        h_in = enthalpy_mass(gas);          % Inlet enthalpy [J/kg]
        s_in = entropy_mass(gas);           % Inlet entropy [J/kg-K]
        % Calculate outlet state properties - isentropic
        set(gas,'P',p_out,'S',s_in);
        h_out_s = enthalpy_mass(gas);       % Outlet entropy (isentropic) [J/kg-K]
        w_comp_s = h_out_s-h_in;            % Isentropic compressor work [J/kg]
        w_comp = w_comp_s*eta_comp;         % Actual compressor work [J/kg]
        h_out = w_comp+h_in;                % Outlet enthalpy [J/kg]
        % Calculate outlet state properties - actual
        set(gas,'P',p_out,'H',h_out);
        % T_out = temperature(gas);               % Outlet temperature [K]
        M_mix = meanMolecularWeight(gas)*1e-3;  % Mean molecular weight [kg/mol]
        % Calculate work
        w_PSA_PF_par_T_red(I) = w_comp*M_mix/((1-xO2_in_PF_para_T_red(I))*eta_PSA);      % PSA work [J/mol-N2]
        % Solve for CF
        x(iN2) = 1-xO2_in_CF_para_T_red(I);    % Set N2 mole fraction
        x(iO2) = xO2_in_CF_para_T_red(I);      % Set O2 mole fraction
        % Calculate inlet state properties
        set(gas,'T',T_in,'P',p_in,'X',x);
        h_in = enthalpy_mass(gas);          % Inlet enthalpy [J/kg]
        s_in = entropy_mass(gas);           % Inlet entropy [J/kg-K]
        % Calculate outlet state properties - isentropic
        set(gas,'P',p_out,'S',s_in);
        h_out_s = enthalpy_mass(gas);       % Outlet entropy (isentropic) [J/kg-K]
        w_comp_s = h_out_s-h_in;            % Isentropic compressor work [J/kg]
        w_comp = w_comp_s*eta_comp;         % Actual compressor work [J/kg]
        h_out = w_comp+h_in;                % Outlet enthalpy [J/kg]
        % Calculate outlet state properties - actual
        set(gas,'P',p_out,'H',h_out);
        % T_out = temperature(gas);               % Outlet temperature [K]
        M_mix = meanMolecularWeight(gas)*1e-3;  % Mean molecular weight [kg/mol]
        % Calculate work
        w_PSA_CF_par_T_red(I) = w_comp*M_mix/((1-xO2_in_CF_para_T_red(I))*eta_PSA);      % PSA work [J/mol-N2]
    end
    % Calculation for omega_red_par
    w_PSA_PF_par_omega_red = zeros(length(omega_red_par),1);
    w_PSA_CF_par_omega_red = zeros(length(omega_red_par),1);
    for I=1:length(omega_red_par)
        % Solve for PF
        x(iN2) = 1-xO2_in_PF_para_omega_red(I);    % Set N2 mole fraction
        x(iO2) = xO2_in_PF_para_omega_red(I);      % Set O2 mole fraction
        % Calculate inlet state properties
        set(gas,'T',T_in,'P',p_in,'X',x);
        h_in = enthalpy_mass(gas);          % Inlet enthalpy [J/kg]
        s_in = entropy_mass(gas);           % Inlet entropy [J/kg-K]
        % Calculate outlet state properties - isentropic
        set(gas,'P',p_out,'S',s_in);
        h_out_s = enthalpy_mass(gas);       % Outlet entropy (isentropic) [J/kg-K]
        w_comp_s = h_out_s-h_in;            % Isentropic compressor work [J/kg]
        w_comp = w_comp_s*eta_comp;         % Actual compressor work [J/kg]
        h_out = w_comp+h_in;                % Outlet enthalpy [J/kg]
        % Calculate outlet state properties - actual
        set(gas,'P',p_out,'H',h_out);
        % T_out = temperature(gas);               % Outlet temperature [K]
        M_mix = meanMolecularWeight(gas)*1e-3;  % Mean molecular weight [kg/mol]
        % Calculate work
        w_PSA_PF_par_omega_red(I) = w_comp*M_mix/((1-xO2_in_PF_para_omega_red(I))*eta_PSA);      % PSA work [J/mol-N2]
        % Solve for CF
        x(iN2) = 1-xO2_in_CF_para_omega_red(I);    % Set N2 mole fraction
        x(iO2) = xO2_in_CF_para_omega_red(I);      % Set O2 mole fraction
        % Calculate inlet state properties
        set(gas,'T',T_in,'P',p_in,'X',x);
        h_in = enthalpy_mass(gas);          % Inlet enthalpy [J/kg]
        s_in = entropy_mass(gas);           % Inlet entropy [J/kg-K]
        % Calculate outlet state properties - isentropic
        set(gas,'P',p_out,'S',s_in);
        h_out_s = enthalpy_mass(gas);       % Outlet entropy (isentropic) [J/kg-K]
        w_comp_s = h_out_s-h_in;            % Isentropic compressor work [J/kg]
        w_comp = w_comp_s*eta_comp;         % Actual compressor work [J/kg]
        h_out = w_comp+h_in;                % Outlet enthalpy [J/kg]
        % Calculate outlet state properties - actual
        set(gas,'P',p_out,'H',h_out);
        % T_out = temperature(gas);               % Outlet temperature [K]
        M_mix = meanMolecularWeight(gas)*1e-3;  % Mean molecular weight [kg/mol]
        % Calculate work
        w_PSA_CF_par_omega_red(I) = w_comp*M_mix/((1-xO2_in_CF_para_omega_red(I))*eta_PSA);      % PSA work [J/mol-N2]
    end
    % Calculate energy - select the minimum energy between PSA and
    % cryogenic separation
    Q_inert_PF_par_T = min(E_inert,w_PSA_PF_par_T_red)*omega_red;               % Specific required separation energy [J/mol-MO] (PF)
    Q_inert_PF_par_omega = min(E_inert,w_PSA_PF_par_omega_red).*omega_red_par'; % Specific required separation energy [J/mol-MO] (PF)
    F_inert_PF_par_T = Q_inert_PF_par_T./Q_fuel_PF_par_T;                       % Ratio of sweep gas separation energy to fuel energy (PF)
    F_inert_PF_par_omega = Q_inert_PF_par_omega./Q_fuel_PF_par_omega;           % Ratio of sweep gas separation energy to fuel energy (PF)
    Q_inert_CF_par_T = min(E_inert,w_PSA_CF_par_T_red)*omega_red;               % Specific required separation energy [J/mol-MO] (CF)
    Q_inert_CF_par_omega = min(E_inert,w_PSA_CF_par_omega_red).*omega_red_par'; % Specific required separation energy [J/mol-MO] (CF)
    F_inert_CF_par_T = Q_inert_CF_par_T./Q_fuel_CF_par_T;                       % Ratio of sweep gas separation energy to fuel energy (CF)
    F_inert_CF_par_omega = Q_inert_CF_par_omega./Q_fuel_CF_par_omega;           % Ratio of sweep gas separation energy to fuel energy (CF)
    if MFR_flag==1
        Q_inert_MFR_par_T = min(E_inert,w_PSA_CF_par_T_red)*omega_red;              % Specific required separation energy [J/mol-MO] (MFR)
        Q_inert_MFR_par_omega = min(E_inert,w_PSA_CF_par_omega_red).*omega_red_par';% Specific required separation energy [J/mol-MO] (MFR)
        F_inert_MFR_par_T = Q_inert_MFR_par_T./Q_fuel_MFR_par_T;                    % Ratio of sweep gas separation energy to fuel energy (MFR)
        F_inert_MFR_par_omega = Q_inert_MFR_par_omega./Q_fuel_MFR_par_omega;        % Ratio of sweep gas separation energy to fuel energy (MFR)
    end
    % Plot F_inert_PF_par_omega
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,F_inert_PF_par_omega'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_inert_max]);
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    levels_ox_h = 0:F_inert_max;
    [Mn,cn] = contour(omega_red_par,omega_ox_par,F_inert_PF_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{inert}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_F_inert_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_inert_CF_par_omega
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,F_inert_CF_par_omega'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_inert_max]);
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    levels_ox_h = 0:F_inert_max;
    [Mn,cn] = contour(omega_red_par,omega_ox_par,F_inert_CF_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{inert}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_F_inert_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_inert_MFR_par_omega
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(omega_red_par,omega_ox_par,F_inert_MFR_par_omega'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_inert_max]);
        xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
        ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        levels_ox_h = 0:F_inert_max;
        [Mn,cn] = contour(omega_red_par,omega_ox_par,F_inert_MFR_par_omega',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'F_{inert}','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_Omega_Map_F_inert_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    % Plot F_inert_PF_par_T
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,F_inert_PF_par_T'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_inert_max]);
    xlabel(ax,'T_{red}','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    levels_ox_h = 0:F_inert_max;
    [Mn,cn] = contour(T_red_par,T_ox_par,F_inert_PF_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{inert}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_F_inert_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_inert_CF_par_T
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,F_inert_CF_par_T'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    clim(ax,[0 F_inert_max]);
    xlabel(ax,'T_{red}','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    levels_ox_h = 0:F_inert_max;
    [Mn,cn] = contour(T_red_par,T_ox_par,F_inert_CF_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'F_{inert}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_F_inert_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot F_inert_MFR_par_T
    if MFR_flag==1
        fig = figure;
        ax = gca;
        ax.FontSize = fontsize_a;
        s = pcolor(T_red_par,T_ox_par,F_inert_MFR_par_T'); 
        hold on;
        shading interp;
        cm = colormap(jet); 
        cb = colorbar;
        clim(ax,[0 F_inert_max]);
        xlabel(ax,'T_{red}','FontSize',fontsize_l);
        ylabel(ax,'T_{ox}','FontSize',fontsize_l);
        ax.XScale = 'log';
        ax.YScale = 'log';
        levels_ox_h = 0:F_inert_max;
        [Mn,cn] = contour(T_red_par,T_ox_par,F_inert_MFR_par_T',levels_ox_h,'LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
        ylabel(cb,'F_{inert}','FontSize',fontsize_l);
        ax.XAxis.Color = 'k';
        ax.XAxis.LineWidth = 1.5;
        ax.YAxis.Color = 'k';
        ax.YAxis.LineWidth = 1.5;
        dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
        legend(dummy,'MFR','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
        if SAVEFLAG
            fig_name = 'Par_T_Map_F_inert_MFR';
            savefig(fig,fullfile(folder_str,fig_name));
            saveas(fig,fullfile(folder_str,fig_name),'epsc');
            saveas(fig,fullfile(folder_str,fig_name),'emf');
            print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
        end
    end
    %% Calculation of sensible heating of MO
    Q_sens_par_omega = integral(cp_s_fun,T_ox,T_red,'ArrayValued',true).*M_MO;  %#ok<*FUNFUN> % Sensible heat in [J/mol-MO] - parametric omega study
    Q_sens_par_T = zeros(length(T_red_par),length(T_ox_par));                   % Sensible heat in [J/mol-MO] - parametric T study
    for I=1:length(T_red_par)
        for J=1:length(T_ox_par)
            Q_sens_par_T(I,J) = integral(cp_s_fun,T_ox_par(J),T_red_par(I),'ArrayValued',true).*M_MO;
        end
    end
    F_sens_PF_par_omega = Q_sens_par_omega./Q_fuel_PF_par_omega;    % Ratio of sensible heat to fuel energy (PF)
    F_sens_PF_par_T = Q_sens_par_T./Q_fuel_PF_par_T;                % Ratio of sensible heat to fuel energy (PF)
    F_sens_CF_par_omega = Q_sens_par_omega./Q_fuel_CF_par_omega;    % Ratio of sensible heat to fuel energy (CF)
    F_sens_CF_par_T = Q_sens_par_T./Q_fuel_CF_par_T;                % Ratio of sensible heat to fuel energy (CF)
    if MFR_flag==1
        F_sens_MFR_par_omega = Q_sens_par_omega./Q_fuel_MFR_par_omega;  % Ratio of sensible heat to fuel energy (MFR)
        F_sens_MFR_par_T = Q_sens_par_T./Q_fuel_MFR_par_T;              % Ratio of sensible heat to fuel energy (MFR)
    end
    %% Calculation of chemical reduction energy
    Q_red_PF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));     % Specific reduction energy [J/mol-MO] (PF)
    Q_red_CF_par_omega = zeros(length(omega_red_par),length(omega_ox_par));     % Specific reduction energy [J/mol-MO] (CF)
    if MFR_flag==1
        Q_red_MFR_par_omega = zeros(length(omega_red_par),length(omega_ox_par));    % Specific reduction energy [J/mol-MO] (MFR)
    end
    for I=1:length(omega_red_par)
        for J=1:length(omega_ox_par)
            Q_red_PF_par_omega(I,J) = integral(dH_fun,delta_ox_PF_par_omega(I,J),delta_red_PF_par_omega(I));
            Q_red_CF_par_omega(I,J) = integral(dH_fun,delta_ox_CF_par_omega(I,J),delta_red_CF_par_omega(I));
            if MFR_flag==1
                Q_red_MFR_par_omega(I,J) = integral(dH_fun,delta_ox_MFR_par_omega(I,J),delta_red_MFR_par_omega(I));
            end
        end
    end
    Q_red_PF_par_T = zeros(length(T_red_par),length(T_ox_par));     % Specific reduction energy [J/mol-MO] (PF)
    Q_red_CF_par_T = zeros(length(T_red_par),length(T_ox_par));     % Specific reduction energy [J/mol-MO] (CF)
    if MFR_flag==1
        Q_red_MFR_par_T = zeros(length(T_red_par),length(T_ox_par));    % Specific reduction energy [J/mol-MO] (MFR)
    end
    for I=1:length(T_red_par)
        for J=1:length(T_ox_par)
            Q_red_PF_par_T(I,J) = integral(dH_fun,delta_ox_PF_par_T(I,J),delta_red_PF_par_T(I));
            Q_red_CF_par_T(I,J) = integral(dH_fun,delta_ox_CF_par_T(I,J),delta_red_CF_par_T(I));
            if MFR_flag==1
                Q_red_MFR_par_T(I,J) = integral(dH_fun,delta_ox_MFR_par_T(I,J),delta_red_MFR_par_T(I));
            end
        end
    end
    F_red_PF_par_omega = Q_red_PF_par_omega./Q_fuel_PF_par_omega;        % Reduction specific energy [J/mol-MO] (PF)
    F_red_PF_par_T = Q_red_PF_par_T./Q_fuel_PF_par_T;                    % Reduction specific energy [J/mol-MO] (PF)
    F_red_CF_par_omega = Q_red_CF_par_omega./Q_fuel_CF_par_omega;        % Reduction specific energy [J/mol-MO] (CF)
    F_red_CF_par_T = Q_red_CF_par_T./Q_fuel_CF_par_T;                    % Reduction specific energy [J/mol-MO] (CF)
    if MFR_flag==1
        F_red_MFR_par_omega = Q_red_MFR_par_omega./Q_fuel_MFR_par_omega;     % Reduction specific energy [J/mol-MO] (MFR)
        F_red_MFR_par_T = Q_red_MFR_par_T./Q_fuel_MFR_par_T;                 % Reduction specific energy [J/mol-MO] (MFR)
    end
    %% Calculation of the exothermic heat of oxidation
    switch K_input
        case 1
            % Calculate enthalpy of oxidation reaction - T_ox
            Delta_H_r_par_omega = Reaction_Enthalpy_WS(T_ox);
            % Calculate enthalpy of oxidation reaction - par_T
            Delta_H_r_par_T = Reaction_Enthalpy_WS(T_ox_par);
        case 2
            % Calculate enthalpy of oxidation reaction - T_ox
            Delta_H_r_par_omega = Reaction_Enthalpy_CDS(T_ox);
            % Calculate enthalpy of oxidation reaction - par_T
            Delta_H_r_par_T = Reaction_Enthalpy_CDS(T_ox_par);
    end
    Q_ox_PF_par_omega = (-Q_red_PF_par_omega+delta_phi_PF_par_omega*Delta_H_r_par_omega)*eps_HR_ox;    % Exothermic heat recovered in [J/mol-MO]
    Q_ox_CF_par_omega = (-Q_red_CF_par_omega+delta_phi_CF_par_omega*Delta_H_r_par_omega)*eps_HR_ox;    % Exothermic heat recovered in [J/mol-MO]
    Q_ox_PF_par_T = (-Q_red_PF_par_T+delta_phi_PF_par_T.*Delta_H_r_par_T)*eps_HR_ox;    % Exothermic heat recovered in [J/mol-MO]
    Q_ox_CF_par_T = (-Q_red_CF_par_T+delta_phi_CF_par_T.*Delta_H_r_par_T)*eps_HR_ox;    % Exothermic heat recovered in [J/mol-MO]
    if MFR_flag==1
        Q_ox_MFR_par_omega = (-Q_red_MFR_par_omega+delta_phi_MFR_par_omega*Delta_H_r_par_omega)*eps_HR_ox;      % Exothermic heat recovered in [J/mol-MO]
        Q_ox_MFR_par_T = (-Q_red_MFR_par_T+delta_phi_MFR_par_T.*Delta_H_r_par_T)*eps_HR_ox;                     % Exothermic heat recovered in [J/mol-MO]
    end
    F_ox_PF_par_omega = Q_ox_PF_par_omega./Q_fuel_PF_par_omega;         % Ratio between oxidation energy and fuel energy (PF)
    F_ox_PF_par_T = Q_ox_PF_par_T./Q_fuel_PF_par_T;                     % Ratio between oxidation energy and fuel energy (PF)
    F_ox_CF_par_omega = Q_ox_CF_par_omega./Q_fuel_CF_par_omega;         % Ratio between oxidation energy and fuel energy (CF)
    F_ox_CF_par_T = Q_ox_CF_par_T./Q_fuel_CF_par_T;                     % Ratio between oxidation energy and fuel energy (CF)
    if MFR_flag==1
        F_ox_MFR_par_omega = Q_ox_MFR_par_omega./Q_fuel_MFR_par_omega;      % Ratio between oxidation energy and fuel energy (MFR)
        F_ox_MFR_par_T = Q_ox_MFR_par_T./Q_fuel_MFR_par_T;                  % Ratio between oxidation energy and fuel energy (MFR)
    end
    %% Calculation of sensible heating of sweep gas
    h_N2_0 = py.CoolProp.CoolProp.PropsSI('HMOLAR','T',T0,'P',p_red,'N2');          % N2 enthalpy at T0 [J/mol]
    h_N2_red = py.CoolProp.CoolProp.PropsSI('HMOLAR','T',T_red,'P',p_red,'N2');     % N2 enthalpy at T_red [J/mol]
    Q_sweep_h_par_omega = omega_red_par'.*(h_N2_red-h_N2_0);                        % Sweep gas required heating [J/mol-MO]
    Q_sweep_h_par_T_red = zeros(length(T_red_par),1);                               % Sweep gas required heating [J/mol-MO]
    for I=1:length(T_red_par)
        Q_sweep_h_par_T_red(I) = omega_red*(py.CoolProp.CoolProp.PropsSI('HMOLAR','T',T_red_par(I),'P',p_red,'N2')-h_N2_0);
    end
    F_sweep_h_PF_par_omega = Q_sweep_h_par_omega./Q_fuel_PF_par_omega;      % Ratio of sweep gas required heating to fuel energy (PF)
    F_sweep_h_PF_par_T = Q_sweep_h_par_T_red./Q_fuel_PF_par_T;              % Ratio of sweep gas required heating to fuel energy (PF)
    F_sweep_h_CF_par_omega = Q_sweep_h_par_omega./Q_fuel_CF_par_omega;      % Ratio of sweep gas required heating to fuel energy (CF)
    F_sweep_h_CF_par_T = Q_sweep_h_par_T_red./Q_fuel_CF_par_T;              % Ratio of sweep gas required heating to fuel energy (CF)
    if MFR_flag==1
        F_sweep_h_MFR_par_omega = Q_sweep_h_par_omega./Q_fuel_MFR_par_omega;    % Ratio of sweep gas required heating to fuel energy (MFR)
        F_sweep_h_MFR_par_T = Q_sweep_h_par_T_red./Q_fuel_MFR_par_T;            % Ratio of sweep gas required heating to fuel energy (MFR)
    end
    %% Calculation of F_total and process efficiency
    F_total_PF_par_omega_no_HR = F_sens_PF_par_omega+F_inert_PF_par_omega+F_red_PF_par_omega+F_ox_h_PF_par_omega+F_sweep_h_PF_par_omega+F_prod_sep_PF_par_omega;
    F_total_PF_par_T_no_HR = F_sens_PF_par_T+F_inert_PF_par_T+F_red_PF_par_T+F_ox_h_PF_par_T+F_sweep_h_PF_par_T+F_prod_sep_PF_par_T;
    F_total_CF_par_omega_no_HR = F_sens_CF_par_omega+F_inert_CF_par_omega+F_red_CF_par_omega+F_ox_h_CF_par_omega+F_sweep_h_CF_par_omega+F_prod_sep_CF_par_omega;
    F_total_CF_par_T_no_HR = F_sens_CF_par_T+F_inert_CF_par_T+F_red_CF_par_T+F_ox_h_CF_par_T+F_sweep_h_CF_par_T+F_prod_sep_CF_par_T;
    if MFR_flag==1
        F_total_MFR_par_omega_no_HR = F_sens_MFR_par_omega+F_inert_MFR_par_omega+F_red_MFR_par_omega+F_ox_h_MFR_par_omega+F_sweep_h_MFR_par_omega+F_prod_sep_MFR_par_omega;
        F_total_MFR_par_T_no_HR = F_sens_MFR_par_T+F_inert_MFR_par_T+F_red_MFR_par_T+F_ox_h_MFR_par_T+F_sweep_h_MFR_par_T+F_prod_sep_MFR_par_T;
    end
    % Efficiency - no heat recovery
    eta_PF_par_omega_no_HR = 1./F_total_PF_par_omega_no_HR;
    eta_PF_par_T_no_HR = 1./F_total_PF_par_T_no_HR;
    eta_CF_par_omega_no_HR = 1./F_total_CF_par_omega_no_HR;
    eta_CF_par_T_no_HR = 1./F_total_CF_par_T_no_HR;
    if MFR_flag==1
        eta_MFRF_par_omega_no_HR = 1./F_total_MFR_par_omega_no_HR;
        eta_MFR_par_T_no_HR = 1./F_total_MFR_par_T_no_HR;
    end
    % Plot efficiency (PF) - omega
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,eta_PF_par_omega_no_HR'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,eta_PF_par_omega_no_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta (no HR)','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_eta_no_HR_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (CF) - omega
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,eta_CF_par_omega_no_HR'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,eta_CF_par_omega_no_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta (no HR)','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_eta_no_HR_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (PF) - T
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,eta_PF_par_T_no_HR'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,eta_PF_par_T_no_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta {no HR}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_eta_no_HR_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (CF) - T
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,eta_CF_par_T_no_HR'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,eta_CF_par_T_no_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta {no HR}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_eta_no_HR_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    %% Efficiency - with heat recovery
    F_total_PF_par_omega_HR = (F_prod_sep_PF_par_omega+F_inert_PF_par_omega+min(F_ox_h_PF_par_omega+F_ox_PF_par_omega,0)*eta_ox_htw)+...
        F_sens_PF_par_omega*(1-eps_HR)+F_red_PF_par_omega+F_sweep_h_PF_par_omega*(1-eps_g)+max(F_ox_h_PF_par_omega+F_ox_PF_par_omega,0);
    F_total_CF_par_omega_HR = (F_prod_sep_CF_par_omega+F_inert_CF_par_omega+min(F_ox_h_CF_par_omega+F_ox_CF_par_omega,0)*eta_ox_htw)+...
        F_sens_CF_par_omega*(1-eps_HR)+F_red_CF_par_omega+F_sweep_h_CF_par_omega*(1-eps_g)+max(F_ox_h_CF_par_omega+F_ox_CF_par_omega,0);
    F_total_PF_par_T_HR = (F_prod_sep_PF_par_T+F_inert_PF_par_T+min(F_ox_h_PF_par_T+F_ox_PF_par_T,0)*eta_ox_htw)+...
        F_sens_PF_par_T*(1-eps_HR)+F_red_PF_par_T+F_sweep_h_PF_par_T*(1-eps_g)+max(F_ox_h_PF_par_T+F_ox_PF_par_T,0);
    F_total_CF_par_T_HR = (F_prod_sep_CF_par_T+F_inert_CF_par_T+min(F_ox_h_CF_par_T+F_ox_CF_par_T,0)*eta_ox_htw)+...
        F_sens_CF_par_T*(1-eps_HR)+F_red_CF_par_T+F_sweep_h_CF_par_T*(1-eps_g)+max(F_ox_h_CF_par_T+F_ox_CF_par_T,0);
    if MFR_flag==1
        F_total_MFR_par_omega_HR = (F_prod_sep_MFR_par_omega+F_inert_MFR_par_omega+min(F_ox_h_MFR_par_omega+F_ox_MFR_par_omega,0)*eta_ox_htw)+...
            F_sens_MFR_par_omega*(1-eps_HR)+F_red_MFR_par_omega+F_sweep_h_MFR_par_omega*(1-eps_g)+max(F_ox_h_MFR_par_omega+F_ox_MFR_par_omega,0);
        F_total_MFR_par_T_HR = (F_prod_sep_MFR_par_T+F_inert_MFR_par_T+min(F_ox_h_MFR_par_T+F_ox_MFR_par_T,0)*eta_ox_htw)+...
            F_sens_MFR_par_T*(1-eps_HR)+F_red_MFR_par_T+F_sweep_h_MFR_par_T*(1-eps_g)+max(F_ox_h_MFR_par_T+F_ox_MFR_par_T,0);
    end
    eta_PF_par_omega_HR = 1./F_total_PF_par_omega_HR;
    eta_PF_par_T_HR = 1./F_total_PF_par_T_HR;
    eta_CF_par_omega_HR = 1./F_total_CF_par_omega_HR;
    eta_CF_par_T_HR = 1./F_total_CF_par_T_HR;
    if MFR_flag==1
        eta_MFR_par_omega_HR = 1./F_total_MFR_par_omega_HR;
        eta_MFR_par_T_HR = 1./F_total_MFR_par_T_HR;
    end
    % Plotting efficiencies with HR
    % Plot efficiency (PF) - omega with HR
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,eta_PF_par_omega_HR'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,eta_PF_par_omega_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta (HR)','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_eta_HR_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (CF) - omega with HR
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(omega_red_par,omega_ox_par,eta_CF_par_omega_HR'); 
    hold on;
    shading interp;
    cm = colormap(jet); 
    cb = colorbar;
    xlabel(ax,'\omega_{red}','FontSize',fontsize_l);
    ylabel(ax,'\omega_{ox}','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(omega_red_par,omega_ox_par,eta_CF_par_omega_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta (HR)','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_Omega_Map_eta_HR_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (PF) - T with HR
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,eta_PF_par_T_HR'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,eta_PF_par_T_HR','LineColor','k','ShowText','on','LabelSpacing',150); %#ok<ASGLU>
    ylabel(cb,'\eta {HR}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'PF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_eta_HR_PF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
    % Plot efficiency (CF) - T with HR
    fig = figure;
    ax = gca;
    ax.FontSize = fontsize_a;
    s = pcolor(T_red_par,T_ox_par,eta_CF_par_T_HR'); 
    hold on;
    shading interp;
    colormap(jet);
    cb = colorbar;
    xlabel(ax,'T_{red}, K','FontSize',fontsize_l);
    ylabel(ax,'T_{ox}, K','FontSize',fontsize_l);
    ax.XScale = 'log';
    ax.YScale = 'log';
    [Mn,cn] = contour(T_red_par,T_ox_par,eta_CF_par_T_HR','LineColor','k','ShowText','on','LabelSpacing',150); 
    ylabel(cb,'\eta {HR}','FontSize',fontsize_l);
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.Color = 'k';
    ax.YAxis.LineWidth = 1.5;
    dummy = line(NaN,NaN,'LineStyle','none','Marker','none','Color','none');
    legend(dummy,'CF','FontWeight','bold','FontSize',fontsize_l,'Location','northwest');
    if SAVEFLAG
        fig_name = 'Par_T_Map_eta_HR_CF';
        savefig(fig,fullfile(folder_str,fig_name));
        saveas(fig,fullfile(folder_str,fig_name),'epsc');
        saveas(fig,fullfile(folder_str,fig_name),'emf');
        print(fig,fullfile(folder_str,fig_name),'-r1000','-dpng');
    end
end
%% Save processed data
if SAVEFLAG
    % Get a list of all variables
    allvars = whos;
    % Identify the variables that ARE NOT graphics handles. This uses a regular
    % expression on the class of each variable to check if it's a graphics object
    tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
    % Pass these variable names to save
    save(fullfile(path_str,strcat(filename_str,'_processed')), allvars(tosave).name);
end
%% Functions
function result = CP_PropsSI(varargin)
    % Shorthand version of CoolProp for MATLAB
    result = py.CoolProp.CoolProp.PropsSI(varargin{:});
end