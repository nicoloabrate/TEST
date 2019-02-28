function []=display_solution(fd_data,fd_solution,plot_opt,layer_label)
% Display solution computed by TEST code. This function can plot:
% 1) Fundamental group-wise solution
% 2) Higher order harmonics and eigenvalues
% 3) Neutron and importance balance histograms
% 4) Fundamental flux/adjoint distribution and nuclear density profile, if
% any

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plotting the flux shapes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Solution structure
k_eff     = fd_solution.k_eff;
k_inf     = fd_solution.k_inf;
k_n       = fd_solution.k_n;
NG        = fd_solution.NG;
NB_phi    = fd_solution.NB_phi;
NB_adj    = fd_solution.NB_adj;
nev       = size(fd_solution.EVd,2);

% Grid structure
fdgrid    = fd_data.FD_grid;
t_reg     = fd_data.t_reg;
geom_type = fd_data.geom;
Nt        = length(fdgrid);

%% --------------------------------- Solution manipulation ------------------------------------------------------------------------------------------------------------

% Flux
phi      = fd_solution.EVd;
phi(:,1) = abs(phi(:,1)); % fundamental harmonics is always positive
% Adjoint
adj      = fd_solution.EVa;
adj(:,1) = abs(adj(:,1)); % fundamental harmonics is always positive

switch plot_opt.NC
    case 1
        phi = phi/max(phi(:,1));
        adj = adj/max(adj(:,1));
    case 2
        phi = phi/NB_phi.TOT_FISS;
        adj = adj/NB_adj.TOT_FISS;
end

% Geometry label definition
switch geom_type
    case 1
        geom_label = 'slab';
    case 2
        geom_label = 'sphere';
    case 3
        geom_label = 'cylinder';
end


%% --------------------------------- Plotting solution ------------------------------------------------------------------------------------------------------------
dmax     = zeros(NG,1);
amax     = zeros(NG,1);
dmin     = zeros(NG,1);
amin     = zeros(NG,1);
leg_strd = cell(1,NG);
leg_stra = cell(1,NG);

% Curve aspect 
col  = brewermap(NG+10,'Spectral');
lin  = {'-','-.'};
mark = {'none','o','^','v','s','x','+','*'};


% Marker indeces
N_reg = length(t_reg)-1;
xmark = zeros(N_reg*5,1);
N_old = 0;
N     = fd_data.N;
for ii=1:N_reg
   xmark(1+(ii-1)*5:5+(ii-1)*5) = ceil(linspace(1+N_old,N(ii)+N_old,5));
   N_old                        = N_old+N(ii);
end
xmark = round(xmark);

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',22)
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','latex');

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'Color', 'w');
hold on
box on

yyaxis left
set(gca,'ycolor','k')

% Plot group-wise fluxes
for gg = 1:NG
    
    ii = gg;
    jj = 1;
    kk = gg;
    
    plot(fdgrid,phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1),'color',col(ii,:),'linestyle',lin{jj},'Marker',mark{kk},...
         'MarkerEdgeColor',col(ii,:),'MarkerFaceColor',col(ii,:),'MarkerSize',5,'MarkerIndices',xmark);                
    leg_strd{gg} = ['$\phi_',num2str(gg),'$'];
    dmax(gg) = max(phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    dmin(gg) = min(phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    
end

ylim([min(dmin) max(dmax)])

if plot_opt.plot_scale == 1
    plot_scale = 'lin'; 
else
    plot_scale = 'log';
end

set(gca,'XMinorTick','on','xminorgrid','on','yminorgrid','on','YScale',plot_scale);

% X and Y labels
if fd_data.geom==1
    ylabel('$\phi(\rm x)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('x [cm]','interpreter','latex')
else
    ylabel('$\phi(\rm r)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('r [cm]','interpreter','latex')
end

yyaxis right
set(gca,'ycolor','k')

% Plot group-wise adjoints
for gg = 1:NG
    
    ii = gg;
    jj = 2;
    kk = gg;
    
    plot(fdgrid,adj(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1),'color',col(end-ii,:),'linestyle',lin{jj},'Marker',mark{kk},...
                    'MarkerEdgeColor',col(end-ii,:),'MarkerFaceColor',col(end-ii,:),'MarkerSize',5,'MarkerIndices',xmark);
                
    leg_stra{gg} = ['$\phi^+_',num2str(gg),'$'];
    amax(gg) = max(adj(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    amin(gg) = min(adj(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    
end

set(gca,'XMinorTick','on','xminorgrid','on','yminorgrid','on','YScale',plot_scale);
ylim([min(amin) max(amax)])

% X and Y labels
if fd_data.geom==1
    ylabel('$\phi^+(\rm x)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('x [cm]','interpreter','latex')
else
    ylabel('$\phi^+(\rm r)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('r [cm]','interpreter','latex')
end

% Plotting region boundaries
if length(t_reg)>2
    for ii=2:length(t_reg)-1
        plot([t_reg(ii) t_reg(ii)],[0  max(amax)],'k--','linewidth',1.5);
        hold on
    end
end


% Set legend specs
hold off
lgd = legend([leg_strd leg_stra],'location','bestoutside');
lgd.FontSize = 20;
grid on

title([geom_label,': $k_{\rm \infty}=',num2str(k_inf),'$, $k_{\rm eff}=',num2str(k_eff),'$'],'interpreter','latex')

% Define figure name to be saved

if plot_opt.plot_save == 1
    model_label = 'flux';
    fig_name = [num2str(model_label),'_',num2str(plot_opt.plot_scale),'_'];
    fig_name2 = 'neut_balance_';
    fig_name3 = 'imp_balance_';
    for ii=1:fd_data.N_reg
        if ii==fd_data.N_reg
            fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
        else
            fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
        end
    end
    export_fig([pwd,'/output/plot/flux_adjoint/slab/',fig_name,'.pdf'])
    export_fig([pwd,'/output/plot/flux_adjoint/slab/',fig_name,'.fig'])
end

%% Plot flux and density profile
% Nuclear density spatial behaviour
if plot_opt.showplot.dens_profile == 1
adens_shape = [];
N_old = 0;
for ii = 1:fd_data.N_reg
   x = 1+N_old:N_old+fd_data.N(ii);
   vect = fd_data.f{ii}(fdgrid(x));
   if length(vect) == length(x)
       adens_shape=[adens_shape vect];
   else
       vect = vect*ones(1,length(x));
       adens_shape=[adens_shape vect];
   end
   N_old = N_old+fd_data.N(ii);
end

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',22)
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','latex');

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf, 'Color', 'w');
hold on
box on

yyaxis left
set(gca,'ycolor','k')

p1 = plot(fdgrid,adens_shape,'r');
aa = area(fdgrid,adens_shape);
set(aa,'FaceColor','r','FaceAlpha',0.15);

hold on
grid on
if fd_data.geom==1
    ylabel('$\rm N(x)~~[nucleii/cm^3]$','interpreter','latex')
else
    ylabel('$\rm N(r)~~[nucleii/cm^3]$','interpreter','latex')
end

title('Atomic density profile')


yyaxis right
set(gca,'ycolor','k')

% Plot group-wise fluxes
for gg = 1:NG
    
    ii = gg;
    jj = 1;
    kk = gg;
    
    h1(gg) = plot(fdgrid,phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1),'color',col(ii,:),'linestyle',lin{jj},'Marker',mark{kk},...
                    'MarkerEdgeColor',col(ii,:),'MarkerFaceColor',col(ii,:),'MarkerSize',5,'MarkerIndices',xmark);
                
    leg_strd{gg} = ['$\phi_',num2str(gg),'$'];
    dmax(gg) = max(phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    dmin(gg) = min(phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,1));
    
end
% Plotting region boundaries
if length(t_reg)>2
    for ii=2:length(t_reg)-1
        plot([t_reg(ii) t_reg(ii)],[min(dmin) max(dmax)],'r--','linewidth',1);
        hold on
    end
end
ylim([min(dmin) max(dmax)])
% X and Y labels
if fd_data.geom==1
    ylabel('$\phi_g(\rm x)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('x [cm]','interpreter','latex')
else
    ylabel('$\phi_g(\rm r)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
    xlabel('r [cm]','interpreter','latex')
end
lgd = legend(h1,leg_strd,'location','bestoutside');
lgd.FontSize = 20;

end

if plot_opt.plot_save == 1
    model_label = 'flux_adens_profile';
    fig_name = [num2str(model_label),'_',num2str(plot_opt.plot_scale),'_'];
    fig_name2 = 'neut_balance_';
    fig_name3 = 'imp_balance_';
    for ii=1:fd_data.N_reg
        if ii==fd_data.N_reg
            fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
        else
            fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
        end
    end
    export_fig([pwd,'\output\plot\slab\',fig_name,'.fig'])
end

%% Plot harmonics 
if plot_opt.showplot.harmonics == 1
    
    col  = brewermap(nev,'Spectral');

    % Marker indeces
    N_reg = length(t_reg)-1;
    xmark = zeros(N_reg*5,1);
    N_old = 0;
    N     = fd_data.N;
    
    for ii=1:N_reg
        xmark(1+(ii-1)*5:5+(ii-1)*5) = linspace(1+N_old,N(ii)+N_old,5);
        N_old                        = N(ii)+N_old;
    end
    
    xmark = round(xmark);
    
    set(0,'DefaultLineLineWidth',2)
    set(0,'DefaultAxesFontSize',22)
    set(groot, 'defaultAxesTickLabelInterpreter','none');
    set(groot, 'defaultLegendInterpreter','latex');
     
    % Plot group-wise fluxes
    for gg = 1:NG
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf, 'Color', 'w');
        hold on
        box on

        for nn = 1:nev
            
            ii = gg;
            jj = 1;
            kk = gg;
            phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,nn) = phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,nn)/max(phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,nn));

            plot(fdgrid,phi(1+(gg-1)*Nt:1:Nt+(gg-1)*Nt,nn),'color',col(nn,:),'linestyle',lin{jj},'Marker',mark{kk},...
                'MarkerEdgeColor',col(nn,:),'MarkerFaceColor',col(nn,:),'MarkerSize',5,'MarkerIndices',xmark);
            
            leg_strd{nn} = ['$\phi_{',num2str(gg),',',num2str(nev-1),'}~~k_',num2str(nn-1),'=',num2str(k_n(nn)),'$'];
        end
        
    if plot_opt.plot_scale == 1
        plot_scale = 'lin';
    else
        plot_scale = 'log';
    end
    
    set(gca,'XMinorTick','on','xminorgrid','on','yminorgrid','on','YScale',plot_scale);
    
    % X and Y labels
    if fd_data.geom==1
        ylabel('$\phi(\rm x)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
        xlabel('x [cm]','interpreter','latex')
    else
        ylabel('$\phi(\rm r)$ [cm \textsuperscript{-2}s \textsuperscript{-1}]','interpreter','latex')
        xlabel('r [cm]','interpreter','latex')
    end
    
    % Plotting region boundaries
    if length(t_reg)>2
        for ii=2:length(t_reg)-1
            plot([t_reg(ii) t_reg(ii)],[-1  1],'k--','linewidth',1.5);
            hold on
        end
    end
    
    
    % Set legend specs
    hold off
    lgd = legend(leg_strd,'location','bestoutside');
    lgd.FontSize = 20;
    grid on
    
    title([geom_label,': $k_{\rm \infty}=',num2str(k_inf),'$, $k_{\rm eff}=',num2str(k_eff),'$'],'interpreter','latex')
    
    % Define figure name to be saved
    
    if plot_opt.plot_save == 1
        model_label = 'flux_harmonics';
        fig_name = [num2str(model_label),'_',num2str(plot_opt.plot_scale),'_'];
        fig_name2 = 'neut_balance_';
        fig_name3 = 'imp_balance_';
        for ii=1:fd_data.N_reg
            if ii==fd_data.N_reg
                fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
                fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
                fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            else
                fig_name  = [fig_name,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
                fig_name2 = [fig_name2,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
                fig_name3 = [fig_name3,layer_label{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            end
        end
        export_fig([pwd,'\output\plot\',fig_name,'.pdf'])
        export_fig([pwd,'\output\plot\',fig_name,'.fig'])
    end
    
    end
    
end



%% ---------------------------------- Neutron balance plot -----------------------------------
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',16)
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','latex');
colmap = brewermap(7,'Spectral');

if NG == 2

if plot_opt.showplot.balance==1

    for ii = 1:fd_data.N_reg
        
        % xlim evaluation
        xmin = min(min([NB_phi.FISS;NB_phi.IN_SCATT;-NB_phi.ABS_FISS;-NB_phi.ABS_CAPT;-NB_phi.OUT_SCATT;-NB_phi.LEAK_R;-NB_phi.LEAK_L(ii,1:NG);]));
        xmax = max(max([NB_phi.FISS;NB_phi.IN_SCATT;-NB_phi.ABS_FISS;-NB_phi.ABS_CAPT;-NB_phi.OUT_SCATT;-NB_phi.LEAK_R;-NB_phi.LEAK_L(ii,1:NG);]));

        figure('units','normalized','outerposition',[0 0 0.7 0.7]);
        set(gcf, 'Color', 'w');
        % IF to distinguish fissile from non-fissile regions
        if sum(NB_phi.FISS(ii,1:NG))>0
            X = [NB_phi.FISS(ii,1:NG);NB_phi.IN_SCATT(ii,1:NG);-NB_phi.ABS_FISS(ii,1:NG);-NB_phi.ABS_CAPT(ii,1:NG);-NB_phi.OUT_SCATT(ii,1:NG);-NB_phi.LEAK_R(ii,1:NG);...
                -NB_phi.LEAK_L(ii,1:NG);];
            % Possible legend entries
            lgd_entries = {'Fission','In scattering','Absorp (fission)','Absorp (captures)','Out scattering',...
                           ['net transfer at x=',num2str(fd_data.t_reg(ii+1)),' cm'],['net transfer at  x=',num2str(fd_data.t_reg(ii))]};
        else % Non-fissile region
            X = [NB_phi.IN_SCATT(ii,1:NG);-NB_phi.ABS_CAPT(ii,1:NG);-NB_phi.OUT_SCATT(ii,1:NG);-NB_phi.LEAK_R(ii,1:NG);...
                -NB_phi.LEAK_L(ii,1:NG);];
            % Possible legend entries
            lgd_entries = {'In scattering','Absorp (captures)','Out scattering',...
                           ['net transfer at x=',num2str(fd_data.t_reg(ii+1)),' cm'],['net transfer at  x=',num2str(fd_data.t_reg(ii))]};
        end
        
        % Exchange legend and array to have same colour on plot for same
        % currents (e.g. J in x0 == J_in in x0+eps)
        if mod(ii,2) == 1
           dummy_lgd          = lgd_entries(end); 
           lgd_entries(end)   = lgd_entries(end-1);
           lgd_entries(end-1) = dummy_lgd;
           dummy_X            = X(end,:);
           X(end,:)           = X(end-1,:);
           X(end-1,:)         = dummy_X;
        end
        % Horizontal bar plot
        barh(X','hist');
        if sum(NB_phi.FISS(ii,1:NG))>0
            colormap(gca,colmap)
        else
            colormap(gca,[colmap(2,:);colmap(4:end,:)])
        end
        lgd2 = legend(lgd_entries,'location','bestoutside');
        lgd2.FontSize = 15;
        grid on
        title([layer_label{ii},' region'])
        ylabel('Energy group')
        yticklabels({'Group 1','Group 2'})
        xlabel('normalized wrt $<\nu\Sigma_{f}\phi>$','interpreter','latex')
        xlim([xmin xmax])
        set(gca,'XMinorTick','on','xminorgrid','on','ygrid','off');
        
        % Saving plots
        if plot_opt.plot_save == 1
            export_fig([pwd,'/output/plot/balance/',geom_label,'/',layer_label{ii},'_',num2str(t_reg(ii)),'_',num2str(t_reg(ii+1)),'.pdf'])
            export_fig([pwd,'/output/plot/balance/',geom_label,'/',layer_label{ii},'_',num2str(t_reg(ii)),'_',num2str(t_reg(ii+1)),'.fig'])
        end
        
    end
    
end




if plot_opt.showplot.balance_adj==1
    
    set(0,'DefaultLineLineWidth',1.5)
    set(0,'DefaultAxesFontSize',16)
    set(groot, 'defaultAxesTickLabelInterpreter','none');
    set(groot, 'defaultLegendInterpreter','latex');
    
    colmap = brewermap(7,'Spectral');
    for ii = 1:fd_data.N_reg
        
        % xlim evaluation
        xmin = min(min([NB_adj.FISS;NB_adj.IN_SCATT;-NB_adj.ABS_FISS;-NB_adj.ABS_CAPT;-NB_adj.OUT_SCATT;-NB_adj.LEAK_R;-NB_adj.LEAK_L(ii,1:NG);]));
        xmax = max(max([NB_adj.FISS;NB_adj.IN_SCATT;-NB_adj.ABS_FISS;-NB_adj.ABS_CAPT;-NB_adj.OUT_SCATT;-NB_adj.LEAK_R;-NB_adj.LEAK_L(ii,1:NG);]));
        
        
        figure('units','normalized','outerposition',[0 0 0.7 0.7]);
        set(gcf, 'Color', 'w');
        % IF to distinguish fissile from non-fissile regions
        if sum(NB_adj.FISS(ii,1:NG))>0
            X = [NB_adj.FISS(ii,1:NG);NB_adj.IN_SCATT(ii,1:NG);-NB_adj.ABS_FISS(ii,1:NG);-NB_adj.ABS_CAPT(ii,1:NG);-NB_adj.OUT_SCATT(ii,1:NG);-NB_adj.LEAK_R(ii,1:NG);...
                -NB_adj.LEAK_L(ii,1:NG);];
            % Possible legend entries
            lgd_entries = {'Fission','In scattering','Absorp (fission)','Absorp (captures)','Out scattering',...
                           ['net transfer at x=',num2str(fd_data.t_reg(ii+1)),' cm'],['net transfer at  x=',num2str(fd_data.t_reg(ii))]};
        else % Non-fissile region
            X = [NB_adj.IN_SCATT(ii,1:NG);-NB_adj.ABS_CAPT(ii,1:NG);-NB_adj.OUT_SCATT(ii,1:NG);-NB_adj.LEAK_R(ii,1:NG);...
                -NB_adj.LEAK_L(ii,1:NG);];
            % Possible legend entries
            lgd_entries = {'In scattering','Absorp (captures)','Out scattering',...
                           ['net transfer at x=',num2str(fd_data.t_reg(ii+1)),' cm'],['net transfer at  x=',num2str(fd_data.t_reg(ii))]};
        end
        
        % Exchange legend and array to have same colour on plot for same
        % currents (e.g. J in x0 == J_in in x0+eps)
        if mod(ii,2) == 1
            dummy_lgd          = lgd_entries(end);
            lgd_entries(end)   = lgd_entries(end-1);
            lgd_entries(end-1) = dummy_lgd;
            dummy_X            = X(end,:);
            X(end,:)           = X(end-1,:);
            X(end-1,:)         = dummy_X;
        end
        
        % Horizontal bar plot
        barh(X','hist');
        if sum(NB_adj.FISS(ii,1:NG))>0
            colormap(gca,colmap)
        else
            colormap(gca,[colmap(2,:);colmap(4:end,:)])
        end
        lgd2 = legend(lgd_entries,'location','bestoutside');
        lgd2.FontSize = 15;
        grid on
        title([layer_label{ii},' region'])
        ylabel('Energy group')
        yticklabels({'Group 1','Group 2'})
        xlabel('normalized wrt $<\nu\Sigma_{f}\phi>$','interpreter','latex')
        xlim([xmin xmax])
        set(gca,'XMinorTick','on','xminorgrid','on','ygrid','off');
        
        
        % Saving plots
        if plot_opt.plot_save == 1
            export_fig([pwd,'/output/plot/balance_adj/',geom_label,'/',layer_label{ii},'_',num2str(t_reg(ii)),'_',num2str(t_reg(ii+1)),'.pdf'])
            export_fig([pwd,'/output/plot/balance_adj/',geom_label,'/',layer_label{ii},'_',num2str(t_reg(ii)),'_',num2str(t_reg(ii+1)),'.fig'])
        end
        
    end
    
end
  
end

end
