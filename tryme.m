clear all
close all
clc

% This script calls the basic TEST functions.
addpath(genpath('main'))  
addpath(genpath('auxiliary'))
%% --------------------------------- Settings ------------------------------------------------------------------------------------------------------------
% Numerical discretization setup
BC = [2 2];                % boundary conditions: set 1 for black, 2 for reflective
NP = 2000;                 % number of spatial points per diffusion length
NG = 2;                    % number of energy group
PN = 1;                    % order of Legendre polynomials expansion
NS = 50;                   % number of points to sample diffusion length when profile is non-constant

% Model problem setup
geom_type  = 3;            % type 1 for slab, 2 for sphere, 3 for cylinder
geom_label = {'slab','sphere','cyl'};

% Post-processing options
flag_plot  =  1;           % type 1 to plot flux/adjoint shapes
NC_type    =  2;           % 2 to normalize flux and adjoint over total fission, 1 to normalize over maximum
plot_scale =  1;           % type 1 for linear scale, 2 for semilogy
plot_save  =  0;           % type 1 to save the plot
mat_save   = -1;           % type 1 save the main variables in workspace
verbosity  =  10;          % type 1 to print output diagnostics

% System parameters
t_reg  = [0,40,60];        % layers/shells coordinates
N_reg  = length(t_reg)-1;   % number of regions of the reactor
t_core = t_reg(end);        % core thickness

% Properties spatial shape
f_def   = {'1','1','1','1','1'};
f_shape = adens_shape(N_reg,t_reg,f_def,geom_type);
%% --------------------------------- Material regions ------------------------------------------------------------------------------------------------------------
% Define each material region assigning a set of material properties and a material id. 
T = 300;            % Evaluation temperature for XS 
%  Core 1st layer
mat_core  = 'UO2';
mat_id    = 1;
core_data = readnucdata(mat_core,NG,mat_id,T);

%  Core 2nd layer
mat_refl  = 'H2O';
mat_id    = 2;
refl_data = serp2struct(mat_refl,mat_id,NG);

% Gathering different region data
multig_data = {core_data,refl_data};  % AUTOMATIZZARE CONCATENAZIONE
layer_label = {mat_core,mat_refl};

%% --------------------------------- Spatial mesh generation ------------------------------------------------------------------------------------------------------------
% Generating the numerical grid for spatial discretization
mesh_opt     = struct('NG',NG,'NP',NP,'BC',BC,'NS',NS);
layer_struct = struct('N_reg',N_reg,'t_reg',t_reg,'geom_type',geom_type);
% The grid needs to be generated on the perturbed system
fd_data      = space_grid(mesh_opt,layer_struct,multig_data,f_shape);
Nt           = sum(fd_data.N,1);
%% --------------------------------- Model numerical approximation ------------------------------------------------------------------------------------------------------------
% Assemblying the operators
[L,F]         = MG(NG,PN,multig_data,fd_data);
[L_adj,F_adj] = MG(-NG,PN,multig_data,fd_data); % -NG set the adjoint operators construction

%% --------------------------------- Numerical solution (IRA) ------------------------------------------------------------------------------------------------------------
nev = 4; % Desired number of spatial modes

% Direct solution
[EVd,KN] = eigs(L,F,nev,'SM');
KN       = diag(KN);
k_n    = 1./KN;
k_eff    = k_n(1);

% Adjoint solution
[EVa,KN]  = eigs(L_adj,F_adj,nev,'SM');
KN        = diag(KN);
k_n_adj = 1./KN;
k_eff_adj = k_n_adj(1);

%% Reference eigenvalues
k_inf = MG_k_inf(NG,multig_data,1);

%% -------------------------------- Diagnostic outputs ------------------------------------------------------------------------------------------------
% Questa parte necessita di revisioni per automatizzare plot e output a
% schermo, ma per il momento, visto che funziona, la lascio cosi!

% Neutron balance
phi = EVd(:,1);
phi = abs(phi); % fundamental harmonics is always positive
adj = EVa(:,1);
adj = abs(adj);

% Call neutron balance evaluation (to check result accuracy and physical conservation of particles)
NB_phi = neutron_balance(fd_data,multig_data,phi,k_eff);
NB_adj = neutron_balance(fd_data,multig_data,adj,-k_eff);


if verbosity>1
    fprintf('Numerical calculation: \r\n')
    fprintf('k_inf = %f \r\n',k_inf)
    fprintf('k_eff = %f \r\n',k_eff)
    fprintf('The difference (k_eff-k_inf) is %f [pcm] \r\n',(k_n(1)-k_inf)*1e5)
    if verbosity>5
        fprintf('-- NEUTRON BALANCE \r\n')
        for ii=1:N_reg
            fprintf('Neutron balance in %s layer: \r\n',layer_label{ii})
            
            for g=1:NG
                fprintf('_________________________________________________ \r\n')
                fprintf('  Energy group %d: \r\n',g)
                fprintf('   Fissions:                       %e \r\n',NB_phi.FISS(ii,g))
                fprintf('   Absorption (fissions):          %e \r\n',-NB_phi.ABS_FISS(ii,g))
                fprintf('   Absorption (capture):           %e \r\n',-NB_phi.ABS_CAPT(ii,g))
                fprintf('   In-group scattering:            %e \r\n',NB_phi.IN_SCATT(ii,g))
                fprintf('   Out-group scattering:           %e \r\n',-NB_phi.OUT_SCATT(ii,g))
                fprintf('   Right boundary current:         %e \r\n',-NB_phi.LEAK_R(ii,g))
                fprintf('   Left boundary current:          %e \r\n',-NB_phi.LEAK_L(ii,g))
                fprintf('   In layer group imbalance:       %e \r\n',NB_phi.NG_IMBAL(ii,g))
                fprintf('________________________________ \r\n')
            end
            fprintf('Total imbalance:   %e \r\n',NB_phi.TOT_IMBAL)
            fprintf('_________________________________________________ \r\n')
        end
        
        if verbosity>10
            fprintf('-- IMPORTANCE BALANCE \r\n')
            for ii=1:N_reg
                fprintf('Importance balance in %s layer: \r\n',layer_label{ii})
                
                for g=1:NG
                    fprintf('_________________________________________________ \r\n')
                    fprintf('  Energy group %d: \r\n',g)
                    fprintf('   Fissions:                       %e \r\n',NB_adj.FISS(ii,g))
                    fprintf('   Absorption (fissions):          %e \r\n',-NB_adj.ABS_FISS(ii,g))
                    fprintf('   Absorption (capture):           %e \r\n',-NB_adj.ABS_CAPT(ii,g))
                    fprintf('   In-group scattering:            %e \r\n',NB_adj.IN_SCATT(ii,g))
                    fprintf('   Out-group scattering:           %e \r\n',-NB_adj.OUT_SCATT(ii,g))
                    fprintf('   Right boundary current:         %e \r\n',-NB_adj.LEAK_R(ii,g))
                    fprintf('   Left boundary current:          %e \r\n',-NB_adj.LEAK_L(ii,g))
                    fprintf('   In layer group imbalance:       %e \r\n',NB_adj.NG_IMBAL(ii,g))
                    fprintf('________________________________ \r\n')
                end
                fprintf('Total imbalance:   %e \r\n',NB_adj.TOT_IMBAL)
                fprintf('_________________________________________________ \r\n')
            end
        end
    end
end

%% --------------------------------- Plot ------------------------------------------------------------------------------------------------------------
if flag_plot > 0
    harmonics    = 0;
    balance      = 0;
    balance_adj  = 0;
    dens_profile = 0;
    showplot     = struct('harmonics',harmonics,'dens_profile',dens_profile,'balance',balance,'balance_adj',balance_adj);
    plot_opt     = struct('plot_save',plot_save,'plot_scale',plot_scale,'verbosity',verbosity,'NC',NC_type,'showplot',showplot);
    fd_solution = struct('EVd',EVd,'EVa',EVa,'k_n',k_n,'k_n_adj',k_n_adj,'k_eff',k_n(1),'k_inf',k_inf,'NG',NG,'NB_phi',NB_phi,'NB_adj',NB_adj);
    display_solution(fd_data,fd_solution,plot_opt,layer_label)
end

if mat_save>0
    if NG>0
        model_label = 'flux';
    else
        model_label = 'importance';
    end

    save([geometry,'_layers_',num2str(N_reg),'.mat'],'k_eff','fd_data','phi','model_label')
end