clear all
close all
clc

% This script calls the basic TEST functions.
addpath(genpath('main'))
addpath(genpath('auxiliary'))
%% --------------------------------- Settings ------------------------------------------------------------------------------------------------------------
% define geometric parameters
H  = 150;  % [cm] domain thickness
dx = 10;   % [cm] perturbation width
x0 = 75;   % [cm] perturbation center coordinate

N = 20;           % perturbation order
M = 100;          % expansion order (number of reference system harmonics)
% pert. [3 2 0 110] a x0=10,50,75 con dx=2.5,5,10

% Perturbation chain options
% pert_mat   = [1 2 0 1;   % fission XS, in-group, out-group, pert. amplitude
%               2 1 0 1;   % capture XS, in-group, out-group, pert. amplitude
%               2 2 0 1;   % capture XS, in-group, out-group, pert. amplitude
%               3 1 0 1;   % nubar, in-group, out-group, pert. amplitude
%               3 2 0 1;   % nubar, in-group, out-group, pert. amplitude
%               4 1 2 1;]; % scattering XS, in-group, out-group, pert. amplitude

pert_mat   = [2 1 1 500]; 

ext        = 'glob';     % perturbation spatial extension
unpert     = 1;
pert       = 1;
store_eig  = 1;
nev        = M;         % Desired number of spatial modes

% Numerical discretization setup
BC = [1 1];                % boundary conditions: set 1 for black, 2 for reflective
NP = 25;                   % number of spatial points per diffusion length
NG = 2;                    % number of energy group
PN = 1;                    % order of Legendre polynomials expansion
NS = 50;                   % number of points to sample diffusion length when profile is non-constant

% Model problem setup
geom_type  = 1;            % type 1 for slab, 2 for sphere, 3 for cylinder
geom_label = {'slab','sphere','cyl'};

% Post-processing options
flag_plot  =  1;           % type 1 to plot flux/adjoint shapes
NC_type    =  2;           % 2 to normalize flux and adjoint over total fission, 1 to normalize over maximum
plot_scale =  1;           % type 1 for linear scale, 2 for semilogy
plot_save  =  0;           % type 1 to save the plot
mat_save   = -1;           % type 1 save the main variables in workspace
verbosity  =  0;          % type 1 to print output diagnostics

% System parameters
t_reg  = [0,x0-dx/2,x0+dx/2,H];           % layers/shells coordinates
N_reg  = length(t_reg)-1;   % number of regions of the reactor
t_core = t_reg(end);        % core thickness

% Properties spatial shape
f_def   = {'1','1','1','1','1'};
f_shape = adens_shape(N_reg,t_reg,f_def,geom_type);
%% --------------------------------- Material regions ------------------------------------------------------------------------------------------------------------
% Define each material region assigning a set of material properties and a material id.
T = 300;            % Evaluation temperature for XS
%  Core 1st layer
mat_core  = 'AGPT';
mat_id    = 1;
core_data = readnucdata(mat_core,NG,mat_id,T);
pert_data = add_perturbation(pert_mat,core_data);
% Unperturbed system
multig_data_unp = {core_data,core_data,core_data};
layer_label_unp = {mat_core,mat_core,mat_core};
% Perturbed system
multig_data_per = {core_data,pert_data,core_data};
layer_label_per = {mat_core,mat_core,mat_core};

%% --------------------------------- Spatial mesh generation ------------------------------------------------------------------------------------------------------------
% Generating the numerical grid for spatial discretization
mesh_opt     = struct('NG',NG,'NP',NP,'BC',BC,'NS',NS);
layer_struct = struct('N_reg',N_reg,'t_reg',t_reg,'geom_type',geom_type);
% The grid needs to be generated on the perturbed system
fd_data      = space_grid(mesh_opt,layer_struct,multig_data_per,f_shape);
Nt           = sum(fd_data.N,1);
%% --------------------------------- Model numerical approximation ------------------------------------------------------------------------------------------------------------
% Assemblying the operators
% Unperturbed operators
[L,F]           = MG(NG,PN,multig_data_unp,fd_data);
[L_adj,F_adj]   = MG(-NG,PN,multig_data_unp,fd_data);
% Perturbed operators
[Lp,Fp]         = MG(NG,PN,multig_data_per,fd_data);
[Lp_adj,Fp_adj] = MG(-NG,PN,multig_data_per,fd_data);
% Operators perturbations
dL = Lp-L;
dF = Fp-F;
%% --------------------------------- Numerical solution (IRA) ------------------------------------------------------------------------------------------------------------
if unpert>0
    % -- Unperturbed system
    
    % Direct solution
    [EVd,KN] = eigs(L,F,nev,'SM');
    KN       = diag(KN);
    k_n      = 1./KN;
    k_eff    = k_n(1);
    SI=MG_SI(k_eff,NG,multig_data_unp,N_reg);
    % get the sign of the second row (to avoid BCs)
    signs = sign(EVd(2, :));
    % multiply all columns by the complex conjugate of sign of the first element:
    EVd = bsxfun(@times, EVd, conj(signs));
    % Normalize both forward and adjoint harmonics
    brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);
    for ii=1:nev
        EVd(:,ii) = EVd(:,ii)/sqrt(braket(EVd(:,ii),EVd(:,ii),brakopt));
    end
    % Adjoint solution
    [EVa,KN]  = eigs(L_adj,F_adj,nev,'SM');
    KN        = diag(KN);
    k_n_adj = 1./KN;
    k_eff_adj = k_n_adj(1);
    SIX=MG_SI(k_eff_adj,-NG,multig_data_unp,N_reg);
    % get the sign of the second row (to avoid BCs)
    signs = sign(EVa(2, :));
    % multiply all columns by the complex conjugate of sign of the first element:
    EVa = bsxfun(@times, EVa, conj(signs));
    % Normalize both forward and adjoint harmonics
    brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);
    for ii=1:nev
        EVa(:,ii) = EVa(:,ii)/sqrt(braket(EVa(:,ii),EVa(:,ii),brakopt));
    end
    
    % Balance evaluation
    NB_phi       = neutron_balance(fd_data,multig_data_unp,(EVd(:,1)),k_eff);
    NB_adj       = neutron_balance(fd_data,multig_data_unp,(EVa(:,1)),-k_eff);
    % Store eigenfunctions and operators
    if store_eig>0
        
        cname = '';
        for ii=1:N_reg
            if ii==N_reg
                cname  = [cname,layer_label_unp{ii},'_',num2str(t_reg(ii+1)-t_reg(ii))];
            else
                cname  = [cname,layer_label_unp{ii},'_',num2str(t_reg(ii+1)-t_reg(ii)),'_'];
            end
        end
        
        % Move to reactor directory
        try
            cd(cname)
        catch
            mkdir(cname)
        end
        save(['P_',num2str([1 2 3],'%g'),'_I_',num2str(pert_mat(4),'%g_'),ext,'_0.mat'],'L','F','L_adj','F_adj','EVd','EVa','KN','fd_data');
        cd ..
    end
end

% -- Perturbed system
% Direct solution
if pert > 0
    
    % Direct solution
    [EVd_per,KN_per] = eigs(Lp,Fp,1,'SM');
    KN_per           = diag(KN_per);
    k_n_per       = KN_per(1);
    k_n_per       = 1/k_n_per;
    % get the sign of the second row (to avoid BCs)
    signs = sign(EVd_per(2, :));
    % multiply all columns by the complex conjugate of sign of the first element:
    EVd_per = bsxfun(@times, EVd_per, conj(signs));
    % Normalize both forward and adjoint harmonics
    brakopt = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);
    for ii=1:1
        EVd_per(:,ii) = EVd_per(:,ii)/sqrt(braket(EVd_per(:,ii),EVd_per(:,ii),brakopt));
    end
    % Adjoint solution
    [EVa_per,KNa] = eigs(Lp_adj,Fp_adj,1,'SM');
    KNa           = diag(KNa);
    k_n_adj_per   = 1./KNa;
    % Normalize both forward and adjoint harmonics
    brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);
    EVa_per(:,1) = EVa_per(:,1)/sqrt(braket(EVa_per(:,1),EVa_per(:,1),brakopt));
    % Balance evaluation
    NB_phi_pert  = neutron_balance(fd_data,multig_data_per,(EVd_per(:,1)),k_n_per(1));
    NB_adj_pert  = neutron_balance(fd_data,multig_data_per,(EVa_per(:,1)),k_n_per(1));
    
    % Move to reactor directory
    try
        cd(cname)
    catch
        mkdir(cname)
    end
    save(['P_',num2str([1 2 3],'%g'),'_I_',num2str(pert_mat(4),'%g_'),ext,'_p.mat'],'Lp','Fp','Lp_adj','Fp_adj','EVd_per','EVa_per','KN_per','fd_data');
    cd ..
end

% Reference eigenvalues
k_inf_reg = zeros(N_reg,1);
for ii=1:N_reg
    k_inf_reg(ii) = MG_k_inf(NG,multig_data_unp,ii);
end
k_inf = max(k_inf_reg);



%% -------------------------------- Perturbation analysis ------------------------------------------------------------------------------------------------

% braket options
brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);

N_pert = N;
M_harm = nev;
SM_opt       = struct('N',N_pert,'M',M_harm);
ref_modes    = struct('EVd',EVd,'EVa',EVa,'F',F,'L',L,...
                      'Fp',Fp,'Lp',Lp,'fd_data',fd_data,'KN',KN);

[phi_pert,lambda,alpha] = GPT(ref_modes,SM_opt);

lambda_pert = 1/k_n(1)+sum(lambda,1);
flux_pert    = EVd(:,1)+sum(phi_pert,2);

fprintf('mu_ref-mu_GPT = %f [pcm] \r\n',1e5*(k_n_per(1)-1/lambda_pert))
%% -------------------------------- Diagnostic outputs ------------------------------------------------------------------------------------------------
% Neutron balance
phi = EVd(:,1);
adj = EVa(:,1);

% Call neutron balance evaluation (to check result accuracy and physical conservation of particles)
NB_phi = neutron_balance(fd_data,multig_data_per,phi,k_eff);

if verbosity>1
    fprintf('Numerical calculation: \r\n')
    fprintf('k_inf = %f \r\n',k_inf)
    fprintf('k_eff = %f \r\n',k_eff)
    fprintf('The difference (k_eff-k_inf) is %f [pcm] \r\n',(k_n(1)-k_inf)*1e5)
    if verbosity>5
        fprintf('-- NEUTRON BALANCE \r\n')
        for ii=1:N_reg
            fprintf('Neutron balance in %s layer: \r\n',layer_label_per{ii})
            
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
% direct perturbation
if flag_plot > 0
    harmonics    = 0;
    balance      = 0;
    balance_adj  = 0;
    dens_profile = 0;
    showplot     = struct('harmonics',harmonics,'dens_profile',dens_profile,'balance',balance,'balance_adj',balance_adj);

% % %     % direct perturbation
% % %     plot_opt     = struct('plot_save',plot_save,'plot_scale',plot_scale,'verbosity',verbosity,'NC',NC_type,'showplot',showplot);
% % %     fd_solution_pert = struct('EVd',EVd_per,'EVa',EVa_per,'k_n',k_n_per,'k_n_adj',k_n_adj_per,'k_eff',k_n_per(1),'k_inf',k_inf,'NG',NG,'NB_phi',NB_phi_pert,'NB_adj',NB_adj_pert);
% % %     display_solution(fd_data,fd_solution_pert,plot_opt,layer_label_per)
% % %     % unperturbed solution    
% % %     plot_opt        = struct('plot_save',plot_save,'plot_scale',plot_scale,'verbosity',verbosity,'NC',NC_type,'showplot',showplot);
% % %     fd_solution_unp = struct('EVd',EVd,'EVa',EVa,'k_n',k_n,'k_n_adj',k_n_adj,'k_eff',k_eff,'k_inf',k_inf,'NG',NG,'NB_phi',NB_phi,'NB_adj',NB_adj);
% % %     display_solution(fd_data,fd_solution_unp,plot_opt,layer_label_per)
    % GPT
    plot_opt        = struct('plot_save',plot_save,'plot_scale',plot_scale,'verbosity',verbosity,'NC',NC_type,'showplot',showplot);
    fd_solution_unp = struct('EVd',flux_pert,'EVa',EVd_per,'k_n',1/lambda_pert,'k_n_adj',1/lambda_pert,'k_eff',1/lambda_pert,'k_inf',k_inf,'NG',NG,'NB_phi',NB_phi,'NB_adj',NB_adj);
    display_solution(fd_data,fd_solution_unp,plot_opt,layer_label_per)
end

if mat_save>0
    if NG>0
        model_label = 'flux';
    else
        model_label = 'importance';
    end
    
    save([geometry,'_layers_',num2str(N_reg),'.mat'],'k_eff','fd_data','phi','model_label')
end

function fshape = define_shape(N_reg,t_reg)
if nargin>3
    shape = varargin{:}; % array of handles
else
    shape = num2str(ones(N_reg,1));
end
% f = cell(N_reg,1);
func = [];
for ii=1:N_reg
    if ii==1
        f = (sprintf('@(x,XS) (%s).*XS(%d).*(x>%f & x<=%f)',shape(ii),ii,t_reg(ii),t_reg(ii+1)));
    else
        f = (sprintf('+(%s).*XS(%d).*(x>%f & x<=%f)',shape(ii),ii,t_reg(ii),t_reg(ii+1)));
    end
    func = [func f]; % strings concatenation
end
fshape = eval(func);
end