function NB = neutron_balance(fd_data,multig_data,phi,k_eff)
% This function takes in input the core data, the flux/adjoint and the effective multiplication
% factor to evaluate the neutron balance terms.

% Define spatial mesh features and number of groups
N_reg   = fd_data.N_reg;
N       = fd_data.N;
dx      = fd_data.dx;
FD_grid = fd_data.FD_grid;
ff      = fd_data.f;
Nt      = sum(N,1);
NG      = length(phi)/Nt;
geom    = fd_data.geom;
N_old   = 0;

% Store flux in group-wise matrix
phi_g = zeros(Nt,NG);
for g = 1:NG
    phi_g(:,g) = phi(1+(g-1)*Nt:Nt+(g-1)*Nt);
end

% Variable allocation
FISS        = zeros(N_reg,NG);
ABS_CAPT    = zeros(N_reg,NG);
LEAK_L      = zeros(N_reg,NG);
LEAK_R      = zeros(N_reg,NG);
NG_IMBAL    = zeros(N_reg,NG);
ABS_FISS    = zeros(N_reg,NG);
IN_SCATT    = zeros(N_reg,NG);
OUT_SCATT   = zeros(N_reg,NG);
ALB_R       = zeros(N_reg,NG);
ALB_L       = zeros(N_reg,NG);
J_inward_L  = zeros(N_reg,NG);
J_outward_L = zeros(N_reg,NG);
J_inward_R  = zeros(N_reg,NG);
J_outward_R = zeros(N_reg,NG);

% Integration weight
switch geom
    case 1
        wi = ones(1,Nt);       % line integral
    case 2
        wi = 4*pi*FD_grid.^2;  % volume integral
    case 3
        wi = 2*pi*FD_grid;     % surface integral
end



% Check if user wants neutron or importance balance
if k_eff>0

    % Loops on regions and energy groups

    for ii=1:N_reg % span all regions
        % Define absorption, capture and fission terms (to avoid statistical
        % difference in Serpent-2 output)
        multig_data{ii}.XS_ABS  = multig_data{ii}.XS_REM-sum(multig_data{ii}.XS_S0,1)+diag(multig_data{ii}.XS_S0)';
        multig_data{ii}.XS_CAPT = multig_data{ii}.XS_ABS-multig_data{ii}.XS_FISS;
        x = 1+N_old:N_old+N(ii);
        
        for g = 1:NG % span all groups
            
            for gp=1:NG
                if gp~=g % skip diagonal terms
                    IN_SCATT(ii,g)  = IN_SCATT(ii,g)+sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_S0(g,gp).*ff{ii}(FD_grid(x)).*phi_g(x,gp)');
                    OUT_SCATT(ii,g) = OUT_SCATT(ii,g)+sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_S0(gp,g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
                end
                % Fission
                FISS(ii,g)  = FISS(ii,g)+1/k_eff*sum(dx(ii)*wi(x).*multig_data{1,ii}.CHIT(g)*multig_data{1,ii}.XS_NSF(gp).*ff{ii}(FD_grid(x)).*phi_g(x,gp)');
            end
            
            ABS_FISS(ii,g)     = sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_FISS(g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
            ABS_CAPT(ii,g)     = sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_CAPT(g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
            % Net leakages
            LEAK_R(ii,g)       = wi(N_old+N(ii))*(multig_data{1,ii}.DIFFCOEF(g)./ff{ii}(FD_grid(N_old+N(ii))).*(phi_g(N(ii)+N_old-1,g)-phi_g(N(ii)+N_old,g))/dx(ii));
            LEAK_L(ii,g)       = wi(1+N_old)*(multig_data{1,ii}.DIFFCOEF(g)./ff{ii}(FD_grid(1+N_old)).*(phi_g(1+N_old+1,g)-phi_g(1+N_old,g))/dx(ii));
            % Group-wise local imbalance
            NG_IMBAL(ii,g)     = FISS(ii,g)+IN_SCATT(ii,g)-OUT_SCATT(ii,g)-ABS_FISS(ii,g)-ABS_CAPT(ii,g)-LEAK_L(ii,g)-LEAK_R(ii,g);
            % Partial currents
            J_inward_L(ii,g)   = (1/4*phi_g(1+N_old,g)-1/2*(multig_data{1,ii}.DIFFCOEF(g)*(phi_g(1+N_old+1,g)-phi_g(1+N_old,g))/dx(ii)));
            J_inward_R(ii,g)   = (1/4*phi_g(N(ii)+N_old,g)-1/2*(multig_data{1,ii}.DIFFCOEF(g)*(phi_g(N(ii)+N_old,g)-phi_g(N(ii)+N_old-1,g))/dx(ii)));
            J_outward_L(ii,g)  = (1/4*phi_g(1+N_old,g)+1/2*(multig_data{1,ii}.DIFFCOEF(g)*(phi_g(1+N_old+1,g)-phi_g(1+N_old,g))/dx(ii)));
            J_outward_R(ii,g)  = (1/4*phi_g(N(ii)+N_old,g)+1/2*(multig_data{1,ii}.DIFFCOEF(g)*(phi_g(N(ii)+N_old,g)-phi_g(N(ii)+N_old-1,g))/dx(ii)));
        end
        
        N_old = N_old+N(ii); % increment to move along the core
        
    end

    % Normalization
    TOT_FISS  = sum(sum(FISS)); % used as normalization constant
    FISS      = FISS/TOT_FISS;
    IN_SCATT  = IN_SCATT/TOT_FISS;
    OUT_SCATT = OUT_SCATT/TOT_FISS;
    ABS_FISS  = ABS_FISS/TOT_FISS;
    ABS_CAPT  = ABS_CAPT/TOT_FISS;
    LEAK_L    = LEAK_L/TOT_FISS;
    LEAK_R    = LEAK_R/TOT_FISS;
    NG_IMBAL  = NG_IMBAL/TOT_FISS;

    % Total balance

    % Production
    PROD = sum(sum(FISS))+sum(sum(IN_SCATT));
    % Losses
    LOSS      = sum(sum(ABS_FISS))+sum(sum(ABS_CAPT))+sum(sum(OUT_SCATT))+sum(sum(LEAK_R))+sum(sum(LEAK_L));
    TOT_IMBAL = PROD-LOSS;
    K_BAL     = PROD/LOSS; % ratio production to losses (should be almost 1)

else
    k_eff = -k_eff;
    % Loops on regions and energy groups

    for ii=1:N_reg % span all regions
        % Define absorption, capture and fission terms (to avoid statistical
        % difference in Serpent-2 output)
        multig_data{ii}.XS_ABS  = multig_data{ii}.XS_REM-sum(multig_data{ii}.XS_S0,1)+diag(multig_data{ii}.XS_S0)';
        multig_data{ii}.XS_CAPT = multig_data{ii}.XS_ABS-multig_data{ii}.XS_FISS;
        x = 1+N_old:N_old+N(ii);
        
        for g = 1:NG % span all groups
            
            for gp=1:NG
                
                if gp~=g % skip diagonal terms
                    IN_SCATT(ii,g)  = IN_SCATT(ii,g)+sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_S0(gp,g).*ff{ii}(FD_grid(x)).*phi_g(x,gp)');
                    OUT_SCATT(ii,g) = OUT_SCATT(ii,g)+sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_S0(gp,g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
                end
                
                % Fission
                FISS(ii,g)  = FISS(ii,g)+1/k_eff*sum(dx(ii)*wi(x).*multig_data{1,ii}.CHIT(gp)*multig_data{1,ii}.XS_NSF(g).*ff{ii}(FD_grid(x)).*phi_g(x,gp)');
            end
            
            ABS_FISS(ii,g)  = sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_FISS(g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
            ABS_CAPT(ii,g)  = sum(dx(ii)*wi(x).*multig_data{1,ii}.XS_CAPT(g).*ff{ii}(FD_grid(x)).*phi_g(x,g)');
            % Net leakages
            LEAK_R(ii,g)    = wi(N_old+N(ii))*(multig_data{1,ii}.DIFFCOEF(g)*ff{ii}(FD_grid(N_old+N(ii)))*(phi_g(N(ii)+N_old-1,g)-phi_g(N(ii)+N_old,g))/dx(ii));
            LEAK_L(ii,g)    = wi(1+N_old)*(multig_data{1,ii}.DIFFCOEF(g)*ff{ii}(FD_grid(1+N_old))*(phi_g(1+N_old+1,g)-phi_g(1+N_old,g))/dx(ii));
            % Group-wise local imbalance
            NG_IMBAL(ii,g)  = FISS(ii,g)+IN_SCATT(ii,g)-OUT_SCATT(ii,g)-ABS_FISS(ii,g)-ABS_CAPT(ii,g)-LEAK_L(ii,g)-LEAK_R(ii,g);
        end
        
        N_old = N_old+N(ii); % increment to move along the core
        
    end

    % Normalization
    TOT_FISS  = sum(sum(FISS)); % used as normalization constant
    FISS      = FISS/TOT_FISS;
    IN_SCATT  = IN_SCATT/TOT_FISS;
    OUT_SCATT = OUT_SCATT/TOT_FISS;
    ABS_FISS  = ABS_FISS/TOT_FISS;
    ABS_CAPT  = ABS_CAPT/TOT_FISS;
    LEAK_L    = LEAK_L/TOT_FISS;
    LEAK_R    = LEAK_R/TOT_FISS;
    NG_IMBAL  = NG_IMBAL/TOT_FISS;

    % Total balance

    % Production
    PROD = sum(sum(FISS))+sum(sum(IN_SCATT));
    % Losses
    LOSS      = sum(sum(ABS_FISS))+sum(sum(ABS_CAPT))+sum(sum(OUT_SCATT))+sum(sum(LEAK_R))+sum(sum(LEAK_L));
    TOT_IMBAL = PROD-LOSS;
    K_BAL     = PROD/LOSS; % ratio production to losses (should be almost 1)


end
% Gathering data in structure
NB = struct('FISS',FISS,...
            'ABS_FISS',ABS_FISS,...
            'ABS_CAPT',ABS_CAPT,...
            'IN_SCATT',IN_SCATT,...
            'OUT_SCATT',OUT_SCATT,...
            'LEAK_L',LEAK_L,...
            'LEAK_R',LEAK_R,...
            'NG_IMBAL',NG_IMBAL,...
            'TOT_IMBAL',TOT_IMBAL,...
            'K_BAL',K_BAL,...
            'TOT_FISS',TOT_FISS,...
            'ALB_R',ALB_R,...
            'ALB_L',ALB_L,...
            'J_inward_L',J_inward_L,...
            'J_inward_R',J_inward_R,...
            'J_outward_L',J_outward_L,...
            'J_outward_R',J_outward_R);
        
if TOT_IMBAL > 0.01
   fprintf('WARNING: Total neutron imbalance >1%% \r\n') 
end
end
