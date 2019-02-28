function k_inf = MG_k_inf(NG,multig_data,N_reg)
%% -- Multi-group approximation for k_infinity.
% This function takes in input the collapsed-energy material
% parameters and computes the infinite multiplication factor for a
% homogeneous material.

if NG<0
    adjoint = 1;
    NG = -NG;
else
    adjoint = 0;
end

% Multi-group operators initialization
MGL = cell(NG,NG); % leakage
MGF = cell(NG,NG); % fission
ii  = N_reg;

for g=1:NG
    % Loading data according to group
    DFC      = multig_data{1,ii}.DIFFCOEF(g);
    XS_REM   = multig_data{1,ii}.XS_REM(g);
    XS_S0    = multig_data{1,ii}.XS_S0(g,:);
    XS_NSF   = multig_data{1,ii}.XS_NSF;
    CHIT     = multig_data{1,ii}.CHIT(g);
    
    XS_S0(:,g) = 0; % eliminating in-group scattering
    % Single group data provided to spatial discretization function
    g_data =  struct('DIFFCOEF',DFC,...
                     'XS_REM',XS_REM,...
                     'XS_S0',XS_S0,...
                     'XS_NSF',XS_NSF,...
                     'CHIT',CHIT);
    %% Leakage
        % Spatial discretization of diffusion operator
        L = g_data.XS_REM;
        % Diagonal blocks
        MGL{g,g} = L;
        % Scattering matrices
        for gp=1:NG
            if gp~=g % skip diagonal terms
                % Non-diagonal blocks
                MGL{g,gp} = -g_data.XS_S0(:,gp);
            end
            
            %% Fission
            F = g_data.CHIT.*g_data.XS_NSF(:,gp);
            MGF{g,gp} = F;
        end
        
    end
    


if adjoint == 1 % adjoint operator
    MGL = MGL';
    MGF = MGF';
end

MGL = cell2mat(MGL);
MGF = cell2mat(MGF);

k_N  = eig(MGL,MGF);
k_inf = 1/k_N(1);
end