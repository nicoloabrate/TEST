function SI = MG_SI(k,NG,multig_data,N_reg)
%% -- Multi-group approximation for k_infinity.
% This function takes in input the collapsed-energy material
% parameters and computes the Spectral Ratio (or Spectral Index, SI)
% assuming that the system is finite and homogeneous (hypothesis of
% Wigner's fundamental theorem of space-energy separability.

if abs(NG) ~= 2
    warning('This function works for a two-group model.')
end

if NG<0
    adjoint = 1;
    NG = -NG;
else
    adjoint = 0;
end

% Multi-group operators initialization
MGA = cell(NG,NG); % multi-group operator
ii  = N_reg;

for g=1:NG
        for gp=1:NG
            MGA{g,gp} = 0;
        end
        
end

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
        L = -g_data.XS_REM;
        % Diagonal blocks
        MGA{g,g} = MGA{g,g}+L;
        % Scattering matrices
        for gp=1:NG
            if gp~=g % skip diagonal terms
                % Non-diagonal blocks
                MGA{g,gp} =  MGA{g,gp}+g_data.XS_S0(:,gp);
            end
            
            %% Fission
            F = g_data.CHIT.*g_data.XS_NSF(:,gp);
            MGA{g,gp} = MGA{g,gp}+F/k;
        end
        
end
    
if adjoint == 1 % adjoint operator
    MGA = MGA';
end

MGA = cell2mat(MGA);
MGA = MGA./multig_data{1,ii}.DIFFCOEF';
[V,D] = eig(MGA);
D = diag(D);
[i,j]=find(D>0);
SI = V(1,i)/V(2,i);
end