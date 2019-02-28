function [MGL,MGF] = MG(NG,PN,multig_data,fd_data)
%% -- Energy groups discretization model.
% This function takes in input the collapsed-energy material
% parameters and builds the leakage and fission operators in
% order to solve the criticality problem.
% It can be applied both to transport and diffusion equations.
% The first is approximated with PN approximation, the second
% can be seen as a P1 transport equation.
% This module calls also a finite difference scheme for spatial
% approximation.
% Provide a negative number of groups to solve adjoint equation.

if NG<0
    adjoint = 1;
    NG = -NG;
else
    adjoint = 0;
end

% Retrieving spatial grid data
N          = fd_data.N;
FD_grid    = fd_data.FD_grid;
N_reg      = fd_data.N_reg; 
f          = fd_data.f;
dx         = fd_data.dx;
Nt         = sum(N,1);

% Macro-group constants allocation
DFC    = zeros(N_reg,1); 
XS_REM = zeros(N_reg,1); 
CHIT   = zeros(N_reg,1); 
XS_S0  = zeros(N_reg,NG); 
XS_NSF = zeros(N_reg,NG); 

% Multi-group operators allocation
MGL = cell(NG,NG); % leakage
MGF = cell(NG,NG); % fission

for g=1:NG
    
    % Loading group-wise data
    for ii = 1:N_reg
        DFC(ii)      = multig_data{1,ii}.DIFFCOEF(g);
        XS_REM(ii)   = multig_data{1,ii}.XS_REM(g);
        XS_S0(ii,:)  = multig_data{1,ii}.XS_S0(g,:);
        XS_NSF(ii,:) = multig_data{1,ii}.XS_NSF;
        CHIT(ii)     = multig_data{1,ii}.CHIT(g);
    end
    
    XS_S0(:,g) = 0; % eliminating in-group scattering
    
    % Single group data provided to spatial discretization function
    g_data =  struct('DIFFCOEF',DFC,...
                     'XS_REM',XS_REM,...
                     'XS_S0',XS_S0,...
                     'XS_NSF',XS_NSF,...
                     'CHIT',CHIT);
    
    if PN == 1 % Diffusion approximation
        %% Leakage
        % Spatial discretization of diffusion operator
        L = FD(g_data,fd_data);
        % Diagonal blocks
        MGL{g,g} = L;
        % Scattering matrices
        for gp=1:NG
            
            if gp~=g % skip diagonal terms
                
                N_int = 0;
                N_old = 0;
                md    = zeros(Nt,1); % Main diagonal definition, incorporating BCs
                
                for ii=1:N_reg
                    
                    p = (N_old-1)*(ii>1); % p and q are indeces used to have a less heavy notation when loading arrays
                    md(2+p:N(ii)+p)   = -g_data.XS_S0(ii,gp)*f{ii}(FD_grid(2+p:N(ii)+p));
                    
                    % Interface condition to impose flux continuity
                    if N_reg>1 && ii<N_reg
                        N_int     = N_int+N(ii);
                        md(N_int) = -weight_avg(g_data.XS_S0(ii,gp)*f{ii}(FD_grid(N_int-1)),g_data.XS_S0(ii+1,gp)*f{ii+1}(FD_grid(N_int+1)),dx(ii),dx(ii+1));
                    end
                    N_old = N_old+N(ii);
                end
                md(end) = 0; % BCs
                % Building matrix
                L = spdiags(md,0,Nt,Nt);
                % Non-diagonal blocks
                MGL{g,gp} = L;
            end % for N_reg (Scattering operator)
            
            N_int = 0;
            N_old  = 0;
            md = zeros(Nt,1); % Main diagonal definition, incorporating BCs
            
            for ii=1:N_reg
                %% Fission
                p = (N_old-1)*(ii>1); % p and q are indeces used to have a less heavy notation when loading arrays
                md(2+p:N(ii)+p)   = g_data.CHIT(ii).*g_data.XS_NSF(ii,gp)*f{ii}(FD_grid(2+p:N(ii)+p));
                % Interface condition to impose flux continuity
                if N_reg>1 && ii<N_reg
                    N_int = N_int+N(ii);
                    md(N_int) = weight_avg(g_data.CHIT(ii).*g_data.XS_NSF(ii,gp)*f{ii}(FD_grid(N_int)),g_data.CHIT(ii+1).*g_data.XS_NSF(ii+1,gp)*f{ii+1}(FD_grid(N_int+1)),dx(ii),dx(ii+1));
                end
                N_old = N_old+N(ii);
            end % for N_reg (Fission operator)
            md(end) = 0; % BCs
            F = spdiags(md,0,Nt,Nt);
            MGF{g,gp} = F;
            
        end % for gp
        
    end % if PN==1
    
end % for g

if adjoint == 1 % adjoint operator
    MGL = MGL';
    MGF = MGF';
end

MGL = cell2mat(MGL);
MGF = cell2mat(MGF);

end

% This function evaluate the weighted average of the material properties
% around the interface.
function C_avg = weight_avg(C1,C2,dr1,dr2)
C_avg = (C1*dr1/2+C2*dr2/2)/(dr1/2+dr2/2);
end