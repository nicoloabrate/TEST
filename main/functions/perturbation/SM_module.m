function [phi,lambda] = SM_module(ref_modes,SM_opt)

% This script approximate the fundamental eigenpair of a  perturbed multi-layer slab with variable material properties
% using the Standard Method (A. Gandini, 1978)

%% -------------------------------- Perturbation analysis ------------------------------------------------------------------------------------------------
N              = SM_opt.N;               % perturbation order
M              = SM_opt.M;               % expansion order (number of reference system harmonics)
varphi         = ref_modes.EVd;
varpsi         = ref_modes.EVa;
L              = ref_modes.L;
F              = ref_modes.F;
Lp             = ref_modes.Lp;
Fp             = ref_modes.Fp;
fd_data        = ref_modes.fd_data;

% array allocation
alpha          = zeros(M,N-1);     % alpha_{1,K} = 0 as normalization constant
lambda         = zeros(N-1,1);     % perturbed eigenvalue
Nt             = sum(fd_data.N,1);
NG             = length(varphi(:,1))/Nt;
phi            = zeros(Nt*NG,N-1); % perturbed system flux perturbations
brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);

% Check normalization to <phi+_m|F phi_m>
if norm(braket(varpsi(:,1),F*varphi(:,1),brakopt)-1)>1e-4
    fprintf('Direct harmonics are not normalized properly.\n')
    % Harmonics ormalization
    C = zeros(size(varphi,2),1);
    for ii=1:size(varphi,2)
        C(ii)     = braket(varpsi(:,ii),(F*varphi(:,ii)),brakopt);
        varphi(:,ii) = varphi(:,ii)/C(ii);
    end
    fprintf('Normalization completed. \n')
else
    C = ones(size(varphi,2),1);
end


mu = ref_modes.KN;  % reference eigenvalues need to be in this form
% Perturbed operators construction
dL = Lp-L;
dF = Fp-F;
% Perturbed eigenproblem definition
dA = dL-mu(1)*dF;

%% -- Perturbation chain
% 1st order perturbation
lambda(1) = braket(varpsi(:,1),dA*varphi(:,1),brakopt)/C(1);
% Expansion coefficient evaluation
for m = 2:M % expansion order
    alpha(m,1) = -braket(varpsi(:,m),dA*varphi(:,1),brakopt)/((mu(m)-mu(1))*C(m));
end
% flux perturbation evaluation
phi(:,1) = varphi(:,1:M)*alpha(:,1);

for n = 2:N-1 % perturbation order
    
    % sum coefficient initialization
    S1_eig   = 0;
    S2_eig   = 0;
    S1_alpha = 0;
    S2_alpha = 0;
    
    % sum over all previous contribution
    for k = 1:n-1
        S1_eig        = S1_eig+lambda(k)*F*phi(:,n-k+1);
        if n>2
            S2_eig    = S2_eig+lambda(k)*dF*phi(:,n-k);
        end
    end
    
    % Eigenvalue perturbation evaluation
    lambda(n) = (braket(varpsi(:,1),(dA*phi(:,n-1)),brakopt)+...
                -braket(varpsi(:,1),S1_eig,brakopt)+...
                -braket(varpsi(:,1),S2_eig,brakopt))/C(1);
    
    for m = 2:M % expansion order
        
        % Expansion coefficient evaluation
        for k = 1:n-1 % loop on over all previous coefficients
            
            S1_alpha      = S1_alpha+lambda(k)*alpha(m,n-k+1);
            
            if n>2
                
                for ii=1:M % loop over all harmonics contributions
                    S2_alpha  = S2_alpha+lambda(k)*braket(varpsi(:,m),dF*varphi(:,n-k),brakopt);
                end
                
            end            
        end
        
        alpha(m,n) =  (-braket(varpsi(:,m),dA*phi(:,n-1),brakopt)+...
                       S1_alpha*C(m)+...
                       S2_alpha)/...
                      ((mu(m)-mu(1))*C(m));
        
    end
    
    phi(:,n) = varphi(:,1:M)*alpha(:,n);
    
end


end
