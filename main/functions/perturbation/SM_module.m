function [phi,lambda] = SM_module(ref_modes,SM_opt)

% This script approximate the fundamental eigenpair of a  perturbed multi-layer slab with variable material properties
% using the Standard Method (A. Gandini, 1978)

%% -------------------------------- Perturbation analysis ------------------------------------------------------------------------------------------------
N              = SM_opt.N;               % perturbation order
M              = SM_opt.M;               % expansion order (number of reference system harmonics)
EVd            = ref_modes.EVd;
EVa            = ref_modes.EVa;
L              = ref_modes.L;
F              = ref_modes.F;
Lp             = ref_modes.Lp;
Fp             = ref_modes.Fp;
fd_data        = ref_modes.fd_data;

% array allocation
alpha          = zeros(M,N);     % alpha_{1,K} = 0 as normalization constant
lambda         = zeros(N,1);     % perturbed eigenvalue
Nt             = sum(fd_data.N,1);
NG             = length(EVd(:,1))/Nt;
phi            = zeros(Nt*NG,N); % perturbed system flux perturbations
brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);

% Check normalization to <phi+_m|F phi_m>
if norm(braket(EVa(:,1),F*EVd(:,1),brakopt)-1)>1e-4
    
    warning('Warning: Direct harmonics are not normalized properly.')
    % Harmonics ormalization
    C = zeros(size(EVd,2),1);
    
    for ii=1:size(EVd,2)
        C(ii)     = braket(EVa(:,ii),(F*EVd(:,ii)),brakopt);
        EVd(:,ii) = EVd(:,ii)/C(ii);
    end  
    
else
    
    C = ones(size(EVd,2),1);
    
end


mu = ref_modes.KN;  % reference eigenvalues need to be in this form
P0 = braket(1,F*EVd(:,1),brakopt);
P1 = braket(1,Fp*EVd(:,1),brakopt);
% Perturbed operators construction
dL = Lp-L;
dF = Fp-F;
% Perturbed eigenproblem definition
dA = dL-mu(1)*dF;

%% -- Perturbation chain

% 1st order perturbation
lambda(1) = braket(EVa(:,1),dA*EVd(:,1),brakopt)/C(1);

% Expansion coefficient evaluation
alpha(1,1) = P0/P1-1;
for m = 2:M % expansion order
    alpha(m,1) = -braket(EVa(:,m),dA*EVd(:,1),brakopt)/((mu(m)-mu(1))*C(m));
end

% flux perturbation evaluation
phi(:,1) = EVd(:,1:M)*alpha(:,1);


for n = 2:N % perturbation order
    
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
    lambda(n) = (braket(EVa(:,1),(dA*phi(:,n-1)),brakopt)+...
                -braket(EVa(:,1),S1_eig,brakopt)+...
                -braket(EVa(:,1),S2_eig,brakopt))/C(1);
    
    for m = 2:M % expansion order
        
        % Expansion coefficient evaluation
        for k = 1:n-1
            S1_alpha      = S1_alpha+lambda(k)*alpha(m,n-k+1);
            if n>2
                S2_alpha  = S2_alpha+lambda(k)*alpha(m,n-k);
            end
        end
        
        alpha(m,n) =  (-alpha(m,n-1)*braket(EVa(:,m),dA*EVd(:,m),brakopt)+...
                       S1_alpha*C(m)+...
                       S2_alpha*braket(EVa(:,m),dF*EVd(:,m),brakopt))/...
                      ((mu(m)-mu(1))*C(m));
                       
    end
    
    alpha(1,n) = (P0-P1-braket(1,Fp*(sum(phi(:,1:n-1),2)),brakopt)-braket(1,Fp*sum(alpha(2:end,n)*EVd(:,2:end),2),brakopt))/P1;
    phi(:,n) = EVd(:,1:M)*alpha(:,n);
    
end


end
