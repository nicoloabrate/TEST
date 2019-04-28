function [phi,lambda,alpha] = GPT(ref_modes,SM_opt)

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
alpha          = NaN*zeros(M,N-1);     % alpha_{1,K} = 0 as normalization constant
lambda         = NaN*zeros(N-1,1);     % perturbed eigenvalue
Nt             = sum(fd_data.N,1);
NG             = length(varphi(:,1))/Nt;
phi            = NaN*zeros(Nt*NG,N-1); % perturbed system flux perturbations
brakopt        = struct('FD_grid',fd_data.FD_grid,'NG',NG,'Nt',Nt);

% Check normalization to <phi+_m|F phi_m>
% if norm(braket(varpsi(:,1),F*varphi(:,1),brakopt)-1)>1e-4
%     fprintf('Direct harmonics are not normalized properly.\n')
%     % Harmonics normalization
    C = zeros(size(varphi,2),1);
    for ii=1:size(varphi,2)
        C(ii)     = braket(varpsi(:,ii),(F*varphi(:,ii)),brakopt);
        %varphi(:,ii) = varphi(:,ii)/C(ii);
    end
%     fprintf('Normalization completed. \n')
% else
%     C = ones(size(varphi,2),1);
% end

mu = ref_modes.KN;  % reference eigenvalues need to be in this form
% Perturbed operators construction
dL = Lp-L;
dF = Fp-F;
% Perturbed eigenproblem definition
dA = dL-mu(1)*dF;

% 0th order perturbed power P'
Pp = braket(1,Fp*varphi(:,1),brakopt);

%% -- Perturbation chain
for n = 1:N % perturbation order
    
    % - eigenvalue perturbations
    
    if n == 1
        S = dA*varphi(:,1); % 0th order perturbation
    else
        S = dA*phi(:,n-1);
    end
    % sum coefficients initialization
    S1_eig   = 0;
    S2_eig   = 0;
    % sum over all previous contribution
    for k = 1:n-1
        if (n-k)>0 % check on (n-k)th order
            S1_eig        = S1_eig+lambda(k)*F*phi(:,n-k);
            if (n-k-1)>0 % check on (n-k-1)th order
                S2_eig    = S2_eig+lambda(k)*dF*phi(:,n-k-1);
            else % evaluate 0th order perturbation (i.e. unperturbed flux)
                S2_eig    = S2_eig+lambda(k)*dF*varphi(:,1);
            end
        end
    end
    
    % nth order perturbation evaluation
    lambda(n) = (braket(varpsi(:,1),S,brakopt)+...
                -braket(varpsi(:,1),S1_eig,brakopt)+...
                -braket(varpsi(:,1),S2_eig,brakopt))/C(1);
    
    % - eigenvector perturbations    
    for m = 2:M % loop over all expansion coefficients (m=1 is evaluated afterwards)
        
        % sum coefficients initialization
        S1_alpha = 0;
        S2_alpha = 0;
        for k = 1:n-1 % sum over all previous perturbation order

            if (n-k)>0 % check on (n-k)th order
                
                S1_alpha = S1_alpha+lambda(k)*alpha(m,n-k);
                
                if (n-k-1)>0 % check on (n-k-1)th order     
                    S2_alpha  = S2_alpha+lambda(k)*dF*phi(:,n-k-1);        
                else % evaluate 0th order perturbation (i.e. unperturbed flux)                    
                    S2_alpha  = S2_alpha+lambda(k)*dF*varphi(:,1);
                end
                
            end
            
        end
        % coefficient evaluation
        alpha(m,n) =  (-braket(varpsi(:,m),S,brakopt)+... % S already evaluated for lambda perturbations 
                      S1_alpha*C(m)+...
                      braket(varpsi(:,m),S2_alpha,brakopt))/...
                      ((mu(m)-mu(1))*C(m));
       
    end
    
    % first order power is equal to the total power
    if n == 1
        alpha(1,n) = -(braket(1,dF*varphi(:,1),brakopt)+braket(1,Fp*(varphi(:,2:M)*alpha(2:M,n)),brakopt))/Pp;
    else
        alpha(1,n) = -braket(1,Fp*(varphi(:,2:M)*alpha(2:M,n)),brakopt)/Pp;
    end
    % perturbation reconstruction
    phi(:,n) = varphi(:,1:M)*alpha(:,n);
    
end


end
