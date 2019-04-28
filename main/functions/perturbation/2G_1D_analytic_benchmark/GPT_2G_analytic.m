clear all
% close all
clc


% define perturbation order
M = 20;
% define expansion order (i.e. number of harmonics)
N = 50;
% define number of spatial points for the evaluation
NP = 5000;
NG = 2;
% reading nuclear data from library
T = 300;
mat_id  = 1;
NG = 2;
mat_name = 'AGPT';
ND = readnucdata(mat_name,NG,mat_id,T);
DL = sqrt(ND.DIFFCOEF./ND.XS_REM);

% adding perturbations
% pert_mat   = [1 2 0 1;   % fission XS, in-group, out-group, pert. amplitude
%               2 1 0 1;   % capture XS, in-group, out-group, pert. amplitude
%               2 2 0 1;   % capture XS, in-group, out-group, pert. amplitude
%               3 1 0 1;   % nubar, in-group, out-group, pert. amplitude
%               3 2 0 1;   % nubar, in-group, out-group, pert. amplitude
%               4 1 2 1;]; % scattering XS, in-group, out-group, pert. amplitude
%         

pert_mat   = [4 2 1 100]; 

% define geometric parameters
H  = 150; % [cm] domain thickness
dx = 5;  % [cm] perturbation width
x0 = 50;  % [cm] perturbation center coordinate
x  = linspace(0,H,NP);

PND = add_perturbation(pert_mat,ND); % perturbed nuclear data
kinf  = ND.XS_NSF(2)*ND.XS_S0(2,1)/(ND.XS_REM(2)*ND.XS_REM(1));

% explicit perturbations (analytic case)
dXSR1  = PND.XS_REM(1)-ND.XS_REM(1);
dXSR2  = PND.XS_REM(2)-ND.XS_REM(2);
dXSM12 = PND.XS_S0(2,1)-ND.XS_S0(2,1);
dD1    = PND.DIFFCOEF(1)-ND.DIFFCOEF(1);
dD2    = PND.DIFFCOEF(2)-ND.DIFFCOEF(2);
dNUBAR = PND.NUBAR(2)-ND.NUBAR(2);
dXSF2  = PND.XS_FISS(2)-ND.XS_FISS(2);

% --- define analytic expressions
Bn   = @(n) pi./H*(n);  % buckling
mu   = @(n) (1+DL(1)^2*Bn(n).^2).*(1+DL(2)^2*Bn(n).^2)/kinf; % eigenvalue sequence
psi  = @(n) ND.XS_S0(2,1)./(ND.DIFFCOEF(2)*(Bn(n).^2+1/DL(2)^2)); % forward spectral index
psiX = @(n) ND.XS_S0(2,1)./(ND.DIFFCOEF(1)*(Bn(n).^2+1/DL(1)^2)); % adjoint spectral index
varphi  = @(x,n) sin(Bn(n)*x); % function shape
AN  = @(n) sqrt(2./(H*(1+psi(n).^2))); % forward eigenfunctions normalization constants (they depend on the spectral index of the harmonics)
ANX = @(n) sqrt(2./(H*(1+psiX(n).^2))); % adjoint eigenfunctions normalization constants (they depend on the spectral index of the harmonics)

% inner product definitions
BK1 = @(a,b,ii) AN(ii).*ANX(ii).*psi(ii).*psiX(ii).*ND.NUBAR(2).*ND.XS_FISS(2).*I1(a,b,ii,ii,Bn); % braket 1 (denominator)

NC  = @(ii) AN(ii).*ANX(ii).*psi(ii).*psiX(ii).*H/2.*ND.NUBAR(2).*ND.XS_FISS(2); % Normalization constant (redundant with I1 expression to check consistency)

BK2 = @(a,b,ii,jj) (dXSR1.*psiX(ii)+dXSR2.*psi(jj)-dXSM12-...       % scattering, absorption
                    mu(1).*(dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+...    % fission
                    dNUBAR.*ND.XS_FISS(2)).*psiX(ii).*psi(jj)).*... % fission
                    ANX(ii).*AN(jj).*I1(a,b,ii,jj,Bn)+...
                    (ANX(ii).*psiX(ii).*AN(jj).*dD1+...             % fast diffusion
                    psi(jj).*AN(jj).*ANX(ii).*dD2).*...             % thermal diffusion
                    Bn(ii).*Bn(jj).*I3(a,b,ii,jj,Bn);
                
BK3 = @(a,b,ii,jj) AN(jj).*psi(jj).*ANX(ii).*psiX(ii).*(dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+dNUBAR.*ND.XS_FISS(2)).*I1(a,b,ii,jj,Bn);
BK4 = @(a,b,ii) AN(ii).*psi(ii).*I2(a,b,ii,Bn);
% reference solutions
PHI  = @(x,ii) [AN(ii)*varphi(x,ii)  AN(ii)*psi(ii)*varphi(x,ii)];
PHIX = @(x,ii) [ANX(ii)*psiX(ii)*varphi(x,ii) ANX(ii)*varphi(x,ii)];
EV  = NaN*ones(NG*NP,N);
EVX = NaN*ones(NG*NP,N);
for ii=1:N
    EV(:,ii)  = PHI(x,ii);
    EVX(:,ii) = PHIX(x,ii);
end
keff  = kinf/((1+DL(1)^2*Bn(1)^2)*(1+DL(2)^2*Bn(1)^2));

% Loop on perturbation order
a = NaN*ones(M,N);      % allocate the set of expansion coefficients
lambda = NaN*ones(M,1); % allocate eigenvalue perturbation array
phi = NaN*ones(NG*NP,M); % allocate eigenvalue perturbation array

for m = 1:M
    % lambda evaluation
    if m == 1
        lambda(m) = BK2(x0-dx/2,x0+dx/2,1,1)/BK1(0,H,1);
    else
        % re-initialize the sums
        S1 = 0; % sum over <F> term
        S2 = 0; % sum over <dF> term
        for kk = 1:m-1 % Loop on previous perturbation order
            p1 = m-kk;
            S1  =S1+lambda(kk)*a(p1,1)*BK1(0,H,1);
            p2 = m-kk-1;
            if p2>0
                S2 = S2+lambda(kk)*a(p2,1:N)*BK3(x0-dx/2,x0+dx/2,1,1:N)';
            else % p2 coincides with the 0-th order term, i.e. the reference harmonic
                S2 = S2+lambda(kk)*BK3(x0-dx/2,x0+dx/2,1,1)';
            end
        end
        lambda(m) = (a(m-1,1:N)*BK2(x0-dx/2,x0+dx/2,1,1:N)'-S1-S2)/BK1(0,H,1);
    end
    
    % exp. coefficient evaluation
    if m == 1
        a(m,2:N) = -BK2(x0-dx/2,x0+dx/2,2:N,1)./((mu(2:N)-mu(1)).*BK1(0,H,2:N));
        % arbitrary constant (1st order perturbed power equal to reference
        % system power)
        a(m,1)   = -(((dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+dNUBAR.*ND.XS_FISS(2))*BK4(x0-dx/2,x0+dx/2,1))+...
                    (a(m,2:N)*(ND.XS_NSF(2)*BK4(0,H,2:N)'+(dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+dNUBAR.*ND.XS_FISS(2))*BK4(x0-dx/2,x0+dx/2,2:N)')))./...
                   (ND.XS_NSF(2)*BK4(0,H,1)+(dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+dNUBAR.*ND.XS_FISS(2))*BK4(x0-dx/2,x0+dx/2,1));
    else
        for n = 2:N % loop on expansion order
            % re-initialize the sums
            S1 = 0; % sum over <F> term
            S2 = 0; % sum over <dF> term
            for kk = 1:m-1 % Loop on previous perturbation order
                p1 = m-kk;
                S1=S1+lambda(kk)*a(p1,n)*BK1(0,H,n);
                p2 = m-kk-1;
                if p2>0
                    S2 = S2+lambda(kk)*a(m-kk-1,1:N)*BK3(x0-dx/2,x0+dx/2,n,1:N)';
                else
                    S2 = S2+lambda(kk)*BK3(x0-dx/2,x0+dx/2,n,1)';
                end
            end
            a(m,n) = (-a(m-1,1:N)*BK2(x0-dx/2,x0+dx/2,n,1:N)'+S1+S2)./((mu(n)-mu(1)).*BK1(0,H,n));
        end
        % coefficients for m=0 (fixing the 1st order, the others require
        % that the power perturbations go to zero.
        a(m,1)   = -a(m,2:N)*(ND.XS_NSF(2)*BK4(0,H,2:N)+(dNUBAR.*dXSF2+ND.NUBAR(2).*dXSF2+dNUBAR.*ND.XS_FISS(2))*BK4(x0-dx/2,x0+dx/2,1:N-1))'./(ND.XS_NSF(2)*BK4(0,H,1)+dNUBAR*dXSF2*BK4(x0-dx/2,x0+dx/2,1));
    end
end

% perturbed eigenvalue
lambda_pert = mu(1)+sum(lambda);

% higher-order harmonics
PHIN  = @(x,n) [AN(n)*varphi(x,n); AN(n)*psi(n)*varphi(x,n)];
PHIP  = 0;
for jj=1:M
    for ii=1:N
        PHIP = PHIP+a(jj,ii)*PHIN(x,ii);
    end
end
flux_pert = PHIN(x,1)+PHIP;

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',22)
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','latex');

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'Color', 'w');
hold on
box on
plot(x,flux_pert(1,:),'b')
plot(x,flux_pert(2,:),'r')
grid on