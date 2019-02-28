function A = FD(g_data,fd_data)
%  This function discretizes second order derivatives for
%  1D cartesian, spherical and cylindrical coordinate systems.
%  It takes in input:  g_data  (struct) --> Material data (diffusion
%                                           coefficient, removal XS).
%                      fd_data (struct) --> Mesh data (grid, number of
%                                           nodes, step-size, material data
%                                           shape.

%% Sub-matrix: leakage, absorption, out-scattering

% Gathering data for discretization
% Material data
DFC     = g_data.DIFFCOEF;
XS_REM  = g_data.XS_REM;
% Mesh data
dx      = fd_data.dx;
N       = fd_data.N;
FD_grid = fd_data.FD_grid;
N_reg   = fd_data.N_reg;
f       = fd_data.f;

% Initialization
N_old = 0;
Nt    = sum(N,1);
xp    = zeros(Nt,1);
dxv   =  zeros(Nt,1);

% Geometrical factor 
s = 0*(fd_data.geom==1)+1*(fd_data.geom==3)+2*(fd_data.geom==2); % useful to keep a unique definition of discretized operator

% Definition of staggered (ghost) grid where material data are interpolated
for ii = 1:length(dx)
    xp(1+(ii>1)*(N_old-1):N(ii)+(ii>1)*((N_old-1))) = FD_grid(1+(ii>1)*(N_old-1):N(ii)+(ii>1)*((N_old-1)))+dx(ii)/2;
    dxv(1+(ii>1)*(N_old):N(ii)+(ii>1)*((N_old))) = dx(ii);
    N_old = N_old+N(ii);
end
xp(end) = xp(end-1)+dx(end)/2;

% ======================================= Diagonals construction ==========================================================================================================
% allocation
md = zeros(Nt,1);
ud = md;
ld = md;
N_old = 0;
for ii=1:N_reg
    
    N_int = N_old+N(ii);
    p = (N_old-1)*(ii>1); % p and q are indeces used to have a less heavy notation when loading arrays
    q = N_old*(ii>1);
    
    % Lower diagonal definition
    ld(1+p:N(ii)-2+p+1)   = -(DFC(ii)./f{ii}(xp(1+p:N(ii)-2+p+1)))./dxv(1+q:N(ii)-2+q+1).^2+...
                             (s.*(DFC(ii)./f{ii}(FD_grid(2+p:N(ii)-1+p+1)))'./(2*dxv(1+q:N(ii)-2+q+1)'.*FD_grid(2+p:N(ii)-1+p+1))');
    % Main diagonal definition
    md(2+p:N(ii)-1+p+1)   =  (DFC(ii)./f{ii}(xp(1+p:N(ii)-2+p+1)))./dxv(1+q:N(ii)-2+q+1).^2+...
                             (DFC(ii)./f{ii}(xp(2+p:N(ii)-1+p+1)))./dxv(2+q:N(ii)-1+q+1).^2+...
                             XS_REM(ii).*f{ii}(FD_grid(2+p:N(ii)-1+p+1))'; 
    % Upper diagonal definition                     
    ud(3+p:N(ii)+p+1)     = -(DFC(ii)./f{ii}(xp(2+p:N(ii)-1+p+1)))./dxv(2+q:N(ii)-1+q+1).^2-...
                             (s.*(DFC(ii)./f{ii}(FD_grid(2+p:N(ii)-1+p+1)))'./(2*dxv(1+q:N(ii)-2+q+1)'.*FD_grid(2+p:N(ii)-1+p+1))');  
                         
    % Interface                     
    if N_reg>1 && ii<N_reg % check interface is not boundary
        % weighting material coefficient around interface (the function lightens the script)
        D   = weight_avg(DFC(ii)./f{ii}(FD_grid(N_int-1)),DFC(ii+1)./f{ii+1}(FD_grid(N_int+1)),  dxv(N_int-1),dxv(N_int+1));
        XSR = weight_avg(XS_REM(ii).*f{ii}(FD_grid(N_int-1)),XS_REM(ii+1).*f{ii+1}(FD_grid(N_int+1)),dxv(N_int-1),dxv(N_int+1));

        % Main diagonal
        md(N_int)   =  (DFC(ii)/f{ii}(xp(N_int-1))/dxv(N_int)+DFC(ii+1)/f{ii+1}(xp(N_int))/dxv(N_int+1))/(dxv(N_int)/2+dxv(N_int+1)/2)+...
                       XSR;
        % Lower diagonal
        ld(N_int-1) =  (-DFC(ii)/f{ii}(xp(N_int-1))/dxv(N_int)) / (dxv(N_int)/2+dxv(N_int+1)/2)+...
                       s*D/(2*(FD_grid(N_int+1)-FD_grid(N_int-1))*FD_grid(N_int));
        % Upper diagonal
        ud(N_int+1) =  (-DFC(ii+1)/f{ii+1}(xp(N_int))/dxv(N_int+1)) / (dxv(N_int)/2+dxv(N_int+1)/2)-...
                       s*D/(2*(FD_grid(N_int+1)-FD_grid(N_int-1))*FD_grid(N_int));
        
    end
    N_old = N_old+N(ii);
    if length(ud)>Nt
        ud(end) = [];
    end
end

% ======================================= Setting BCs ==========================================================================================================

% Left boundary
if fd_data.geom==1 % Cartesian geometry
    % 1 for vacuum, 2 for reflective BC
    md(1) = 1;  
    ud(2) = 0*(fd_data.BC(1)==1)-1*(fd_data.BC(1)==2);
else % Spherical or cylindrical geometry
    % Finiteness BC
    md(1) = 1;
    ud(2) = -1; 
end
% Right boundary, 1 for vacuum, 2 for reflective BC
md(Nt) = 1;
ld(Nt-1) = 0*(fd_data.BC(2)==1)-1*(fd_data.BC(2)==2);

% ======================================= Matrix construction ==========================================================================================================
A = spdiags([ld,md,ud],-1:1,Nt,Nt);
end

function C_avg = weight_avg(C1,C2,dr1,dr2)
C_avg = (C1*dr1/2+C2*dr2/2)/(dr1/2+dr2/2);
end