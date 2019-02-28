function fd_data = space_grid(mesh_opt,layer_struct,multig_data,f_shape)
% Sample points to compute system diffusion length and discretizes taking NP points per each 
% characteristic length.

% retrieving data from structures
% mesh options
NG = abs(mesh_opt.NG);
NP = mesh_opt.NP;
BC = mesh_opt.BC;
NS = mesh_opt.NS;
% reactor structure
N_reg     = layer_struct.N_reg;
t_reg     = layer_struct.t_reg;
geom_type = layer_struct.geom_type;

% preallocation
dx        = zeros(N_reg,1); 
N         = dx;
prev_grid = [];
L         = zeros(N_reg,NG);

if NP<0
    NP     = -NP;
    nonuni = 1;
else
    nonuni = 0;
end


for ii = 1:N_reg
    
    for g=1:NG
        DFC     = multig_data{1,ii}.DIFFCOEF(g);
        XS_REM  = multig_data{1,ii}.XS_REM(g);
        xs      = linspace(t_reg(ii),t_reg(ii+1),NS); % sampling grid
        L(ii,g) = min(sqrt((DFC./f_shape{ii}(xs))./(XS_REM.*f_shape{ii}(xs))));
    end
    
end
% Grid definition for each region according to its diffusion length (kind of adaptive grid)

%% ======================================= Spatial mesh ==========================================================================================================
% build a mesh per each region
for ii = 1:N_reg
    if nonuni==1
        dx(ii)    = min(L(ii,:))/NP; % NP>10 should guarantee an adequate spatial grid
    else
        dx(ii)    = min(min(L))/NP;
    end
    
    N(ii)     = ceil((t_reg(ii+1)-t_reg(ii))/dx(ii));   % Number of grid spatial points
    new_grid  = linspace(t_reg(ii),t_reg(ii+1),N(ii));
    FD_grid   = [prev_grid new_grid];
    prev_grid = FD_grid;
    dx(ii)    = new_grid(2)-new_grid(1);
    
    if length(unique(FD_grid))<length(FD_grid)
        N(ii) = N(ii)-1;
        FD_grid = unique(FD_grid);
    end
    
end

FD_grid = unique(FD_grid);
pos = find(FD_grid == 0);
if pos
   FD_grid(pos)=eps;
end
% Mesh struct definition

fd_data = struct('FD_grid',FD_grid,'dx',dx,'N',N,'BC',BC,'geom',geom_type,'N_reg',N_reg,'t_reg',t_reg);
fd_data.f = f_shape; % to avoid Matlab bug that divides the cell into many structs
end
