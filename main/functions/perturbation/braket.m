function I = braket(psi,phi,brakopt)
% This function evaluates the bra-ket of psi and phi vectors in phase
% space. The argument 'brakopt' contains the spatial grid and the
% information needed to perform the integral.

prod = psi.*phi;       % Element-wise product
NG   = brakopt.NG;     % Number of energy-groups
Nt   = brakopt.Nt;     % Number of spatial points in the grid
multi_grid = zeros(Nt*NG,1);   % Integral initialization

for g = 1:NG
    multi_grid(1+(g-1)*Nt:Nt+(g-1)*Nt) = brakopt.FD_grid;
end

I = trapz(multi_grid,prod);

end