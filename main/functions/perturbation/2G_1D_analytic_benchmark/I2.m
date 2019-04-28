function  O2 = I2(a,b,i,Bn)
% evaluate the integral of one sinusoidal function
O2 = -1./Bn(i).*(cos(Bn(i)*b)-cos(Bn(i)*a));
end