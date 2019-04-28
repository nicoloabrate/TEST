function  O1 = I3(a,b,i,j,Bn)
% evaluate the integral of the product of two sinusoidal functions
N = 1*(length(i)>length(j))+2*(length(i)==length(j))+3*(length(i)<length(j));

switch N
    case 1 % many evaluation along i are requested
        O1 = NaN*ones(1,length(i));
        for k=1:length(i)
            if (i(k)~=j)
                O1(k) = (sin((Bn(i(k))-Bn(j))*b)./(2*(Bn(i(k))-Bn(j)))+sin((Bn(i(k))+Bn(j))*b)./(2*(Bn(i(k))+Bn(j))))-...,
                        (sin((Bn(i(k))-Bn(j))*a)./(2*(Bn(i(k))-Bn(j)))+sin((Bn(i(k))+Bn(j))*a)./(2*(Bn(i(k))+Bn(j))));
            elseif (i(k)==j)
                O1(k) = 1/2*(b+sin(Bn(i(k))*b).*cos(Bn(j)*b)/Bn(i))-1/2*(a+sin(Bn(i(k))*a).*cos(Bn(j)*a)/Bn(i));
            end
        end
    case 2
        if (i~=j)
            O1 = (sin((Bn(i)-Bn(j))*b)./(2*(Bn(i)-Bn(j)))+sin((Bn(i)+Bn(j))*b)./(2*(Bn(i)+Bn(j))))-...,
                 (sin((Bn(i)-Bn(j))*a)./(2*(Bn(i)-Bn(j)))+sin((Bn(i)+Bn(j))*a)./(2*(Bn(i)+Bn(j))));
        elseif (i==j)
            O1 = 1/2*(b+sin(Bn(i)*b).*cos(Bn(j)*b)/Bn(i))-1/2*(a+sin(Bn(i)*a).*cos(Bn(j)*a)/Bn(i));
        end
    case 3 % many evaluation along j are requested
        O1 = NaN*ones(1,length(j));
        for k=1:length(j)
            if (i~=j(k))
                O1(k) = (sin((Bn(i)-Bn(j(k)))*b)./(2*(Bn(i)-Bn(j(k))))+sin((Bn(i)+Bn(j(k)))*b)./(2*(Bn(i)+Bn(j(k)))))-...,
                        (sin((Bn(i)-Bn(j(k)))*a)./(2*(Bn(i)-Bn(j(k))))+sin((Bn(i)+Bn(j(k)))*a)./(2*(Bn(i)+Bn(j(k)))));
            elseif (i==j(k))
                O1(k) = 1/2*(b+sin(Bn(i)*b).*cos(Bn(j(k))*b)/Bn(i))-1/2*(a+sin(Bn(i)*a).*cos(Bn(j(k))*a)/Bn(i));
            end
        end
end

end
