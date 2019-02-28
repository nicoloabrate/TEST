function fshape = define_shape(N_reg,t_reg)
% This function generates the function handles to evaluate local nuclear
% properties.

func = [];
for ii=1:N_reg
    if ii==1
        f = (sprintf('@(x,XS) (%s).*XS(%d).*(x>%f & x<=%f)',shape(ii),ii,t_reg(ii),t_reg(ii+1)));
    else
        f = (sprintf('+(%s).*XS(%d).*(x>%f & x<=%f)',shape(ii),ii,t_reg(ii),t_reg(ii+1)));
    end
    func = [func f]; % strings concatenation
end
fshape = eval(func);
end