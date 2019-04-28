function f_shape = adens_shape(N_reg,t_reg,f_def,geom_type)
%% This function takes the density function shapes and computes its normalization
% constants. Its output is an array of function handles ready to be used.

% Evaluation of region volumes
V  = zeros(length(N_reg),1);
NP = 1e3;

% preallocation
x         = zeros(N_reg,NP);
f_int     = zeros(N_reg,1);
f_handle  = cell(N_reg,1);
f_shape   = cell(N_reg,1);

% Evaluation of normalization constants (volume and shape function integral)
switch geom_type
    
    case 1 % slab
        
        for ii=1:N_reg
            
            V(ii)        = (t_reg(ii+1)-t_reg(ii));            % volume of the i-th layer
            x(ii,:)      = linspace(t_reg(ii),t_reg(ii+1),NP); % grid definition
            f_handle{ii} = eval(sprintf('@(x) %s',f_def{ii})); % build function handle
            fx           = f_handle{ii}(x(ii,:));              % function handle evaluation
            
            if length(fx)==length(x(ii,:))
                
                f_int(ii)    = trapz(x(ii,:),fx);              % integrate the function
                
            elseif fx>0
                
                f_int(ii)    = fx*(x(ii,end)-x(ii,1));         % integrate the function
                
            else
                
                f_int(ii) =  1;
                V(ii)     = -1;
            end
            
            f_shape{ii}  = eval(sprintf('@(x) (%s)*(%f)/(%f)',f_def{ii},V(ii),f_int(ii)));
            
            if f_int(ii)<=0
                error('The atomic density function needs to have a positive integral.')
            end
            
        end
        
    case 2 % sphere
        for ii=1:N_reg
            
            V(ii)        = (4/3*pi*(t_reg(ii+1)^3-t_reg(ii)^3));  % volume of the i-th layer
            x(ii,:)      = linspace(t_reg(ii),t_reg(ii+1),NP);    % grid definition
            f_handle{ii} = eval(sprintf('@(x) %s',f_def{ii}));    % build function handle
            fx           = f_handle{ii}(x(ii,:));                 % function handle evaluation
            
            if length(fx)==length(x(ii,:))
                
                f_int(ii)    = trapz(x(ii,:),f_handle{ii}(x(ii,:))*4*pi.*x.^2); % integrate the function
                
            elseif fx>0
                
                f_int(ii)    = fx*4/3*pi.*(x(ii,end).^3-x(ii,1).^3);         % integrate the function
                
            else
                
                f_int(ii) = 1;
                V(ii)     = -1;
            end
            
            f_shape{ii}  = eval(sprintf('@(x) (%s)*(%f)/(%f)',f_def{ii},V(ii),f_int(ii)));
            
            if f_int(ii)<=0
                error('The atomic density function needs to have a positive integral. \r\n')
            end
            
        end
        
    case 3 % cylinder
        
        for ii=1:N_reg
            
            V(ii)        = (pi*(t_reg(ii+1)^2-t_reg(ii)^2));     % volume of the i-th layer
            x(ii,:)      = linspace(t_reg(ii),t_reg(ii+1),NP);   % grid definition
            f_handle{ii} = eval(sprintf('@(x) %s',f_def{ii})); % build function handle
            fx           = f_handle{ii}(x(ii,:));                 % function handle evaluation
            
            if length(fx)==length(x(ii,:))
                
                f_int(ii) = trapz(x(ii,:),f_handle{ii}(x(ii,:))*2*pi.*x(ii,:));   % integrate the function
                
            elseif fx>0
                
                f_int(ii) = fx*pi.*(x(ii,end).^2-x(ii,1).^2);         % integrate the function
                
            else
                
                f_int(ii) = 1;
                V(ii)     = -1;
                
            end
            
            f_shape{ii}  = eval(sprintf('@(x) (%s)*(%f)/(%f)',f_def{ii},V(ii),f_int(ii)));
            
            if f_int(ii)<=0
                error('The atomic density function needs to have a positive integral. \r\n')
            end
            
        end
        
end

end
