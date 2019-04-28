function output_data = add_perturbation(pert_mat,input_data)
% List of reactions that can be perturbed for a given set of data
% XS_FISS = 1
% XS_CAPT = 2
% NUBAR   = 3
% XS_S0   = 4
XS     = pert_mat(:,1); % data to be perturbed
NG_in  = pert_mat(:,2); % perturbed in groups
NG_out = pert_mat(:,3); % perturbed out groups
delta  = pert_mat(:,4); % perturbation amplitude

N_pert = size(pert_mat,1);

for ii = 1:N_pert
    
    if XS(ii) == 1
        % Evaluate XS correction
        XS_corr = delta(ii)*1e-2*input_data.XS_FISS(NG_in(ii));
        % Perturb fission XS
        input_data.XS_FISS(NG_in(ii)) = (1+delta(ii)*1e-2)*input_data.XS_FISS(NG_in(ii));
        % Correct nuSf XS
        input_data.XS_NSF(NG_in(ii))  = input_data.NUBAR(NG_in(ii))*input_data.XS_FISS(NG_in(ii));
        NG = NG_in(ii);
        
    elseif XS(ii) == 2
        % Evaluate XS correction
        XS_corr = delta(ii)*1e-2*input_data.XS_CAPT(NG_in(ii));
        input_data.XS_CAPT(NG_in(ii)) = (1+delta(ii)*1e-2)*input_data.XS_CAPT(NG_in(ii));
        NG = NG_in(ii);
        
    elseif XS(ii) == 3
        
        input_data.NUBAR(NG_in(ii))  = (1+delta(ii)*1e-2)*input_data.NUBAR(NG_in(ii));
        input_data.XS_NSF(NG_in(ii)) = input_data.NUBAR(NG_in(ii))*input_data.XS_FISS(NG_in(ii));
        NG = NG_in(ii);
        
    elseif XS(ii) == 4
        % Evaluate XS correction
        XS_corr = delta(ii)*1e-2*input_data.XS_S0(NG_in(ii),NG_out(ii));
        input_data.XS_S0(NG_in(ii),NG_out(ii)) = (1+delta(ii)*1e-2)*input_data.XS_S0(NG_in(ii),NG_out(ii));
        % Group number to correct REM XS
        NG = NG_out(ii);
    end
    
    % Edit XS to ensure normalization constraint
    if XS(ii) ~= 3
        input_data.XS_ABS(NG)       = input_data.XS_ABS(NG)+XS_corr;
        input_data.XS_TOT(NG)       = input_data.XS_TOT(NG)+XS_corr;
        input_data.XS_REM(NG)       = input_data.XS_REM(NG)+XS_corr;
        input_data.XS_TRANSPXS(NG)  = input_data.XS_TRANSPXS(NG)+XS_corr;
        input_data.DIFFCOEF(NG)     = 1/(3.*(input_data.XS_TRANSPXS(NG)));
    end
end
output_data = input_data;
end