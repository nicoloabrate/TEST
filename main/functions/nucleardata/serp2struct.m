function [data] = serp2struct(fname,mat_id,NG)
% This file stores nuclear constants evaluated with a macro-group structure
% from Serpent-2 Monte Carlo code to Matlab structures.
% fname ====> Atomic symbol and Mass number for nuclide, Common name for
% molecules and alloys;
% mat_id   ====> to be removed from future releases

NG = abs(NG);
path = pwd;
run([pwd,filesep,'auxiliary',filesep,'data_library',filesep,num2str(NG),'G/serpent_output/',fname,'_res.m'])
cd([pwd,filesep,'auxiliary',filesep,'data_library',filesep,num2str(NG),'G/XS'])


S = 0;
if sum(INF_CHIT(mat_id,1:2:end)) > 1
    for ii=1:length(INF_CHIT(mat_id,1:2:end))
        S = S+INF_CHIT(mat_id,ii);
        if S==1
            INF_CHIT(mat_id,ii+1)=0;
        end
    end
end


clear NUBAR
NUBAR = INF_NSF(mat_id,1:2:end)./INF_FISS(mat_id,1:2:end); 


if length(FISSE(mat_id,1:2:end))<NG
    EF = NaN*ones(NG,1);
    for ii = 1:NG
       EF(ii) = FISSE(mat_id,1:2:end); 
    end
end

data = struct('DIFFCOEF',INF_DIFFCOEF(mat_id,1:2:end),...
              'XS_TOT',INF_TOT(mat_id,1:2:end),...
              'XS_TRANSPXS',INF_TRANSPXS(mat_id,1:2:end),...
              'XS_ABS',INF_ABS(mat_id,1:2:end),...
              'XS_CAPT',INF_CAPT(mat_id,1:2:end),...
              'XS_FISS',INF_FISS(mat_id,1:2:end),...
              'XS_REM',INF_REMXS(mat_id,1:2:end),...
              'XS_S0',reshape(INF_S0(mat_id,1:2:end), NG, NG),...
              'XS_NSF',INF_NSF(mat_id,1:2:end),...
              'CHIT',INF_CHIT(mat_id,1:2:end),...
              'NUBAR',NUBAR,...
              'EF',EF);

% Come back to older position
cd(path)
end
