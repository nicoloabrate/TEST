function data = readnucdata(mat_name,NG,mat_id,T)
% This function reads data from small the library stored in
% data_library.

if nargin == 4
    % Export serpent output data to file
    GCU_name = serp2file(mat_name,NG,mat_id,T);
    % Load data from formatted file
    XS = load([pwd,filesep,'auxiliary',filesep,'data_library',filesep,num2str(NG),'G',filesep,'XS',filesep,GCU_name]);
    data = struct('XS_TOT',XS(1,:),'XS_TRANSPXS',XS(2,:),'DIFFCOEF',XS(3,:),...
                  'XS_ABS',XS(4,:),'XS_CAPT',XS(5,:),'XS_FISS',XS(6,:),...
                  'XS_REM',XS(7,:),'CHIT',XS(8,:),'XS_NSF',XS(9,:),...
                  'XS_S0',XS(10:10+NG-1,:));
elseif nargin == 3
    % Load data from serpent output
    core_data = serp2struct(mat_name,NG,mat_id);
else
    % Load data from formatted file
    XS = load([pwd,'/TEST_data_library/',num2str(NG),'G/XS/',mat_name]);
    data = struct('XS_TOT',XS(1,:),'XS_TRANSPXS',XS(2,:),'DIFFCOEF',XS(3,:),...
                  'XS_ABS',XS(4,:),'XS_CAPT',XS(5,:),'XS_FISS',XS(6,:),...
                  'XS_REM',XS(7,:),'CHIT',XS(8,:),'XS_NSF',XS(9,:),...
                  'XS_S0',XS(10:10+NG-1,:));

end

end
