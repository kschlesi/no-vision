% a script for loading in a run, generating sample synchronization
% matrices, and adding them to the run files

% specify existing runs
runList = {'nrun_g1o05',...
           'nrun_g1o05_rem20',...
           'nrun_g1o05_rem8',...
           'nrun_g1o1',...
           'nrun_g1o1_rem20',...
           'nrun_g1o1_rem8',...
           'arun20_g1o1',...
           'arun20_g1o5',...
           'aruns8_g1o1',...
           'arunss8_g1o1',...
           's8nrun_g1o1',...
           's8nrun_g1o1_rem20',...
           's8nrun_g1o1_rem20_1',...
           's8nrun_g1o1_rem20_2',...
           's8nrun_g1o1_rem8',...
           's8nrun_g1o1_rem8_1',...
           'ss8nrun_g1o1',...
           'ss8nrun_g1o1_rem20',...
           'ss8nrun_g1o1_rem8'};

for i=1:numel(runList)
    paramString = char(runList{i});
    disp(['Running for ' paramString '...']);
    
    % open existing run
    load(['Results/' paramString '.mat']);
    Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

    % kappa matrix
    if all(size(kappa_)==size(M))
        kappa_ = make_kappa_mat(kappa_,N,M,m);
    end

    % run dynamics
    if ~exist('sigma_','var')
        sigma_ = 1; 
    end
    if ~exist('ts','var')
        ts = 0.1;
    end
    [~,A_ens,C_ens] = ksims_ens(sims,Cgenfun,kappa_,sigma_,ts,endtime);

    % append synchs to existing .mat file
    save(['Results/' paramString '.mat'],'A_ens','sigma_','ts','-append');
end
