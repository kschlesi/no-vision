%% add subdirectories to path

addpath(genpath('Kuramoto/'));
addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));
addpath(genpath('Analysis_Utils/'));

%% Try with multiple time windows. First.
clear all
newSims = 10;

% for each saveString
strList = {'s8nrun_g1o1',...
           'arun20_g1o1',...
           's8nrun_g1o1_rem8_1',...
           's8nrun_g1o1_rem8',...
           'aruns8_g1o1',...
           's8nrun_g1o1_rem20_1',...
           's8nrun_g1o1_rem20_2',...
           's8nrun_g1o1_rem20',...
           };

for i=1:numel(strList)
    % set paramString
    paramString = char(strList{i});

    % load all structural params, dynamic params, CD params
    % leave out A_ens, Cens, Rass, Rg
    load(['Results/OrigFigs10/' paramString '.mat'],...
                                '-regexp','^((?!A_ens|Cens|Rass|Rg).)*$');
    totSims = sims + newSims;

    % setup networks
    Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');
    Cdeepfun = @()modcoupler(N,M,m,0,1,0);

    % kappa matrix
    if all(size(kappa_)==size(M))
        kappa_ = make_kappa_mat(kappa_,N,M,m);
    end

    % try full runs with given parameters
    [RassN,~,~,~,RgN,~,A_ensN,CensN] = remove_CD_comp(Cgenfun,toRemove,newSims,...
                'endtime',endtime,'T',T,'kappa_',kappa_,...
                'gamma_',gamma_,'omega_',omega_,...
                'makePlot',false,'doFullRun',false...
                );

    % load old results and concatenate
    load(['Results/OrigFigs10/' paramString '.mat'],'A_ens','Cens','Rass','Rg');
    % A_ens - 4 % Cens - 3 % Rg - 4 
    % Rass - recompute: (old Rass * old sims + new Rass*newSims) / totSims
    A_ens = cat(4,A_ens,A_ensN);
    Cens = cat(3,Cens,CensN);
    Rg = cat(4,Rg,RgN);
    Rass = (Rass*sims + RassN*newSims)/totSims;
    sims = totSims;

    % save all plus A_ens, ts, and sigma_
    save(['Results/' paramString '.mat'],'endtime','gamma_','omega_',...
                                            'm','M','N','p','nR','toRemove',...
                                            'pbase','pin','pout','sims',...
                                            'ts','sigma_','T','kappa_',...
                                            'Cens','A_ens','Rass','Rg');

    % clean up for next run
    clearvars -except strList newSims i
end                                    
