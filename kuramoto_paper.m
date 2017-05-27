%% kuramoto runs for paper

%% add subdirectories to path

addpath(genpath('Kuramoto/'));
addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));

%% Try with multiple time windows. First.

sims = 8;
endtime = 40;
T = 8;
%kappa_ = [0.2 0.2 0.3 0.3 0.5 0.5 0.2 0.2]; % change over time
kappa_ = 0.2;
%kappa_ = [0.2;0.5];

% structural parameters
N = 100;
M = [20;8];   % unique community sizes present 
m = [3;5];    % number of communities of each size
pin = [0.9;0.9]; % coupling probability for in-community nodes
pout = 0.01;    % coupling probability for out-of-community nodes
pbase = pout;

% CD parameters
gamma_ = 1;
omega_ = 0.5;
saveString = 'ss8run_g1o05';

Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');
Cdeepfun = @()modcoupler(N,M,m,0,1,0);

% kappa matrix
if all(size(kappa_)==size(M))
    kappa_ = make_kappa_mat(kappa_,N,M,m);
end

% try full runs
toRemove = [];
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');

%%
% try removing 20 comms...
toRemove = 1:60;
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '_rem20.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');

% now try removing 8 comms
toRemove = 61:100;
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '_rem8.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');

%% calculate Truepos and Falsepos
[tPos,fPos] = calculate_single_acc(R,Cens);

figure; scatter(mean(fPos),mean(tPos),'filled'); hold on;
        scatter(mean(mean(fPos)),mean(mean(tPos)),'filled');
        xlabel('false positive rate');
        ylabel('true positive rate');

%% ROC curve
[tPosRates,fPosRates] = calculate_thresh_acc(Rass,Cdeepfun()+eye(N),sims);

figure; plot(mean(fPosRates,2),mean(tPosRates,2),'-k'); hold on;
      for t=1:T
        plot(fPosRates(:,t),tPosRates(:,t),'--');
      end
        xlabel('false positive rate');
        ylabel('true positive rate');
        axis([0 1 0 1]);

%% QaD plot...
doFullRun = 0;
Rcons = Rg;
nR = N - numel(toRemove);
p = 100;

for s=1:5
      % create nx(p*T) matrix
    if doFullRun  
      pFmat = zeros(n,p*T);
    end
      pRmat = zeros(nR,p*T);
      for t=1:T
        if doFullRun  
          pFmat(:,(t-1)*p+1:t*p) = Fcons(:,:,t,s)';
        end
          pRmat(:,(t-1)*p+1:t*p) = Rcons(:,:,t,s)';
      end
    if doFullRun  
      ncommsF = numel(unique(pFmat));
      figure; bcolor(pFmat); 
              colormap(lines(ncommsF)); caxis([1 ncommsF]); colorbar;
    end
      ncommsR = numel(unique(pRmat));
      figure; bcolor(pRmat); 
              colormap(lines(ncommsR)); caxis([1 ncommsR]); colorbar;        
end
        
        
%% Same-size runs. Goal: Test influence of intra-comm connectivity

sims = 10;
endtime = 40;

% define parameters: structural
N = 100;        % number of nodes
M = [20;20];   % unique community sizes present 
m = [2;3];    % number of communities of each size
pin = [0.9;0.85]; % coupling probability for in-community nodes
pout = 0.01;    % coupling probability for out-of-community nodes
pbase = pout;
toRemove = 1:40;

Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

[R,F] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,'makePlot',true);

%% Same-size runs, a new size

sims = 10;
endtime = 40;

% define parameters: structural
N = 100;        % number of nodes
M = [8;8];   % unique community sizes present 
m = [6;6];    % number of communities of each size
pin = [0.9;0.7]; % coupling probability for in-community nodes
pout = 0.01;    % coupling probability for out-of-community nodes
pbase = pout;
toRemove = 1:48;

Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

[R,F] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,'makePlot',true);

%% Runs for two community sizes...

sims = 10;
endtime = 40;

% define parameters: structural
N = 100;         % number of nodes
M = [20;8];      % unique community sizes present 
m = [3;5];       % number of communities of each size
pin = [0.9,0.7]; % coupling probability for in-community nodes
pout = 0.01;     % coupling probability for out-of-community nodes
pbase = pout;    % coupling probability for all singletons

% define parameters: community detection
p = 100;
gamma_ = 1;
omega_ = 0;



Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

%toRemove = 1:40;
%[Ra2,Fa2,R2,~,Rcons2,~,A_ens2,C_ens2] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
%                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true);
toRemove = 1:60;
[Ra3,Fa3,R3,~,Rcons3,Fcons3,A_ens3,C_ens3] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true);%,...
                       %'doFullRun',false);                   
%toRemove = 61:100;
%[Ra1,Fa1,R1,F1,Rcons1,Fcons1,A_ens1,C_ens1] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
%                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true,...
%                       'doFullRun',false);
