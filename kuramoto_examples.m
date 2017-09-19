%% example code to run kuramoto oscillator model

%% add subdirectories to path

addpath(genpath('Kuramoto/'));
addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));
addpath(genpath('Analysis_Utils/'));

%% Example 1: Influence matrix with single sizes (5-node comms).
%  Uses some of the utility functions in the appropriate subdiretories.

% simulation parameters
sims = 10;      % number of simulations of network
endtime = 40;   
T = 8;          % number of time slices
kappa_ = 0.2;   % time slice connection change
%kappa_ = [0.2 0.2 0.3 0.3 0.5 0.5 0.2 0.2]; % change over time

% structural parameters
N = 100;
M = [20;8];  % unique community sizes present 
m = [3;5];   % number of communities of each size
pin = [0.9;0.7]; % coupling probability for in-community nodes
pout = 0.01;     % coupling probability for out-of-community nodes
pbase = pout;

% CD parameters
gamma_ = 1;
omega_ = 1;
saveString = 's8nrun_g1o1_rem8_1';

% adjacency matrix-generating fcns
Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');
Cdeepfun = @()modcoupler(N,M,m,0,1,0);

% kappa matrix
if all(size(kappa_)==size(M))
    kappa_ = make_kappa_mat(kappa_,N,M,m);
end

% run full networks and save
toRemove = [];
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');


%% try removing 20 comms...
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

%% try removing ONE and TWO 20 comms...

% remove one 20-community
toRemove = 1:20;
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '_rem20_1.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');

% remove second 20-community                  
toRemove = 1:40;
[Rass,~,~,~,Rg,~,~,Cens] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',T,...
            'gamma_',gamma_,'omega_',omega_,'makePlot',true,'doFullRun',false,...
            'kappa_',kappa_);
p = 100; nR = N-numel(toRemove);
save(['Results/' saveString '_rem20_2.mat'],'endtime','gamma_','omega_','m','M','N','p','nR','toRemove',...
                      'Rass','Rg','pbase','pin','pout','sims','T','kappa_','Cens');                  
                  
%% calculate True positive rate and False positive rate
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
%  (not in paper)

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

