%% kuramoto runs for paper

%% add subdirectories to path

addpath(genpath('Kuramoto/'));
addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));

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
pin = [0.9,0.6]; % coupling probability for in-community nodes
pout = 0.01;     % coupling probability for out-of-community nodes
pbase = pout;    % coupling probability for all singletons

% define parameters: community detection
p = 100;
gamma_ = 1;
omega_ = 0;



Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

toRemove = 1:40;
[Ra2,Fa2,R2,~,Rcons2,~,A_ens2,C_ens2] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true);
toRemove = 1:60;
[Ra3,Fa3,R3,~,Rcons3,~,A_ens3,C_ens3] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true,...
                       'doFullRun',false);                   
toRemove = 61:100;
[Ra1,Fa1,R1,F1,Rcons1,Fcons1,A_ens1,C_ens1] = remove_CD_comp(Cgenfun,toRemove,sims,'endtime',endtime,'T',1,...
                       'p',p,'gamma_',gamma_,'omega_',omega_,'makePlot',true,...
                       'doFullRun',false);
