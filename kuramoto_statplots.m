%% load things for plots: gamma=1, omega=0.5
clear all;

load('run_g1o05_rem20.mat');
R20 = Rg1o05;
Rass20 = Rassg1o05;
Cens20 = Cens;
nR20 = nR;
toRem20 = toRemove;

load('run_g1o05_rem8.mat');
R8 = Rg1o05;
Rass8 = Rassg1o05;
Cens8 = Cens;
nR8 = nR;
toRem8 = toRemove;

load('run_g1o05.mat');
R = Rg1o05;
Rass = Rassg1o05;

open('g1o05_Fass_t8_20-3_8-5.fig');
open('g1o05_Fass_t8_20-3_8-5_rem8.fig');
open('g1o05_Fass_t8_20-3_8-5_rem20.fig');

%% load things for plots: gamma=2, omega=0.5
clear all;

load('run_g2o05_rem20.mat');
R20 = Rg2o05;
Rass20 = Rassg2o05;
Cens20 = Cens;
nR20 = nR;
toRem20 = toRemove;

load('run_g2o05_rem8.mat');
R8 = Rg2o05;
Rass8 = Rassg2o05;
Cens8 = Cens;
nR8 = nR;
toRem8 = toRemove;

load('run_g2o05.mat');
R = Rg2o05;
Rass = Rassg2o05;

open('g2o05_Fass_t8_20-3_8-5.fig');
open('g2o05_Fass_t8_20-3_8-5_rem8.fig');
open('g2o05_Fass_t8_20-3_8-5_rem20.fig');

%% load things for plots: gamma=1, omega=0.05
clear all;

load('run_g1o005_rem20.mat');
R20 = Rg1o005;
Rass20 = Rassg1o005;
%Cens20 = Cens;
nR20 = nR;
toRem20 = toRemove;

load('run_g1o005_rem8.mat');
R8 = Rg1o005;
Rass8 = Rassg1o005;
%Cens8 = Cens;
nR8 = nR;
toRem8 = toRemove;

load('run_g1o005.mat');
R = Rg1o005;
Rass = Rassg1o005;

open('g1o005_Fass_t8_20-3_8-5.fig');
open('g1o005_Fass_t8_20-3_8-5_rem8.fig');
open('g1o005_Fass_t8_20-3_8-5_rem20.fig');

%% load things for plots: gamma=1, omega=0.1
clear all;

load('run_g1o01_rem20.mat');
R20 = Rg1o01;
Rass20 = Rassg1o01;
%Cens20 = Cens;
nR20 = nR;
toRem20 = toRemove;

load('run_g1o01_rem8.mat');
R8 = Rg1o01;
Rass8 = Rassg1o01;
%Cens8 = Cens;
nR8 = nR;
toRem8 = toRemove;

load('run_g1o01.mat');
R = Rg1o01;
Rass = Rassg1o01;

open('g1o01_Fass_t8_20-3_8-5.fig');
open('g1o01_Fass_t8_20-3_8-5_rem8.fig');
open('g1o01_Fass_t8_20-3_8-5_rem20.fig');

%% calculate Truepos and Falsepos
[tPos,fPos] = calculate_single_acc(R,Cens);
[tPos20w,fPos20w] = calculate_single_acc(R(:,toRem20,:,:),Cens(toRem20,toRem20,:));
[tPos8w,fPos8w] = calculate_single_acc(R(:,toRem8,:,:),Cens(toRem8,toRem8,:));
[tPos20r,fPos20r] = calculate_single_acc(R8,Cens8(toRem20,toRem20,:));
[tPos8r,fPos8r] = calculate_single_acc(R20,Cens20(toRem8,toRem8,:));

figure; scatter(mean(fPos),mean(tPos),'k','filled'); hold on;
        scatter(mean(fPos20w),mean(tPos20w),'b','filled'); 
        scatter(mean(fPos8w),mean(tPos8w),'r','filled'); 
        scatter(mean(fPos20r),mean(tPos20r),'bx'); 
        scatter(mean(fPos8r),mean(tPos8r),'rx');
        %scatter(mean(mean(fPos)),mean(mean(tPos)),'filled');
        xlabel('false positive rate');
        ylabel('true positive rate');
        legend('whole network',...
               '20-node comms in whole network',...
               '8-node comms in whole network',...
               '20-node comms alone',...
               '8-node comms alone',...
               'location','southeast');

%% ROC curve
Cdeepfun = @()modcoupler(N,M,m,0,1,0);
Cdeep = Cdeepfun()+eye(N);

[tPosRates,fPosRates] = calculate_thresh_acc(Rass,Cdeep,sims);
[tPosRates20w,fPosRates20w] = calculate_thresh_acc(Rass(toRem20,toRem20,:),Cdeep(toRem20,toRem20),sims);
[tPosRates8w,fPosRates8w] = calculate_thresh_acc(Rass(toRem8,toRem8,:),Cdeep(toRem8,toRem8),sims);
[tPosRates20r,fPosRates20r] = calculate_thresh_acc(Rass8,Cdeep(toRem20,toRem20),sims);
[tPosRates8r,fPosRates8r] = calculate_thresh_acc(Rass20,Cdeep(toRem8,toRem8),sims);

figure; plot(mean(fPosRates,2),mean(tPosRates,2),'-'); hold on;
        plot(mean(fPosRates20w,2),mean(tPosRates20w,2),'-');
        plot(mean(fPosRates8w,2),mean(tPosRates8w,2),'-');
        plot(mean(fPosRates20r,2),mean(tPosRates20r,2),'--');
        plot(mean(fPosRates8r,2),mean(tPosRates8r,2),'--');
%       for t=1:T
%         plot(fPosRates(:,t),tPosRates(:,t),'--');
%       end
        xlabel('false positive rate');
        ylabel('true positive rate');
        axis([0 1 0 1]);
        legend('whole network',...
               '20-node comms in whole network',...
               '8-node comms in whole network',...
               '20-node comms alone',...
               '8-node comms alone',...
               'location','southeast');

%% QaD plot...
doFullRun = 0;
Rcons = Rg1o05;
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