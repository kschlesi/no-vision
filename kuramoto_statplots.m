%% load things for plots
clear all;

paramString = 'arun20_g09o1';

load(['Results/' paramString '_rem20.mat']);
R20 = Rg;
Rass20 = Rass;
Cens20 = Cens;
nR20 = nR;
toRem20 = toRemove;

load(['Results/' paramString '_rem8.mat']);
R8 = Rg;
Rass8 = Rass;
Cens8 = Cens;
nR8 = nR;
toRem8 = toRemove;

load(['Results/' paramString '.mat']);
R = Rg;

%open('g1o05_Fass_t8_20-3_8-5.fig');
%open('g1o05_Fass_t8_20-3_8-5_rem8.fig');
%open('g1o05_Fass_t8_20-3_8-5_rem20.fig');

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
        xlabel('false positive rate');
        ylabel('true positive rate');
        legend('whole network',...
               '20-node comms in whole network',...
               '8-node comms in whole network',...
               '20-node comms alone',...
               '8-node comms alone',...
               'location','southeast');
        title(paramString);

%% ROC curve
Cdeepfun = @()modcoupler(N,M,m,0,1,0);
Cdeep = Cdeepfun()+eye(N);

[tPosRates,fPosRates] = calculate_thresh_acc(Rass,Cdeep,sims);
[tPosRates20w,fPosRates20w] = calculate_thresh_acc(Rass(toRem20,toRem20,:),Cdeep(toRem20,toRem20),sims);
[tPosRates8w,fPosRates8w] = calculate_thresh_acc(Rass(toRem8,toRem8,:),Cdeep(toRem8,toRem8),sims);
[tPosRates20r,fPosRates20r] = calculate_thresh_acc(Rass8,Cdeep(toRem20,toRem20),sims);
[tPosRates8r,fPosRates8r] = calculate_thresh_acc(Rass20,Cdeep(toRem8,toRem8),sims);

figure; plot(mean(fPosRates,2),mean(tPosRates,2),'-o','MarkerFaceColor',next_color()); hold on;
        plot(mean(fPosRates20w,2),mean(tPosRates20w,2),'-o','MarkerFaceColor',next_color());
        plot(mean(fPosRates8w,2),mean(tPosRates8w,2),'-o','MarkerFaceColor',next_color());
        plot(mean(fPosRates20r,2),mean(tPosRates20r,2),'--o','MarkerFaceColor',next_color());
        plot(mean(fPosRates8r,2),mean(tPosRates8r,2),'--o','MarkerFaceColor',next_color());
%       for t=1:T
%         plot(fPosRates(:,t),tPosRates(:,t),'--');
%       end
        xlabel('false positive rate');
        ylabel('true positive rate');
        axis([-0.05 1 0 1.05]);
        legend('whole network',...
               '20-node comms in whole network',...
               '8-node comms in whole network',...
               '20-node comms alone',...
               '8-node comms alone',...
               'location','southeast');
        title(paramString);
        
%% PR curve
Cdeepfun = @()modcoupler(N,M,m,0,1,0);
Cdeep = Cdeepfun()+eye(N);

[recall,~,~,~,precision] = calculate_thresh_acc(Rass,Cdeep,sims);
[recall20w,~,~,~,precision20w] = calculate_thresh_acc(Rass(toRem20,toRem20,:),Cdeep(toRem20,toRem20),sims);
[recall8w,~,~,~,precision8w] = calculate_thresh_acc(Rass(toRem8,toRem8,:),Cdeep(toRem8,toRem8),sims);
[recall20r,~,~,~,precision20r] = calculate_thresh_acc(Rass8,Cdeep(toRem20,toRem20),sims);
[recall8r,~,~,~,precision8r] = calculate_thresh_acc(Rass20,Cdeep(toRem8,toRem8),sims);

figure; plot(mean(recall,2),mean(precision,2),'-o','MarkerFaceColor',next_color()); hold on;
        plot(mean(recall20w,2),mean(precision20w,2),'-o','MarkerFaceColor',next_color());
        plot(mean(recall8w,2),mean(precision8w,2),'-o','MarkerFaceColor',next_color());
        plot(mean(recall20r,2),mean(precision20r,2),'--o','MarkerFaceColor',next_color());
        plot(mean(recall8r,2),mean(precision8r,2),'--o','MarkerFaceColor',next_color());
%       for t=1:T
%         plot(fPosRates(:,t),tPosRates(:,t),'--');
%       end
        xlabel('recall (# true pos / all true)');
        ylabel('precision (# true pos / all pos)');
        axis([0 1.05 0 1.05]);
        legend('whole network',...
               '20-node comms in whole network',...
               '8-node comms in whole network',...
               '20-node comms alone',...
               '8-node comms alone',...
               'location','southwest');
        title(paramString);

%% plots for all individual runs

strList = {'arun20_g1o1',...
           'aruns8_g1o1',...
           's8nrun_g1o1',...
           's8nrun_g1o1_rem8',...
           's8nrun_g1o1_rem20',...
           };
plotList = {'tpfp',...
            'PR'};

% plot setup
figNos = zeros(size(plotList));
for i=1:numel(plotList)
    figNos(i) = next_fig;
    figure(figNos(i)); hold on; hold off;
end
        
for i = 1:numel(strList)
    paramString = char(strList{i});
    disp(['Plotting ' paramString]);
    load(['Results/' paramString '.mat']);
    rix = removeval(1:N,toRemove);
    
    if any(ismemvar(plotList,'tpfp'))
    disp('...tpfp');    
    % create tpos/fpos plot for all selected runs
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        [tPos,fPos] = calculate_single_acc(Rg,Cens(rix,rix,:));
        scatter(mean(fPos),mean(tPos),'o','MarkerFaceColor',next_color);
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
    disp('...ROC');
    % create ROC curve for all selected runs
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        Cdeepfun = @()modcoupler(N,M,m,0,1,0);
        Cdeep = Cdeepfun()+eye(N);
        [tPosRates,fPosRates] = calculate_thresh_acc(Rass,Cdeep(rix,rix,:),sims);
        plot(mean(fPosRates,2),mean(tPosRates,2),'-o','MarkerFaceColor',next_color);
      hold off;
    end
    
    if any(ismemvar(plotList,'PR'))
    disp('...PR');
    % create PR curve for all selected runs
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        Cdeepfun = @()modcoupler(N,M,m,0,1,0);
        Cdeep = Cdeepfun()+eye(N);
        [recall,~,~,~,precision] = calculate_thresh_acc(Rass,Cdeep(rix,rix,:),sims);
        plot(mean(recall,2),mean(precision,2),'-o','MarkerFaceColor',next_color);
      hold off;
    end
end

% make legend list
%legList = cellfun(@(x)insertBefore(x,'_','\'),strList,'UniformOutput',false);
legList = {'20-node single-size ntwk','8-node single-size ntwk',...
           'full multi-res ntwk','20-node comms within full',...
           '8-node comms within full'};

% add plot details
for i=1:numel(plotList)
    
    if any(ismemvar(plotList,'tpfp'))
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        xlabel('false positive rate');
        ylabel('true positive rate');
        l = legend(legList,'location','southeast');
        title('CD performance by run');
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        xlabel('false positive rate');
        ylabel('true positive rate');
        axis([-0.05 1 0 1.05]);
        l = legend(legList,'location','southeast');
        title('ROC');
      hold off;  
    end
    
    if any(ismemvar(plotList,'PR'))
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        xlabel('recall (# true pos / all true)');
        ylabel('precision (# true pos / all pos)');
        axis([0 1.05 0 1.05]);
        l = legend(legList,'location','southwest');
        title('precision-recall curves');
      hold off;  
    end
    
end
        
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

%% create sysintro plot for some runs

runList = {'arun20_g1o5',...
           'aruns8_g1o1',...
           's8nrun_g1o1'
           };
for i=1:numel(runList)       
    paramString = char(runList{i});
    load(['Results/' paramString '.mat']);
    sys_intro_figure(Cens,Rass,'title',paramString);
end