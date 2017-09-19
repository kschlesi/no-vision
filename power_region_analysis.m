% does community alignment with 'power' functional systems improve 
% after visual cortex removal? (Figures 8 and S1)

%% add subdirectories to path

addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));
addpath(genpath('Analysis_Utils/'));

%% load data
load('lifespan_rv.mat');
load('officer_rv.mat');
[n,t_OF,nsubs_OF] = size(partnS_rv_OF);
[~,t_LS,nsubs_LS] = size(partnS_rv_LS);

% load power region classifications for n=194 hybrid atlas
fid = fopen('hybrid_labels_long.txt');
regName = textscan(fid,'%s',n,'Delimiter','\n');
regName = regName{:};
fclose(fid);
[regClass,regCN] = region_class_map(regName,[pwd '/Region_Info/']);
nodeix = ~ismemvar(regClass,'Visual');
visix = ismemvar(regClass,'Visual');
motix = ismemvar(regClass,'Somatosensory');
nvmix = ~(visix+motix);

%% calculate recruitments for each node, time window, and subject...
recruS_v_OF = zeros(n,t_OF,nsubs_OF);
recruS_rv_OF = zeros(n,t_OF,nsubs_OF);
recruS_v_LS = zeros(n,t_LS,nsubs_LS);
recruS_rv_LS = zeros(n,t_LS,nsubs_LS);

% officer
for t=1:t_OF
  for k=1:nsubs_OF
    recruS_v_OF(:,t,k) = corecruitment(mod_allegiance(partnS_v_OF(:,t,k)'),regCN);
    recruS_rv_OF(:,t,k) = corecruitment(mod_allegiance(partnS_rv_OF(:,t,k)'),regCN);
  end
end

%lifespan
for t=1:t_LS
  for k=1:nsubs_LS
    recruS_v_LS(:,t,k) = corecruitment(mod_allegiance(partnS_v_LS(:,t,k)'),regCN);
    recruS_rv_LS(:,t,k) = corecruitment(mod_allegiance(partnS_rv_LS(:,t,k)'),regCN);
  end
end

%% overall non-visual recru v. tw, before and after vision removal
meanvOF = mean(recruS_v_OF(regCN~=10,:,:),1);
meanrvOF = mean(recruS_rv_OF(regCN~=10,:,:),1);

% plot: recruitment v. subject ID, for all 4 tws (tasks)
figure; subplot(1,2,1); plot(squeeze(meanvOF)'); axis([0 80 0.2 0.6]);
        subplot(1,2,2); plot(squeeze(meanrvOF)'); axis([0 80 0.2 0.6]);

meanvLS = mean(recruS_v_LS(regCN~=10,:,:),1);
meanrvLS = mean(recruS_rv_LS(regCN~=10,:,:),1);

% plot: recruitment v. subject ID, for all 3 tws (same task)
figure; subplot(1,2,1); plot(squeeze(meanvLS)'); axis([0 105 0.1 0.5]);
        subplot(1,2,2); plot(squeeze(meanrvLS)'); axis([0 105 0.1 0.5]);
        
%% compute system recruitment for each fcnal system, for each tw, for each subject
nC = numel(unique(regCN));
sysrecruS_v_OF = zeros(nC,t_OF,nsubs_OF);
sysrecruS_rv_OF = zeros(nC,t_OF,nsubs_OF);
sysrecruS_v_LS = zeros(nC,t_LS,nsubs_LS);
sysrecruS_rv_LS = zeros(nC,t_LS,nsubs_LS);

for c=unique(regCN)'
        sysrecruS_v_OF(c,:,:) = mean(recruS_v_OF(regCN==c,:,:),1);
        sysrecruS_rv_OF(c,:,:) = mean(recruS_rv_OF(regCN==c,:,:),1);
        sysrecruS_v_LS(c,:,:) = mean(recruS_v_LS(regCN==c,:,:),1);
        sysrecruS_rv_LS(c,:,:) = mean(recruS_rv_LS(regCN==c,:,:),1);
end

sysrecruS_rv_OF(10,:,:) = NaN;
sysrecruS_rv_LS(10,:,:) = NaN;

%% plot changes in each sys recru with removed vision (Figs. 8, S1)
figure;
for t=1:t_OF
    subplot(ceil(t_OF/2),2,t);
    bar_error([mean(sysrecruS_v_OF(:,t,:),3),...
                       mean(sysrecruS_rv_OF(:,t,:),3)],...
                      [std(sysrecruS_v_OF(:,t,:),[],3),...
                       std(sysrecruS_rv_OF(:,t,:),[],3)],...
                      'ExtFigureCmd','XTickLabel',unique(regClass),...
                      'legend',{'with vision','without vision'},...
                      'BarArgs',{'FaceColor','LineStyle','FaceAlpha';...
                                 'b','-',1;...
                                 'y','-.',0.5});
    title(['OF tw ' num2str(t)]); hold on;
end
figure;
for t=1:t_LS
    subplot(ceil(t_LS/2),2,t);
    bar_error([mean(sysrecruS_v_LS(:,t,:),3),...
                       mean(sysrecruS_rv_LS(:,t,:),3)],...
                      [std(sysrecruS_v_LS(:,t,:),[],3),...
                       std(sysrecruS_rv_LS(:,t,:),[],3)],...
                      'ExtFigureCmd','XTickLabel',unique(regClass),...
                      'legend',{'with vision','without vision'},...
                      'BarArgs',{'FaceColor','LineStyle','FaceAlpha';...
                                 'b','-',1;...
                                 'y','-.',0.5});
    title(['LS tw ' num2str(t)]); hold on;
end

%% in each time window, does removing vision correspond to a significant 
% increase (one-tailed test) in system recruitment distribution across subjects?

alpha = 0.05; % sig level

% compute p-vals for each system / time window
sys_pvals_OF = zeros(nC,t_OF);
sys_pvals_LS = zeros(nC,t_LS);
for c=1:nC
  for t=1:t_OF
    [~,p] = ttest(squeeze(sysrecruS_v_OF(c,t,:)),...
                  squeeze(sysrecruS_rv_OF(c,t,:)),...
                  'Tail','left','Alpha',alpha...
                  );
    sys_pvals_OF(c,t) = p;
  end
  for t=1:t_LS
    [~,p] = ttest(squeeze(sysrecruS_v_LS(c,t,:)),...
                  squeeze(sysrecruS_rv_LS(c,t,:)),...
                  'Tail','left','Alpha',alpha...
                  );    
    sys_pvals_LS(c,t) = p;
  end
end

% visualize signficant increases (using bonferroni correction)
plotLS = NaN.*zeros(nC,t_LS);
plotOF = NaN.*zeros(nC,t_OF);
plotLS(sys_pvals_LS<alpha/(nC*2-1)) = sys_pvals_LS(sys_pvals_LS<alpha/(nC*2-1));
plotOF(sys_pvals_OF<alpha/(nC*2-1)) = sys_pvals_OF(sys_pvals_OF<alpha/(nC*2-1));
figure; bcolor(-1*log(plotLS*(nC*2-1))'); colormap(winter); caxis([8 40]); colorbar;
figure; bcolor(-1*log(plotOF*(nC*2-1))'); colormap(winter); caxis([8 40]); colorbar;
