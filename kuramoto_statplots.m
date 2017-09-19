
%% add subdirectories to path

addpath(genpath('Kuramoto/'));
addpath(genpath('CD_Utils/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('Region_Info/'));
addpath(genpath('Analysis_Utils/'));

colorOrder = get(gca,'ColorOrder');

%% plots for single v multi runs (Fig. 3)

strList = {'s8nrun_g1o1',...
           'arun20_g1o1',...
           's8nrun_g1o1_rem8',...
           'aruns8_g1o1',...
           's8nrun_g1o1_rem20',...
           };
plotList = {'tpfp'
            };
addAvg = false;

% plot setup
figNos = zeros(size(plotList));
for i=1:numel(plotList)
    figNos(i) = next_fig;
    figure(figNos(i)); hold on; hold off;
    colorOrder = get(gca,'ColorOrder');
end

colorList = {[0 0 0];
             colorOrder(1,:); % blue
             colorOrder(1,:); % blue
             [1 69/255 0]; % red/orange
             [1 69/255 0]; % red/orange
            };
dotSize = 12;        

for i = 1:numel(strList)
    paramString = char(strList{i});
    disp(['Plotting ' paramString]);
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
    rix = removeval(1:N,toRemove);
    if i==1
        RgM = Rg;
        CensM = Cens;
        RassM = Rass;
        RgAvgM = RgAvg;
        RassAvgM = RassAvg;
        continue;
    end
    
    if any(ismemvar(plotList,'tpfp'))
    disp('...tpfp');
    % create tpos/fpos plot for all selected runs
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
       if (i==3 || i==5)
        tpfp_figure(RgM(rix,rix,:,:),CensM(rix,rix,:),'ExtFigureCmd','DataOnly',...
                        'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.2},...
                        'plotargs',{'o','MarkerSize',dotSize,...
                                        'MarkerFaceColor',colorList{i},...
                                        'MarkerEdgeColor',colorList{i}});
        if addAvg                            
        tpfp_figure(RgAvgM(rix,rix,:,:),mean(CensM(rix,rix,:),3),'ExtFigureCmd','DataOnly',...
                        'plotargs',{'p','MarkerSize',10,...
                                        'MarkerFaceColor',colorList{i},...
                                        'MarkerEdgeColor',colorList{i}});
        end
       else
        tpfp_figure(Rg,Cens(rix,rix,:),'ExtFigureCmd','DataOnly',...
                        'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.2},...
                        'plotargs',{'o','MarkerSize',dotSize,...
                                     ...'MarkerFaceColor',last_color,...
                                        'MarkerEdgeColor',colorList{i}});
        if addAvg
        tpfp_figure(RgAvg,mean(Cens(rix,rix,:),3),'ExtFigureCmd','DataOnly',...
                        'plotargs',{'p','MarkerSize',10,...
                                     ...'MarkerFaceColor',last_color,...
                                        'MarkerEdgeColor',colorList{i}});                            
        end
       end
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
    disp('...ROC');
    % create ROC curve for all selected runs
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure(Rass,N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly');
      hold off;
    end
    
    if any(ismemvar(plotList,'PR'))
    disp('...PR');
    % create PR curve for all selected runs
      figure(figNos(strcmp(plotList,'PR'))); hold on;
       if (i==3 || i==5)
        prec_rec_figure(RassM(rix,rix,:),N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly',...
                        'plotargs',{':o','MarkerSize',8,...
                                    'Color',colorList{i}});...
                                    %'MarkerFaceColor',last_color(2)});
       else
        prec_rec_figure(Rass,N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly',...
                         'plotargs',{'-o','MarkerSize',8,...
                                     'Color',colorList{i},...
                                     'MarkerFaceColor',colorList{i}});
       end
      hold off;
    end

end

% make legend list
%legList = cellfun(@(x)insertBefore(x,'_','\'),strList,'UniformOutput',false);
legList = {'20-node comms in single-res network',...
           '20-node comms in multi-res network'..., no rem.','20-node comms in multi-res, 8-node comms rem.',...
           '8-node comms in single-res network',...
           '8-node comms in multi-res network'};%, no rem.','8-node comms in multi-res, 20-node comms rem.'};

% add plot details
for i=1:numel(plotList)
    
    if any(ismemvar(plotList,'tpfp'))
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        tpfp_figure([],[],'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
    if any(ismemvar(plotList,'PR'))
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        prec_rec_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
end
        
%% plots for removal runs: remove large comms (Fig. 4)

strList = {'s8nrun_g1o1_rem20',...      %  remove all 20s
           's8nrun_g1o1_rem20_2',...    %  remove 2 20s
           's8nrun_g1o1_rem20_1',...    %  remove 1 20
           's8nrun_g1o1',...           %  full multi-resolution network
           };
plotList = {'tpfp'};
msizes = [14,12,10,8];
tpfps = zeros(20,2,numel(strList));

% plot setup
figNos = zeros(size(plotList));
for i=1:numel(plotList)
    figNos(i) = next_fig;
    figure(figNos(i)); hold on; hold off;
    colorOrder = get(gca,'ColorOrder');
end

colorList = {%[1 20/255 147/255]; % pink
             %[1 0 0.2] % red
             %colorOrder(4,:); % purple
             [0.7 0.2 0.55] % light purple
             %colorOrder(5,:); % green
             %[0.8 0.1 0.2] % deep red
             colorOrder(5,:); % green
             colorOrder(3,:); % yellow
             %colorOrder(7,:); % maroon
             [1 69/255 0]; % red/orange
            };

% plotting loop        
for i = 1:numel(strList)
    paramString = char(strList{i});
    disp(['Plotting ' paramString]);
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
    rix = removeval(1:N,toRemove);
    
    if any(ismemvar(plotList,'tpfp'))
    disp('...tpfp');    
    % create tpos/fpos plot for all selected runs
      msize = msizes(i);
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        plotargs = {'o','MarkerSize',msize,...
                        'MarkerEdgeColor',colorList{i},...
                        'MarkerFaceColor',colorList{i},...
                        };
        [tps,fps] = tpfp_figure(Rg(:,end-40+1:end,:,:),...
                    Cens(end-40+1:end,end-40+1:end,:),...
                    'plotargs',plotargs,...
                    'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.25},...
                    'ExtFigureCmd','DataOnly');
        tpfps(:,:,i) = [tps(:),fps(:)];
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
    disp('...ROC');
    % create ROC curve for all selected runs
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure(Rass(end-40+1:end,end-40+1:end),...
                   N,M,m,1:60,sims,'ExtFigureCmd','DataOnly');
      hold off;
    end
    
    if any(ismemvar(plotList,'PR'))
    disp('...PR');
    % create PR curve for all selected runs
      msize = msizes(i);
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        plotargs = {'o:','MarkerSize',msize,...
                         'MarkerEdgeColor',colorList{i},...
                         'MarkerFaceColor',colorList{i},...
                         };
        prec_rec_figure(Rass(end-40+1:end,end-40+1:end),...
                        N,M,m,1:60,sims,...
                        'plotargs',plotargs,...
                        'ExtFigureCmd','DataOnly');
      hold off;
    end
end

% make legend list
%legList = cellfun(@(x)insertBefore(x,'_','\'),strList,'UniformOutput',false);
legList = {'8-node communities in full multi-res network',...
           '8-node communities with one 20-node comm. removed',...
           '8-node communities with two 20-node comms. removed',...
           '8-node communities with all 20-node comms. removed'...
           };

% add plot details
for i=1:numel(plotList)
    
    if any(ismemvar(plotList,'tpfp'))
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        tpfp_figure([],[],'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title');%,'legend',legList);
      hold off;  
    end
    
    if any(ismemvar(plotList,'PR'))
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        prec_rec_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
end

%% plots for removal runs: remove small comms (Fig. 5)

strList = {'s8nrun_g1o1_rem8',...      %  remove all 8s
           's8nrun_g1o1_rem8_1',...    %  remove 1 8
           's8nrun_g1o1',...           %  full multi-resolution network
           };
plotList = {'tpfp'};
msizes = [13,11,9];
tpfps = zeros(20,2,numel(strList));

% plot setup
figNos = zeros(size(plotList));
for i=1:numel(plotList)
    figNos(i) = next_fig;
    figure(figNos(i)); hold on; hold off;
    colorOrder = get(gca,'ColorOrder');
end

colorList = {
             colorOrder(4,:); % purple
             colorOrder(5,:); % green
             colorOrder(1,:); % blue
            };

for i = 1:numel(strList)
    paramString = char(strList{i});
    disp(['Plotting ' paramString]);
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
    rix = removeval(1:N,toRemove);
    
    if any(ismemvar(plotList,'tpfp'))
    disp('...tpfp');    
    % create tpos/fpos plot for all selected runs
      msize = msizes(i);
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        plotargs = {'o','MarkerSize',msize,...
                        'MarkerEdgeColor',colorList{i},...
                        'MarkerFaceColor',colorList{i},...
                        };
        [tps,fps] = tpfp_figure(Rg(:,1:60,:,:),...
                    Cens(1:60,1:60,:),...
                    'plotargs',plotargs,...
                    'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.25},...
                    'ExtFigureCmd','DataOnly');
        tpfps(:,:,i) = [tps(:),fps(:)];
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
    disp('...ROC');
    % create ROC curve for all selected runs
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure(Rass(1:60,1:60),...
                   N,M,m,60:end,sims,'ExtFigureCmd','DataOnly');
      hold off;
    end
    
    if any(ismemvar(plotList,'PR'))
    disp('...PR');
    % create PR curve for all selected runs
      msize = msizes(i);
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        plotargs = {'o:','MarkerSize',msize,...
                         'MarkerEdgeColor',colorList{i},...
                         'MarkerFaceColor',colorList{i},...
                         };
        prec_rec_figure(Rass(1:60,1:60),...
                        N,M,m,60:end,sims,...
                        'plotargs',plotargs,...
                        'ExtFigureCmd','DataOnly');
      hold off;
    end
end

% make legend list
%legList = cellfun(@(x)insertBefore(x,'_','\'),strList,'UniformOutput',false);
legList = {'20-node comms with all 8-node comms removed',...
           '20-node comms with two 8-node comms removed',...
           '20-node comms in full ntwk',...
           };

% add plot details
for i=1:numel(plotList)
    
    if any(ismemvar(plotList,'tpfp'))
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        tpfp_figure([],[],'ExtFigureCmd','AccOnly','title');%,'legend',legList);
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title');%,'legend',legList);
      hold off;  
    end
    
    if any(ismemvar(plotList,'PR'))
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        prec_rec_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
end

%% statistical tests of the null hypothesis that means of TPR and FPR are 
%  the same between different configurations

ndists = 3;             % number of dists to compare
alpha = 0.05;           % significance level
k2mat = zeros(ndists);  
tpmat = zeros(ndists);
fpmat = zeros(ndists);
for i=1:ndists
    for j=1:ndists
        [~,pValue] = kstest_2s_2d(tpfps(:,:,i),tpfps(:,:,j));
        k2mat(i,j) = pValue;
        [~,pValue] = ttest2(tpfps(:,1,i),tpfps(:,1,j));
        tpmat(i,j) = pValue;
        [~,pValue] = ttest2(tpfps(:,2,i),tpfps(:,2,j));
        fpmat(i,j) = pValue;
    end
end

disp(tpmat);
disp(fpmat);
disp(k2mat);

% plot significant pairs (after bonferroni correction for mult comps)
figure; bcolor(tpmat<(alpha/(ndists*(ndists-1)/2))); title('true pos t-test');
figure; bcolor(fpmat<(alpha/(ndists*(ndists-1)/2))); title('false pos t-test');
figure; bcolor(k2mat<(alpha/(ndists*(ndists-1)/2))); title('2d KS test');

%% plots for all individual runs: Summary Figure (not in paper)

strList = {'s8nrun_g1o1',...
           'arun20_g1o1',...
           's8nrun_g1o1_rem8',...
           'aruns8_g1o1',...
           's8nrun_g1o1_rem20',...
           };
plotList = {'tpfp'
            };

% plot setup
figNos = zeros(size(plotList));
for i=1:numel(plotList)
    figNos(i) = next_fig;
    figure(figNos(i)); hold on; hold off;
    colorOrder = get(gca,'ColorOrder');
end

% color setup
colorList = {[0 0 0];
             colorOrder(1,:); % blue
             colorOrder(1,:); % blue
             [1 69/255 0]; % red/orange
             [1 69/255 0]; % red/orange
            };

% plotting loop
for i = 1:numel(strList)
    paramString = char(strList{i});
    disp(['Plotting ' paramString]);
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
    rix = removeval(1:N,toRemove);
    if i==1
        RgM = Rg;
        CensM = Cens;
        RassM = Rass;
        continue;
    end
    
    if any(ismemvar(plotList,'tpfp'))
    disp('...tpfp');
    % create tpos/fpos plot for all selected runs
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
       if (i==3 || i==5)
        tpfp_figure(RgM(rix,rix,:,:),CensM(rix,rix,:),'ExtFigureCmd','DataOnly',...
                        'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.2},...
                        'plotargs',{'o','MarkerFaceColor',colorList{i},...
                                        'MarkerEdgeColor',colorList{i},...
                                        'MarkerSize',6});        
        tpfp_figure(Rg,Cens(rix,rix,:),'ExtFigureCmd','DataOnly',...
                        'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.2},...
                        'plotargs',{'^',...%'MarkerFaceColor',colorList{i},...
                                        'MarkerEdgeColor',colorList{i},...
                                        'MarkerSize',8});
       else
        tpfp_figure(Rg,Cens(rix,rix,:),'ExtFigureCmd','DataOnly',...
                        'addEllipse','ellipseArgs',{'NoLine',...
                                     'ellipseColor',colorList{i},...
                                     'ellipseAlpha',0.2},...
                        'plotargs',{'o',...'MarkerFaceColor',last_color,...
                                        'MarkerEdgeColor',colorList{i},...
                                        'MarkerSize',7});
       end
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
    disp('...ROC');
    % create ROC curve for all selected runs
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure(Rass,N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly');
      hold off;
    end
    
    if any(ismemvar(plotList,'PR'))
    disp('...PR');
    % create PR curve for all selected runs
      figure(figNos(strcmp(plotList,'PR'))); hold on;
       if (i==3 || i==5)
        prec_rec_figure(RassM(rix,rix,:),N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly',...
                        'plotargs',{':o','Color',colorList{i},...
                                    'MarkerFaceColor',colorList{i},...
                                        'MarkerSize',8});
        prec_rec_figure(Rass,N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly',...
                         'plotargs',{'--*','Color',colorList{i},...
                                     'MarkerFaceColor',colorList{i},...
                                        'MarkerSize',8});
       else
        prec_rec_figure(Rass,N,M,m,toRemove,sims,'ExtFigureCmd','DataOnly',...
                         'plotargs',{'-o','Color',colorList{i},...
                                        'MarkerSize',8});
                                     ...'MarkerFaceColor',last_color});
       end
      hold off;
    end

end

% make legend list
%legList = cellfun(@(x)insertBefore(x,'_','\'),strList,'UniformOutput',false);
legList = {'20-node single-res ntwk',...
           '20-node comms in multi-res, no rem.','20-node comms in multi-res, 8-node comms rem.',...
           '8-node single-res ntwk',...
           '8-node comms in multi-res, no rem.','8-node comms in multi-res, 20-node comms rem.'};

% add plot details
for i=1:numel(plotList)
    
    if any(ismemvar(plotList,'tpfp'))
      figure(figNos(strcmp(plotList,'tpfp'))); hold on;
        tpfp_figure([],[],'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;
    end
    
    if any(ismemvar(plotList,'ROC'))
      figure(figNos(strcmp(plotList,'ROC'))); hold on;
        ROC_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
    if any(ismemvar(plotList,'PR'))
      figure(figNos(strcmp(plotList,'PR'))); hold on;
        prec_rec_figure([],[],[],[],[],[],...
                   'ExtFigureCmd','AccOnly','title','legend',legList);
      hold off;  
    end
    
end

%% network visualization (for Figs. 3, 4, 5)

runList = {'s8nrun_g1o1',...         % full multi-resolution network
          's8nrun_g1o1_rem8_1',...    %  remove 1 20
          's8nrun_g1o1_rem8',...    %  remove 2 20s
          };      
      
for i=1:numel(runList)       
    paramString = char(runList{i});
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
    
    Cens_mask = mask_remove(Cens,toRemove);
    Rass_addback = addback_remove(addback_remove(Rass,N,toRemove,1),N,toRemove,2);
    
    figure(next_fig);
    structmat_figure(Cens_mask,'ExtFigureCmd');
    
    figure(next_fig);
    load(['Results/' paramString '.mat'],'A_ens');
    A = make_synchs(A_ens,'T',8);
    clear A_ens;
    synchmat_figure(mask_remove(A{8,1},toRemove),'ExtFigureCmd');
    
    figure(next_fig);
    det_prob_tAvg_figure(Rass_addback,'ExtFigureCmd','colorbar');

    %sys_intro_figure(Cens_mask,Rass_addback,'title',paramString);
end
