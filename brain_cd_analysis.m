%% removing vision from brain networks: analysis
% (Figs. 6, 7, and 9)

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

% load atlas classifications for n=194 hybrid atlas
fid = fopen('hybrid_labels_long.txt');
regName = textscan(fid,'%s',n,'Delimiter','\n');
regName = regName{:};
fclose(fid);
[regClass,regCN] = region_class_map(regName,[pwd '/Region_Info/']);
nodeix = ~ismemvar(regClass,'Visual');
visix = ismemvar(regClass,'Visual');
motix = ismemvar(regClass,'Somatosensory');
nvmix = ~(visix+motix);

%% visualize overall node-specific flexibility (Fig. 7)
figure;
    colorOrder = get(gca,'ColorOrder');
    correlate(mean(flex_v_LS(~~visix,:),2),mean(flex_v_OF(~~visix,:),2),...
        colorOrder(1,:),'fill','ExtFigureCmd'); hold on;
    correlate(mean(flex_v_LS(~~motix,:),2),mean(flex_v_OF(~~motix,:),2),...
        colorOrder(2,:),'fill','ExtFigureCmd');
    correlate(mean(flex_v_LS(~~nvmix,:),2),mean(flex_v_OF(~~nvmix,:),2),...
        colorOrder(3,:),'fill','ExtFigureCmd');
    filled_ellipse(std(mean(flex_v_LS(~~visix,:),2),[],1),...
                   std(mean(flex_v_OF(~~visix,:),2),[],1),...
                   0,...
                   mean(mean(flex_v_LS(~~visix,:),2),1),...
                   mean(mean(flex_v_OF(~~visix,:),2),1),...
                   'ExtFigureCmd','ellipseAlpha',0.13,...
                   'ellipseColor',colorOrder(1,:),'NoLine');
    filled_ellipse(std(mean(flex_v_LS(~~motix,:),2),[],1),...
                   std(mean(flex_v_OF(~~motix,:),2),[],1),...
                   0,...
                   mean(mean(flex_v_LS(~~motix,:),2),1),...
                   mean(mean(flex_v_OF(~~motix,:),2),1),...
                   'ExtFigureCmd','ellipseAlpha',0.13,...
                   'ellipseColor',colorOrder(2,:),'NoLine');
    filled_ellipse(std(mean(flex_v_LS(~~nvmix,:),2),[],1),...
                   std(mean(flex_v_OF(~~nvmix,:),2),[],1),...
                   0,...
                   mean(mean(flex_v_LS(~~nvmix,:),2),1),...
                   mean(mean(flex_v_OF(~~nvmix,:),2),1),...
                   'ExtFigureCmd','ellipseAlpha',0.13,...
                   'ellipseColor',colorOrder(3,:),'NoLine');           
    xlabel('mean node flexibility, single-task data');
    ylabel('mean node flexibility, multi-task data');
    legend('visual cortex','motor cortex','other','location','southeast');
    
%% histograms of mean flexibility (Fig. 9)
nbins = 15;
figure;
    subplot(2,1,1);
        histogram(mean(flex_v_LS(~~nodeix,:),1),nbins,...
            'EdgeColor','k','FaceColor','b',...
            'LineStyle','-'); hold on;
            %'EdgeColor','b','FaceColor','none','DisplayStyle','stairs',...
        histogram(mean(flex_rv_LS(~~nodeix,:),1),nbins,...
            'EdgeColor','k','FaceColor','y',...
            'LineStyle','-.','FaceAlpha',0.5);
            %'EdgeColor','b','FaceColor','none','DisplayStyle','stairs',...
        axis([0.1,0.8,0,26]);
        xlabel('whole-brain flexibility'); ylabel('number of subjects');
        title('Single-task Data');
        legend('vision','no vision',...
               'location','northeast');
    subplot(2,1,2); 
        histogram(mean(flex_v_OF(~~nodeix,:),1),nbins,...
            'EdgeColor','k','FaceColor','b',...
            'LineStyle','-'); hold on;
            %'EdgeColor','r','FaceColor','none','DisplayStyle','stairs',...
        histogram(mean(flex_rv_OF(~~nodeix,:),1),nbins+5,...
            'EdgeColor','k','FaceColor','y',...
            'LineStyle','-.','FaceAlpha',0.5);
            %'EdgeColor','r','FaceColor','none','DisplayStyle','stairs',...
        axis([0.1,0.8,0,16]); 
        xlabel('whole-brain flexibility'); ylabel('number of subjects');
        title('Multi-task Data');
        legend('vision','no vision',...
               'location','northeast');
    
%% T-TESTS of differences in mean...
% (1) between vis and novis for officer
[t1,p1] = ttest(mean(flex_rv_OF(~~nodeix,:),1),mean(flex_v_OF(~~nodeix,:),1));
% (2) between vis and novis for lifespan
[t2,p2] = ttest(mean(flex_rv_LS(~~nodeix,:),1),mean(flex_v_LS(~~nodeix,:),1));
% (3) between vis-novis officer and vis-novis lifespan
[t3,p3] = ttest2( mean(flex_v_OF(~~nodeix,:),1)-mean(flex_rv_OF(~~nodeix,:),1),...
                 mean(flex_v_LS(~~nodeix,:),1)-mean(flex_rv_LS(~~nodeix,:),1));

%% visualize changes in VAR flexibility (figure not in paper)
[sortedflex_LS,sortix_LS] = sort(var(flex_rv_LS(~~nodeix,:)));
[sortedflex_OF,sortix_OF] = sort(var(flex_rv_OF(~~nodeix,:)));
[sortedvflex_LS,sortvix_LS] = sort(var(flex_v_LS(~~nodeix,:)));
[sortedvflex_OF,sortvix_OF] = sort(var(flex_v_OF(~~nodeix,:)));
figure; plot(sortedflex_LS); hold all;
        plot(sortedflex_OF,'r'); hold all;
        plot(sortedvflex_LS,'--b');
        plot(sortedvflex_OF,'--r');
legend('lifespan, no vision','officer, no vision',...  
    'lifespan, vision','officer, vision',...
    'Location','SouthEast');
xlabel('subject ID (sorted)'); ylabel('variance in flexibility');

%% T-TESTS of differences in variance...
% (1) between vis and novis for officer
[t1v,p1v] = ttest(var(flex_rv_OF(~~nodeix,:)),var(flex_v_OF(~~nodeix,:)));
% (2) between vis and novis for lifespan
[t2v,p2v] = ttest(var(flex_rv_LS(~~nodeix,:)),var(flex_v_LS(~~nodeix,:)));
% (3) between vis-novis officer and vis-novis lifespan
[t3v,p3v] = ttest2(var(flex_v_OF(~~nodeix,:))-var(flex_rv_OF(~~nodeix,:)),...
                 var(flex_v_LS(~~nodeix,:))-var(flex_rv_LS(~~nodeix,:)));

%% visualize communities found with and without vision removal (Fig. 6)
x = 1;
ncomms = 11;
figure; subplot(2,2,1);
        bcolor(partnS_v_LS(:,:,x)); 
        colormap([1,1,1;parula(ncomms)]); caxis([0 ncomms]); colorbar;
      subplot(2,2,2);
        bcolor(partnS_v_OF(:,:,x)); 
        colormap([1,1,1;parula(ncomms)]); caxis([0 ncomms]); colorbar;
      subplot(2,2,3);
        bcolor(partnS_rv_LS(:,:,x)); 
        colormap([1,1,1;parula(ncomms)]); caxis([0 ncomms]); colorbar;
      subplot(2,2,4);
        bcolor(partnS_rv_OF(:,:,x)); 
        colormap([1,1,1;parula(ncomms)]); caxis([0 ncomms]); colorbar;

%% same community visualization with separate plots
x = 1;
ncomms_LS = max(max(max(partnS_v_LS(:,:,x))),max(max(partnS_rv_LS(:,:,x))));
figure; subplot(2,1,1);
        bcolor(partnS_v_LS(:,:,x)); 
        colormap([1,1,1;lines(ncomms_LS)]); caxis([0 ncomms_LS]); colorbar;
      subplot(2,1,2);
        bcolor(partnS_rv_LS(:,:,x)); 
        colormap([1,1,1;lines(ncomms_LS)]); caxis([0 ncomms_LS]); colorbar;
ncomms_OF = max(max(max(partnS_v_OF(:,:,x))),max(max(partnS_rv_OF(:,:,x))));        
figure;        
      subplot(2,1,1);
        bcolor(partnS_v_OF(:,:,x)); 
        colormap([1,1,1;lines(ncomms_OF)]); caxis([0 ncomms_OF]); colorbar;
      subplot(2,1,2);
        bcolor(partnS_rv_OF(:,:,x)); 
        colormap([1,1,1;lines(ncomms_OF)]); caxis([0 ncomms_OF]); colorbar;
        
%% visualize 9 random subjects' communities
figure; j=0;
for i=randi(nsubs_LS,[9,1])'
    j = j+1;
    subplot(3,3,j); bcolor(partnS_rv_LS(~~nodeix,:,i)); colorbar; 
end
        