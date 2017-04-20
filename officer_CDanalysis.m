%%%% officerCD analysis

%% load necessary parameters and partition data
taskcell = {'rest';'attn';'face';'word'};
filepath = '/Users/kimberly/Documents/traffic_light/';
%filepath = '/Users/kimberly/Documents/MATLAB/traffic_light/';

% set list of subject_IDs and nsubs
load([filepath 'subject_names.mat']);
subject_IDs = behavior_subjs(ismember(behavior_subjs,hyperedge_subjs));
nsubs = length(subject_IDs);

% set n and t, load all partitions into partnS
partn = load([filepath '../officer_CD/CommAll1.05_' num2str(ib2(1))...
    '.mat'],'CommAll');
[n,t] = size(partn.CommAll);
partnS = zeros(n,t,nsubs);
for k=1:nsubs
    partn = load([filepath '../officer_CD/CommAll1.05_' num2str(ib2(k))...
        '.mat'],'CommAll');
    partnS(:,:,k) = partn.CommAll;
end

% load node positions
pos = load('hybrid_centroids.mat');
pos = squeeze(mean(pos.hybrid_c,2));
pos = pos(ib2,:,:);

fprintf('n = %u, t = %u, %u subjects loaded\n',n,t,nsubs);

%% calculate categorical flexibility for each subject and node

% allocate space in flexS
flexS = zeros(n,nsubs);
for k=1:nsubs
    flexS(:,k) = flexibility(partnS(:,:,k)',1);
end
%figure; bcolor([flexS,mean(flexS,2)]); colorbar;

% create an average flexibility per node figure
BrainView(squeeze(pos(1,:,:))',mean(flexS,2),0,'View','sagittal');
colorbar('Location','WestOutside');
title('Mean Node Flexibility across Subjects');

% create a variance in flexibility per node figure
BrainView(squeeze(pos(1,:,:))',var(flexS,0,2),0,'View','sagittal');
colorbar('Location','WestOutside');
title('Variance in Node Flexibility across Subjects');

% create a scatter plot of means v. vars in node flexibility
figure; scatter(mean(flexS,2),var(flexS,0,2));
xlabel('mean flexibility across subjects'); 
ylabel('variance in flexibility across subjects');
title('mean v. variance of node flexibility across subjs');
legend(['Pearson''s r = ' num2str(corr(mean(flexS,2),var(flexS,0,2)))]);

% create an average flexibility per subject distribution
overavg = mean(mean(flexS));
[nels,cens] = hist(mean(flexS),nsubs);
figure; bar(cens,nels); hold all;
        bar(overavg,max(nels),0.5/nsubs,'FaceColor','r');
xlabel('mean flexibility across nodes'); 
ylabel('number of individual subjects');
title('Mean Flexibility across Nodes: Subject Distribution');
legend(['Subjects (n = ' num2str(nsubs) ')'],'Overall Mean');


meanflex = mean(flexS,2);
varflex = var(flexS,0,2);
figure;scatter(meanflex(meanTb==1),meanTp(meanTb==1),'r');
hold all;scatter(meanflex(meanTb==0),meanTp(meanTb==0),'b');
figure;scatter(varflex(meanTb==1),meanTp(meanTb==1),'r');
hold all;scatter(varflex(meanTb==0),meanTp(meanTb==0),'b');
