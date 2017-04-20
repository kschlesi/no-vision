% visualize officer communities

n = 194; % number of nodes (brain regions)
p = 100; % number of genlouvain optimizations to perform
t = 4;   % number of slices in multislice (do not change without modifying
                          % code below; this code is set up for rest, attn,
                          % word, face (t=4) by the switch statements below
tosave = 1;  % save result? 1=yes, 0=no
file_path = '/Users/kimberly/Google Drive/traffic_light/'; % string giving 
                                    % location to find saved node node 
                                    % matrices & save output partitions

% list of gamma,omega pairs to test (omega~=0 gives categorical multislice)
goRange = [1,0.14];

plot_brainruns = [1,2];

load('missing.mat'); % binary vectors indicating subjs w full data per task
% ib = vector containing indices of subjects with full data for all tasks
i1 = find(att1_miss_adj>0);
i2 = find(att2_miss_adj>0);
i3 = find(face_miss_adj>0);
i4 = find(word_miss_adj>0);
i5 = find(rest_miss_adj>0);
ib = intersect(intersect(intersect(intersect(i1,i2),i3),i4),i5);

% override ib (if desired) to restrict analysis to only a few subjects
ib = 68;
disp(ib);

% remove missing regions
toremove = 0;
missing = [90:96,188:194];
if toremove
    str1='nosub';
    oldN = n;
    n = n-length(missing);
else
    str1=[];
end

for goix=1:size(goRange,1)
    gamma = goRange(goix,1);
    omega = goRange(goix,2);

for k=1:numel(ib)
    person = ib(k);

% read in data
C = csvread([file_path 'subj' num2str(person) '_' num2str(gamma) '_' num2str(omega) str1 '_Call.txt']);
Cnew = csvread([file_path 'subj' num2str(person) '_' num2str(gamma) '_' num2str(omega) str1 '_Cnewall.txt']);
Cmode = zeros(n,t);
for T=1:t
    Cmode(:,T) = mode(Cnew((T-1)*p+1:T*p,:)',2);
end

figure
bcolor(C');
title(['orig. partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(Cnew');
title(['CC partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);

% positions of regions
pos = load('hybrid_centroids.mat');
posPass = squeeze(mean(pos.hybrid_c(ib(k),:,:,:),2))'; % Nx3 positions
if toremove
    newix = removeval(1:oldN,missing)';
    posPass = posPass(newix,:);
end
for T=plot_brainruns
    BrainGraph([], ib(k), zeros(n,1), posPass, Cmode(:,T));
    title(['Partition, subject ' num2str(ib(k)) ', run' num2str(T) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
    view(90,0); title('Sagittal View',  'FontWeight','bold');
    %  view(180,0); title('Axial View',  'FontWeight','bold');
    %  view(0,90); title('Coronal View',  'FontWeight','bold');
end

if toremove
    newix = removeval(1:oldN,missing)';
    C2 = add_missings_dim(Cmode,newix,oldN,1,0);
else
    C1 = Cmode;
end   

end

end


%% run this after running above for toremove=1 and toremove=0

[overlap,C1new] = comm_overlap(C1,C2);

xlabels={...%'orig subc',
    'matched subc','orig no-subc','overlap'};
figure;
for T=1:t;
    subplot(2,2,T);
    [mat,theix]=sortentry([...%C1(newix,T),
        C1new(newix,T),C2(newix,T),overlap(newix,T)],'col',...%4
        3);
    bcolor(mat);
    title(['task' num2str(T)]);
    ylabel('brain region ID');
    set(gca,'XTick',1.5:1:...%4.5
        3.5,'XTickLabel',xlabels,'YTick',1:20:length(newix),'YTickLabel',theix(1:20:length(newix)));
end;

A = cell(t,1);
Asmall = cell(t,1);
    for T=1:t
        switch T  % task info for reading in node-node matrices
            case 1, task = 'rest';
            case 2, task = 'attn';
            case 3, task = 'word';
            case 4, task = 'face';
        end
        A{T} = csvread([file_path 'N' task '_' num2str(person) '.dat']);
        Asmall{T} = remove_missings(A{T},missing);
    end
newix = removeval(1:oldN,missing)';
C2small = C2(newix,:);
pcoefNS = add_missings_dim(partcoef(C2small,Asmall),newix,oldN,1,0);
pcoefWS = partcoef(C2,A);

figure;
plot1 = mean(pcoefWS(newix,:),2)-mean(mean(pcoefWS(newix,:)));
plot2 = mean(pcoefNS(newix,:),2)-mean(mean(pcoefNS(newix,:)));
[plot1s,six1] = sort(plot1); 
[plot2s,six2] = sort(plot2);
%subplot(2,1,1)
plot(plot1)
hold all;
plot(plot2,'--');
legend(['with subcort. included, top five: ' num2str(six1(end)) ', ' num2str(six1(end-1)) ', ' num2str(six1(end-2)) ', ' num2str(six1(end-3)) ', ' num2str(six1(end-4)) ],...
        ['subcort. NOT included, top five: ' num2str(six2(end)) ', ' num2str(six2(end-1)) ', ' num2str(six2(end-2)) ', ' num2str(six2(end-3)) ', ' num2str(six2(end-4)) ],...
        'Location','SouthEast');
% subplot(2,1,2)
% plot(plot1s);
% hold all;
% plot(plot2(six),'--');

figure;
title('participation coefficient with and without subcortical regions');
for T=1:t
    switch T  % task info for reading in node-node matrices
            case 1, task = 'rest';
            case 2, task = 'attn';
            case 3, task = 'word';
            case 4, task = 'face';
    end
    subplot(2,2,T);
    plot(pcoefWS(newix,T)-mean(pcoefWS(newix,T)));
    hold all;
    title(task);
    plot(pcoefNS(newix,T)-mean(pcoefNS(newix,T)),'--');
    [~,ix1] = sort(pcoefWS(newix,T));
    [~,ix2] = sort(pcoefNS(newix,T));
    legend(['with subcort. included, top five: ' num2str(ix1(end)) ', ' num2str(ix1(end-1)) ', ' num2str(ix1(end-2)) ', ' num2str(ix1(end-3)) ', ' num2str(ix1(end-4)) ],...
        ['subcort. NOT included, top five: ' num2str(ix2(end)) ', ' num2str(ix2(end-1)) ', ' num2str(ix2(end-2)) ', ' num2str(ix2(end-3)) ', ' num2str(ix2(end-4)) ],...
        'Location','SouthEast');
end





