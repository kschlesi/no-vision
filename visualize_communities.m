% FOR A PARTICULAR CHOKING ENVIRONMENT:

goPair = [1.1,0.05];
plot_people = [11,12];
plot_brainruns = [];

goix = find(goRange(goRange(:,1)==goPair(:,1),2)==goPair(:,2));
gamma = goPair(1); omega = goPair(2);
kix = find(ismember(ib,plot_people));
% node positions for choking data
posPass = read_positions('center_of_mass.txt'); % Nx3 positions
posPass = [posPass(:,1),-1*posPass(:,3),posPass(:,2)]; 

for k=1:numel(ib)

C = zeros(n,p*t);
Call = zeros(n,p*t);
Cmode = zeros(n,t);
CmodeZ = zeros(p,n*t);
for T=1:t
    C(:,(T-1)*p+1:T*p) = squeeze(Cplot(goix,(k-1)*p+1:k*p,:,T))';
    Call(:,(T-1)*p+1:T*p) = squeeze(Cplotall(goix,(k-1)*p+1:k*p,:,T))';
    Cmode(:,T) = mode(squeeze(Cplotall(goix,(k-1)*p+1:k*p,:,T))',2);
    CmodeZ(:,(T-1)*n+1:T*n) = squeeze(Cplot(goix,(k-1)*p+1:k*p,:,T));
end
[CmodeZ,zrands] = choose_Z_partition(CmodeZ);
CmodeZ = reshape(CmodeZ,n,t);
if ~all([p,p]==size(zrands)); warning('Z comparison ignores identical partitions'); end; 
clear zrands;
[~,CmodeZ,~] = comm_overlap(CmodeZ,Cmode);

if ismember(k,kix)

% visualize unsorted communities
figure
bcolor(C);
title(['orig. partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(Call);
title(['CC partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(Cmode);
title(['final partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(CmodeZ);
title(['final Z partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);


% visualize sorted communities
% figure
% bcolor(sortentry(C,'col',1));
% title(['SORTED orig. partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(Call,'col',1));
title(['SORTED CC partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(Cmode,'col',1));
title(['SORTED final partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(CmodeZ,'col',1));
title(['SORTED final Z partitions, subject ' num2str(ib(k)) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);

% visualize communities IN THE BRAIN
for T=plot_brainruns
    BrainGraph([], ib(k), zeros(n,1), posPass, Cmode(:,T));
    title(['Partition, subject ' num2str(ib(k)) ', run' num2str(T) ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
    view(90,0); title('Sagittal View',  'FontWeight','bold');
    %  view(180,0); title('Axial View',  'FontWeight','bold');
    %  view(0,90); title('Coronal View',  'FontWeight','bold');
end

end % end visualization over kix people

end

%% plot code (out of date): 2 subjects from 3 different p=100 CC runs
% subject 11 CC
figure;
subplot(1,3,3);bcolor(Cmode3_11);
subplot(1,3,2);bcolor(Cmode2_11);
subplot(1,3,1);bcolor(Cmode1_11);

[~,Cmode2_11n,~] = comm_overlap(Cmode2_11,Cmode1_11); 
[~,Cmode3_11n,~] = comm_overlap(Cmode3_11,Cmode1_11); 

figure;
subplot(1,3,1);bcolor(Cmode1_11);  caxis([1,numel(unique(Cmode3_11n))]);
subplot(1,3,2);bcolor(Cmode2_11n); caxis([1,numel(unique(Cmode3_11n))]);
subplot(1,3,3);bcolor(Cmode3_11n); caxis([1,numel(unique(Cmode3_11n))]);

figure;
subplot(1,3,1);bcolor(sortentry(Cmode1_11,'col',1));
caxis([1,numel(unique(Cmode3_11n))]);
subplot(1,3,2);bcolor(sortentry(Cmode2_11n,'col',0,Cmode1_11(:,1)));
caxis([1,numel(unique(Cmode3_11n))]);
subplot(1,3,3);bcolor(sortentry(Cmode3_11n,'col',0,Cmode1_11(:,1)));
caxis([1,numel(unique(Cmode3_11n))]);

% subject 12 CC
figure;
subplot(1,3,3);bcolor(Cmode3_12);
subplot(1,3,2);bcolor(Cmode2_12);
subplot(1,3,1);bcolor(Cmode1_12);

[~,Cmode2_12n,~] = comm_overlap(Cmode2_12,Cmode1_12); 
[~,Cmode3_12n,~] = comm_overlap(Cmode3_12,Cmode1_12); 

figure;
subplot(1,3,1);bcolor(Cmode1_12);  caxis([1,numel(unique(Cmode3_12n))]);
subplot(1,3,2);bcolor(Cmode2_12n); caxis([1,numel(unique(Cmode3_12n))]);
subplot(1,3,3);bcolor(Cmode3_12n); caxis([1,numel(unique(Cmode3_12n))]);

figure;
subplot(1,3,1);bcolor(sortentry(Cmode1_12,'col',1));
caxis([1,numel(unique(Cmode3_12n))]);
subplot(1,3,2);bcolor(sortentry(Cmode2_12n,'col',0,Cmode1_12(:,1)));
caxis([1,numel(unique(Cmode3_12n))]);
subplot(1,3,3);bcolor(sortentry(Cmode3_12n,'col',0,Cmode1_12(:,1)));
caxis([1,numel(unique(Cmode3_12n))]);

% subject 11 ZSCORE
figure;
subplot(1,3,3);bcolor(CmodeZ3_11);
subplot(1,3,2);bcolor(CmodeZ2_11);
subplot(1,3,1);bcolor(CmodeZ1_11);

[~,CmodeZ2_11n,~] = comm_overlap(CmodeZ2_11,CmodeZ1_11); 
[~,CmodeZ3_11n,~] = comm_overlap(CmodeZ3_11,CmodeZ1_11); 

figure;
subplot(1,3,1);bcolor(CmodeZ1_11);  caxis([1,numel(unique(CmodeZ3_11n))]);
subplot(1,3,2);bcolor(CmodeZ2_11n); caxis([1,numel(unique(CmodeZ3_11n))]);
subplot(1,3,3);bcolor(CmodeZ3_11n); caxis([1,numel(unique(CmodeZ3_11n))]);

figure;
subplot(1,3,1);bcolor(sortentry(CmodeZ1_11,'col',1));
caxis([1,numel(unique(CmodeZ3_11n))]);
subplot(1,3,2);bcolor(sortentry(CmodeZ2_11n,'col',0,CmodeZ1_11(:,1)));
caxis([1,numel(unique(CmodeZ3_11n))]);
subplot(1,3,3);bcolor(sortentry(CmodeZ3_11n,'col',0,CmodeZ1_11(:,1)));
caxis([1,numel(unique(CmodeZ3_11n))]);

% subject 12 ZSCORE
figure;
subplot(1,3,3);bcolor(CmodeZ3_12);
subplot(1,3,2);bcolor(CmodeZ2_12);
subplot(1,3,1);bcolor(CmodeZ1_12);

[~,CmodeZ2_12n,~] = comm_overlap(CmodeZ2_12,CmodeZ1_12); 
[~,CmodeZ3_12n,~] = comm_overlap(CmodeZ3_12,CmodeZ1_12); 

figure;
subplot(1,3,1);bcolor(CmodeZ1_12);  caxis([1,numel(unique(CmodeZ3_12n))]);
subplot(1,3,2);bcolor(CmodeZ2_12n); caxis([1,numel(unique(CmodeZ3_12n))]);
subplot(1,3,3);bcolor(CmodeZ3_12n); caxis([1,numel(unique(CmodeZ3_12n))]);

figure;
subplot(1,3,1);bcolor(sortentry(CmodeZ1_12,'col',1));
caxis([1,numel(unique(CmodeZ3_12n))]);
subplot(1,3,2);bcolor(sortentry(CmodeZ2_12n,'col',0,CmodeZ1_12(:,1)));
caxis([1,numel(unique(CmodeZ3_12n))]);
subplot(1,3,3);bcolor(sortentry(CmodeZ3_12n,'col',0,CmodeZ1_12(:,1)));
caxis([1,numel(unique(CmodeZ3_12n))]);

%%

% turn C into a px(n*t) matrix
Ct1_11 = zeros(p,n*t);
Ct1_12 = zeros(p,n*t);
Ct2_11 = zeros(p,n*t);
Ct2_12 = zeros(p,n*t);
Ct3_11 = zeros(p,n*t);
Ct3_12 = zeros(p,n*t);
for T=1:t
    Ct1_11(:,(T-1)*n+1:T*n) = Corig1_11(:,(T-1)*p+1:T*p)';
    Ct1_12(:,(T-1)*n+1:T*n) = Corig1_12(:,(T-1)*p+1:T*p)';
    Ct2_11(:,(T-1)*n+1:T*n) = Corig2_11(:,(T-1)*p+1:T*p)';
    Ct2_12(:,(T-1)*n+1:T*n) = Corig2_12(:,(T-1)*p+1:T*p)';
    Ct3_11(:,(T-1)*n+1:T*n) = Corig3_11(:,(T-1)*p+1:T*p)';
    Ct3_12(:,(T-1)*n+1:T*n) = Corig3_12(:,(T-1)*p+1:T*p)';
end

% now make a NEW P (pnew=300)x(n*t) matrix
Ct11 = [Ct1_11;Ct2_11;Ct3_11];
Ct12 = [Ct1_12;Ct2_12;Ct3_12];
p=300;

% feed that into community consensus AND z-partitioning
Cmult11 = reshape(consensus_comm_GL2(Ct11),p,n,t);
Cmult12 = reshape(consensus_comm_GL2(Ct12),p,n,t);

Call11 = zeros(t*p,n);
Cmode11 = zeros(n,t);
for T=1:t
    Call11((T-1)*p+1:T*p,:) = Cmult11(:,:,T);
    Cmode11(:,T) = mode(Cmult11(:,:,T),1)';

end
Call12 = zeros(t*p,n);
Cmode12 = zeros(n,t);
for T=1:t
    Call12((T-1)*p+1:T*p,:) = Cmult12(:,:,T);
    Cmode12(:,T) = mode(Cmult12(:,:,T),1)';
end

[CmodeZ11,zrands] = choose_Z_partition(Ct11);
CmodeZ11 = reshape(CmodeZ11,n,t);
if ~all([p,p]==size(zrands)); warning('Z comparison ignores identical partitions'); end; 
clear zrands;
[~,CmodeZ11,~] = comm_overlap(CmodeZ11,Cmode11);

[CmodeZ12,zrands] = choose_Z_partition(Ct12);
CmodeZ12 = reshape(CmodeZ12,n,t);
if ~all([p,p]==size(zrands)); warning('Z comparison ignores identical partitions'); end; 
clear zrands;
[~,CmodeZ12,~] = comm_overlap(CmodeZ12,Cmode12);

figure
bcolor(sortentry(Call11','col',1));
title(['SORTED final partitions, subject 11' ]);%', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(CmodeZ11,'col',1));
title(['SORTED final Z partitions, subject 11' ]);%', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(Call12','col',1));
title(['SORTED final partitions, subject 12']);%  ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
figure
bcolor(sortentry(CmodeZ12,'col',1));
title(['SORTED final Z partitions, subject 12']);%  ', \gamma = ' num2str(gamma) ', \omega = ' num2str(omega)]);
