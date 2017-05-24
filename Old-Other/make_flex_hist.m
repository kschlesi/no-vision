% flexS is N x nsubj
X = [[mean(flexS,1)';nan*zeros(size(flexS,1)-size(flexS,2),1)],mean(flexS,2)];
hist3(X,[20 20]);
xlabel('subject flexibility'); ylabel('node flexibility');
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
