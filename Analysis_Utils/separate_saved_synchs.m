% script for separately saving each A type thing...

%% first, test the loading time with and without A_ens

paramString = 's8nrun_g1o1';

tic;
load(['Results/' paramString '.mat']);
toc;

clearvars -except paramString

tic;
load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens).)*$');
toc;

clearvars -except paramString

tic;
load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens|Cens|Rass|Rg).)*$');
toc;
