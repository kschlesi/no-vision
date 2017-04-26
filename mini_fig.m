% mini network script for figure

%% add path vars
addpath(genpath('Kuramoto/'));
addpath(genpath('MATLAB_Utils/'));
addpath(genpath('CD_Utils/'));

%% define connections
% module 1 (nodes 1-4)
mod1 = [0,1,1,1;
        1,0,1,1;
        1,1,0,1;
        1,1,1,0];
    
% module 2 (nodes 5-9)
mod2 = [0,1,0,0,1;
        1,0,1,1,1;
        0,1,0,1,1;
        0,1,1,0,1;
        1,1,1,1,0];

% module 3 (nodes 10-15)
mod3 = [0,1,0,0,1,1;
        1,0,1,1,1,1;
        0,1,0,1,1,1;
        0,1,1,0,1,0;
        1,1,1,1,0,1;
        1,1,1,0,1,0];

% module 4 (nodes 16-21)
mod4 = [0,1,0,1,1,1;
        1,0,1,1,0,1;
        0,1,0,1,1,0;
        1,1,1,0,1,0;
        1,0,1,1,0,1;
        1,1,0,0,1,0];
    
% module 5 (nodes 22-26)
mod5 = [0,1,1,1,1;
        1,0,1,1,1;
        1,1,0,1,0;
        1,1,1,0,1;
        1,1,0,1,0];

% other connections
others = [2,9;
          8,10;
          14,18;
          15,23;
          3,22;
          27,1;
          27,6;
          28,13;
          28,19;
          29,26;
          30,2;
          30,9;
          30,15;
          30,23];

%% create adjacency matrix

adj = blkdiag(mod1,mod2,mod3,mod4,mod5,zeros(4));
for i=1:length(others)
    adj(others(i,1),others(i,2)) = 1;
    adj(others(i,2),others(i,1)) = 1;
end
figure; bcolor(adj); colormap([1,1,1;0,0.31,0.93]);

%% run one simulation

kappa = 0.15;   % constant coupling strength
sigma = 3;     % standard deviation of intrinsic frequencies (norm dist)
[theta,A] = ksims(adj,kappa,sigma,0.1,50);
theta = mod(theta,2*pi);
n = size(theta,2);
%figure; bcolor(cos(theta)); colormap(jet); colorbar;
h=figure; plot( cos(theta) + ...
              repmat(3.5*(1:n),[size(theta,1),1]),'k');
        ylabel('oscillator');
        xlabel('time step');
        h.CurrentAxes.XLim = [0,500];
        h.CurrentAxes.YLim = [0,110];
        h.CurrentAxes.YTickLabel = 5:5:n; 
        h.CurrentAxes.YTick = 3.5*(5:5:n);
%%
bluemap = @(n_) flipud( ...
                [linspace(0,1,n_)',...
                linspace(0.31,1,n_)',...
                linspace(0.93,1,n_)'] ...
               );
figure; bcolor(squeeze(A(end,:,:))); colormap(bluemap(100)); colorbar;

T = 3;
fcell = cell(T,1);
for i=1:T
    fcell{i} = squeeze(mean(A(((i-1)*150)+1:i*150,:,:)));
    figure; bcolor(fcell{i});
            colormap(bluemap(100)); caxis([0 1]); colorbar;
end

%% apply CD and get dynamics

p = 100;
C = genlouvainREPs(fcell,p,1,0.01);
Ccons = zeros(p,n,T);
for i=1:T
    Ccons(:,:,i) = consensus_comm_GL2(C(:,:,i));
end

figure;
for i=1:T
    subplot(T,2,i*2-1);
    bcolor(C(:,:,i)'); 
    colormap(lines(5)); caxis([1 5]); colorbar;
    subplot(T,2,i*2);
    bcolor(Ccons(:,:,i)'); 
    colormap(lines(5)); caxis([1 5]); colorbar;
end

figure; bcolor(squeeze(Ccons(1,:,:))); 
        colormap(lines(5)); caxis([1 5]); colorbar;
        disp(colormap);

%% remove nodes 1-4 and apply CD again

ixRemove = [22:26,29];
fcell_r = cell(T,1);
for i=1:T
    fcell_r{i} = fcell{i}(removeval(1:n,ixRemove),removeval(1:n,ixRemove));
end

C_r = genlouvainREPs(fcell_r,p,1,0.01);
Ccons_r = zeros(p,n-numel(ixRemove),T);
for i=1:T
    Ccons_r(:,:,i) = consensus_comm_GL2(C_r(:,:,i));
end

ncomms = numel(unique(C_r));
figure;
for i=1:T
    subplot(T,2,i*2-1);
    bcolor(C_r(:,:,i)'); 
    colormap(lines(ncomms)); caxis([1 ncomms]); colorbar;
    subplot(T,2,i*2);
    bcolor(Ccons_r(:,:,i)'); 
    colormap(lines(ncomms)); caxis([1 ncomms]); colorbar;
end

figure; bcolor(squeeze(Ccons_r(1,:,:)),(1:T)',removeval(1:n,ixRemove)); 
        colormap(lines(ncomms)); caxis([1 ncomms]); colorbar;
        disp(colormap);