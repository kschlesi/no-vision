%%% Simulator for Kuramoto Oscillator Example

sims = 100;
endtime = 40;
ts = 0.1;
period = 2*pi;

% define parameters: functional
kappa_ = 0.2;   % constant coupling strength
sigma_ = 1;     % standard deviation of intrinsic frequencies (norm dist)

% define parameters: structural
N = 100;        % number of nodes
M = [20;8];  % unique community sizes present 
m = [2;5];    % number of communities of each size
pin = [0.7,0.6]; % coupling probability for in-community nodes
pout = 0.01;    % coupling probability for out-of-community nodes
%pbase = pout;   % coupling probability for all singletons
pbase = 0.01;

% define parameters: community detection
p = 100;
gamma_ = 1;
omega_ = 0;

%%
% runs:
% M = [20;8]; m = [2;5]; pin = [0.9,0.5]; pout = 0.01; run for 50; works
% M = [10;8]; m = [2;5]; pin = [0.9,0.5]; pout = 0.01; works
% M = [25;15;8]; m = [1;2;4]; pin = [0.7;0.6;0.9]; pout = 0.01; works
% M = [25;15;8]; m = [1;2;4]; pin = [0.7;0.6;0.9]; pout = 0.1; run for 100; no 
% M = [25;15;8]; m = [1;2;4]; pin = [0.7;0.6;0.9]; pout = 0.1; run for 100; no 
% M = [25;15;8]; m = [1;2;4]; pin = [0.7;0.7;0.7]; pout = 0.05; run for 100; 2comms
% M = [25;15;8]; m = [1;2;4]; pin = [0.7;0.7;0.7]; pout = 0.025; run for 100; remove2 works
% M = [1;25;15;8]; m = [3;1;2;4]; pin = [1,0.7,0.7,0.7]; pout = [0.4,0.025,0.025,0.025]; pbase = 0.025; run for 100; remove2 works
%%
% create coupling matrix function:
% takes no arguments, returns NxN binary matrix
Cgenfun = @()modcoupler(N,M,m,pbase,pin,pout,'EnsureConnect');

% figure: single instance of underlying connectivity matrix
figure; bcolor(Cgenfun());

%%
% run sims
[theta_ens,A_ens,C_ens] = ksims_ens(sims,Cgenfun,kappa_,sigma_,ts,endtime);

% figure: mean synchronization over time, within-comm v. global
tot_sync = zeros(size(A_ens,1),length(M)+1);
pstarts = [0;cumsum(M.*m)];
tsync = sum(sum(mean(A_ens,4) - repmat(shiftdim(eye(N),-1),[size(A_ens,1),1,1]),2),3);
for i = 1:length(M)  % iterate thru comm. sizes
  csize = M(i);  
  for j = 1:m(i)     % iterate thru comms of each size
      cstart = pstarts(i) + M(i)*(j-1) + 1;
      cend = pstarts(i) + M(i)*j;
      csync = mean(A_ens(:,cstart:cend,cstart:cend,:),4) - ...
              repmat(shiftdim(eye(csize),-1),[size(A_ens,1),1,1]);
      tot_sync(:,i) = tot_sync(:,i) + mean(mean(csync,2),3);
      tsync = tsync - sum(sum(csync,2),3);
  end
  tot_sync(:,i) = tot_sync(:,i) / m(i);
end
%tsync = mean(A_ens,4) - repmat(shiftdim(eye(N),-1),[size(A_ens,1),1,1]);
%cfrac = M.*(M-1).*m/(N*(N-1));
%wted_sync = bsxfun(@times,tot_sync(:,1:end-1),cfrac');
%tot_sync(:,end) = mean(mean(csync,2),3) - sum(wted_sync(:,1:end-1),2);
tot_sync(:,end) = tsync/(N*N);
figure; plot(tot_sync(1:50,:)); legend([M;0]);

%%
% figure: position v. time (example)
mtheta = mod(theta_ens,period);
figure; colormap jet; bcolor(mtheta(1:50,:,end)'./period); colorbar;

% figure: final synchronization (example)
%figure; bcolor(squeeze(A_ens(end,:,:,end))); colorbar;

% figure: final synchronization (mean)
figure; colormap parula; bcolor(squeeze(mean(A_ens(end,:,:,:),4))); colorbar;

% figure: mean synchronization (over sims and ts)
Ameans = squeeze(mean(mean(A_ens,4),1));
%figure; colormap parula; bcolor(Ameans);

%%
% implement basic community detection
% network of interest: Ameans (mean sync over time and over networks)
[Corig,Qorig] = genlouvainREPs(Ameans,p,gamma_,omega_);
[Ccons,Qcons] = consensus_comm_GL2(Corig);

% figures: CD

% figure: original p community finds v. consensus p community finds
figure; subplot(1,2,1); bcolor(Corig'); colorbar;
        subplot(1,2,2); bcolor(Ccons'); colorbar;
        
% presentation notes: differences for different CD alg runs
% versus differences for different instances of network noise

%%
% do a sweep over gammas
gammas = [1,20,500];

% initialize matrices
Corigs = zeros(p,N,numel(gammas));
Cconss = zeros(p,N,numel(gammas));
for g = 1:numel(gammas)
    [Co,Qo] = genlouvainREPs(Ameans,p,gamma_,omega_);
    [Cc,Qc] = consensus_comm_GL2(Co);
    Corigs(:,:,g) = Co;
    Cconss(:,:,g) = Cc;
end

% figure: plots of success at different gammas
figure; 
for g=1:numel(gammas)
    subplot(numel(gammas),2,2*g-1); bcolor(Corigs(:,:,g)); colorbar;
    subplot(numel(gammas),2,2*g);   bcolor(Cconss(:,:,g)); colorbar;
end

% presentation notes: this can always see the biggest two, and can never
% see the smallest communities, regardless of the value of gamma... (??)

%%
% try removing the first few.
Ameans_RV = Ameans([41:end],[41:end]);
Ameans_RV2 = Ameans([40:end],[40:end]);
[Corig,Qorig1] = genlouvainREPs(Ameans_RV,p,gamma_,omega_);
[Ccons,Qcons1] = consensus_comm_GL2(Corig);

% give the figure...
figure; subplot(1,2,1); bcolor(Corig'); colorbar;
        subplot(1,2,2); bcolor(Ccons'); colorbar;

% presentation notes: IT WORKS
% with 20,20,8,8,8,8,8 AND 10,10,10,10,8,8,8,8,8
% removing first gives different results from removing both.
% now try.... ONE VISION. it works still.

%%
% do a sweep over pins
pinlist = 0:0.2:1;
A_fin = zeros(N,N,numel(pinlist),numel(pinlist));
C_sample = zeros(N,N,numel(pinlist),numel(pinlist));
f=0;
for p1 = pinlist
    for p2 = pinlist
        f = f+1;
        [~,A_ens,C_ens] = ksims(sims,N,M,m,[p1;p2],pout,pbase,kappa_,sigma_,...
                                            ts,endtime);
        A_fin(:,:,pinlist==p1,pinlist==p2) = squeeze(mean(A_ens(end,:,:,:),4));
        C_sample(:,:,pinlist==p1,pinlist==p2) = C_ens(:,:,end) ;
        figure(1); hold on;
        subplot(numel(pinlist),numel(pinlist),f);
        bcolor(A_fin(:,:,pinlist==p1,pinlist==p2)); 
        title(['p1 = ' num2str(p1) ', p2 = ' num2str(p2)]);
        hold off;
        figure(2); hold on;
        subplot(numel(pinlist),numel(pinlist),f);
        bcolor(C_sample(:,:,pinlist==p1,pinlist==p2));
        title(['p1 = ' num2str(p1) ', p2 = ' num2str(p2)]);
        hold off;
    end
end

%%
% do a sweep over community sizes
pinlist = 0.2:0.2:1;
Mlist = [5,10,20,30,50];
siglist = [0.5,1,5,10,50];
As = zeros(numel(pinlist),numel(siglist),numel(Mlist));
Ad = zeros(numel(pinlist),numel(siglist),numel(Mlist));
Cs = zeros(N,N,numel(pinlist),numel(Mlist));
Aend = zeros(N,N,numel(pinlist),numel(siglist),numel(Mlist));
for M = Mlist
  for pin = pinlist
    for sigma_ = siglist
        m = floor(N/M);
        [~,A_ens,C_ens] = ksims(sims,N,M,m,pin,pout,pbase,...
                                kappa_,sigma_,ts,endtime);
        Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M) = squeeze(mean(A_ens(end,:,:,:),4));
        % compute mean total synchronization in each community (assumes all same size)
        As(pinlist==pin,siglist==sigma_,Mlist==M) = ...
            mean(mean(comm_sync(squeeze(A_ens(end,:,:,:)),M,m)));
        % compute sumtotal synchronization between non-community pairs
        Ad(pinlist==pin,siglist==sigma_,Mlist==M) = ...
            sum(sum(squeeze(mean(A_ens(end,:,:,:),4)))) - m*As(pinlist==pin,siglist==sigma_,Mlist==M);
    end
    Cs(:,:,pinlist==pin,Mlist==M) = C_ens(:,:,end);
  end
end

%%
% plot things from previous run
f2 = 0;
for M = Mlist
  f = 0;
  for pin = pinlist
    for sigma_ = siglist
        f = f+1;
        figure(find(Mlist==M)); hold on;
        subplot(numel(pinlist),numel(siglist),f);
        bcolor(Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M));
        title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]); 
        %suptitle(['M = ' num2str(M)]); 
        hold off;
    end
    f2 = f2+1;
    figure(2*numel(Mlist)+1); hold on;
    subplot(numel(pinlist),numel(Mlist),f2);
    bcolor(Cs(:,:,pinlist==pin,Mlist==M)); 
    title(['M = ' num2str(M) ', pin = ' num2str(pin)]); hold off;
  end
  % plot average sync per community, for each M, pin v sig
  figure(numel(Mlist)+find(Mlist==M)); hold on;
  subplot(1,2,1); bcolor(As(:,:,Mlist==M)');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean comm sync, M = ' num2str(M)]); hold off;
  subplot(1,2,2); bcolor(Ad(:,:,Mlist==M)');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean noncomm sync, M = ' num2str(M)]); hold off;
end
%%
%%%%%%%% STATIC community detection %%%%%%%%%
p = 100; gamma = 1; omega = 0;
Mlist = Mlist(1:3);
siglist = siglist(1:3);
C = zeros(p,N,numel(pinlist),numel(siglist),numel(Mlist));
Cnew = zeros(p,N,numel(pinlist),numel(siglist),numel(Mlist));
Q = zeros(p,numel(pinlist),numel(siglist),numel(Mlist));
Qnew = zeros(p,numel(pinlist),numel(siglist),numel(Mlist));
for M=Mlist
  for pin = pinlist
     for sigma_ = siglist
       % static community detection
       [C(:,:,pinlist==pin,siglist==sigma_,Mlist==M),Q(:,pinlist==pin,siglist==sigma_,Mlist==M)] = ...
           genlouvainREPs(Aend(:,:,pinlist==pin,siglist==sigma_,Mlist==M),p,gamma,omega);
       % community consensus
       [Cnew(:,:,pinlist==pin,siglist==sigma_,Mlist==M),Qnew(:,pinlist==pin,siglist==sigma_,Mlist==M)] = ...
           consensus_comm_GL2(C(:,:,pinlist==pin,siglist==sigma_,Mlist==M));  % uses "genlouvainREPs"
     end
  end
end


%%
% community detection results

for M = Mlist
  f = 0;
  f1 = 0;
  for pin = pinlist
    for sigma_ = siglist
      % plot original communities
      f = f+1;
      figure(find(Mlist==M)); hold on;
      subplot(numel(pinlist),numel(siglist),f);
      bcolor(C(:,:,pinlist==pin,siglist==sigma_,Mlist==M)');
      title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]);
      hold off;
      % plot CC communities
      f1 = f1+1;
      figure(find(Mlist==M)+numel(Mlist)); hold on;
      subplot(numel(pinlist),numel(siglist),f1);
      bcolor(Cnew(:,:,pinlist==pin,siglist==sigma_,Mlist==M)');
      title(['M = ' num2str(M) ', pin = ' num2str(pin) ', sig = ' num2str(sigma_)]);
      hold off;
    end
  end
  % make a sig v pin plot of Q's
  figure(2*numel(Mlist)+find(Mlist==M)); hold on;
  bcolor(squeeze(mean(Q(:,:,:,Mlist==M),1))');
  xlabel('in-comm connection prob'); ylabel('natural freq variance');
  set(gca,'XTick',(0:1:numel(pinlist)-1)+1.5,'YTick',(0:1:numel(siglist)-1)+1.5);
  set(gca,'XTickLabel',pinlist,'YTickLabel',siglist);
  title(['mean Q, M = ' num2str(M)]); colorbar; hold off;
end

%%
%%%%%%%% plots %%%%%%%%%

% from single realisation: phase of each oscillator over time
figure
plot(tspan,theta(:,1))
hold on
hold all
for i=2:N
    plot(tspan,theta(:,i))
end

% single realisation: final phase of each oscillator
mtheta = mod(theta,period);
figure
plot(1:N,sort(mtheta(end,:))./period)

% synchronization of osc 1 with all others over time (avg of sims realzns)
figure
plot(tspan,A(:,:,1))

