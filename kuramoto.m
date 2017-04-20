%%% Simulator for Kuramoto Oscillator Example

sims = 20;
endtime = 100;
ts = 0.1;
period = 2*pi;

% define parameters
kappa_ = 0.2;  % constant coupling strength
N = 128;       % number of nodes
m = 8;         % number of communities
pin = 0.9;     % percentage of in-community nodes connected to node i
pout = 0.1;    % percentage of out-of-community nodes connected to node i


A = zeros(numel(tspan),N,N);
for s=1:sims
    disp(s);

% create initial condition vectors
w = randn(N,1);               % frequencies (normally distributed)
C = modcoupler(N,m,0.8,0.2);  % NxN coupling matrix

theta0 = ones(N,1);
tspan = 0:ts:endtime;

% solve discretised equation & keep results in theta
theta = zeros(numel(tspan),numel(theta0));
for t=1:numel(tspan)
    
  if tspan(t)==0
      % initial condition
      theta(t,:) = theta0;
  else
      sinji = sin(repmat(theta(t-1,:),N,1)-repmat(theta(t-1,:)',1,N));
      theta(t,:) = theta(t-1,:)' + ts.*w + kappa_.*sum(C.*sinji,2);      
  end
   
  % compute synchronization matrix
  cosij = cos(repmat(theta(t,:)',1,N)-repmat(theta(t,:),N,1));
  A(t,:,:) = A(t,:,:) + shiftdim(abs(cosij)./sims,-1);
  
end  % end solve over tspan

end  % end loop over simulations

%%%%%%%% community detection %%%%%%%%%

[C,Q] = genlouvainREPs(squeeze(A(end,:,:)),100,1);
[Cnew,Qnew,Asig,Qcons] = consensus_comm_GL2(C);  % uses "genlouvainREPs"


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

% community detection on FINAL synchronization
figure
h = pcolor(sortentry(C,'row',1)');
set(h, 'EdgeColor', 'none');
colorbar('location','EastOutside');
title('Original');
figure
h = pcolor(sortentry(Cnew,'row',1)');
set(h, 'EdgeColor', 'none');
title('Community Consensus');
colorbar('location','EastOutside');