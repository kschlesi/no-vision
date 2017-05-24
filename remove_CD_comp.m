function [Rassoc,Fassoc,R,F,Rcons,Fcons,A_ens,C_ens] = ...
               remove_CD_comp(genFun,toRemove,sims,varargin)
% This function completes an experiment on an input generating function for
% a Kuramoto network structure, given the varargins for ksims_ens as well
% as subset of nodes to remove. it runs the experiment on both the full
% network and the nodes-removed network (default) and compares the results.
% It returns the simulation-averaged CD results from both networks, as well
% as (optionally) making comparison figures for them.
% *****
% Inputs. genFun   : function handle to generator of structure. (no args)
%         toRemove : vector of indices of nodes to remove.
%         sims     : number of simulations in ensemble.
% Optional inputs (precede with input name string).
%         kappa_   : (optional) arg passed to ksims (default = 0.2)
%         sigma_   : (optional) arg passed to ksims (default = 1)
%         ts       : (optional) arg passed to ksims (default = 0.1)
%         endtime  : (optional) arg passed to ksims (default = 50)
%         gamma_   : (optional) arg passed to genlouvain (default = 1)
%         omega_   : (optional) arg passed to genlouvain (default = 0)
%         p        : (optional) arg passed to genlouvain (default = 100)
%         T        : (optional) number of time windows (default = 1)
%         tLength  : (optional)
%         transient: (optional) number initial ts to ignore (default = 0)
%         makePlot : (optional) bool, true makes figs (default = false)
%         doFullRun: (optional) bool, true runs on full network (default = true)

% assign varargs
kappa_ = varargAssign('kappa_',0.2,varargin{:});
sigma_ = varargAssign('sigma_',1,varargin{:});
ts = varargAssign('ts',0.1,varargin{:});
endtime = varargAssign('endtime',50,varargin{:});
gamma_ = varargAssign('gamma_',1,varargin{:});
omega_ = varargAssign('omega_',0,varargin{:});
p = varargAssign('p',100,varargin{:});
T = varargAssign('T',1,varargin{:});
transient = varargAssign('transient',0,varargin{:});
tLength = varargAssign('tLength',floor((numel(0:ts:endtime)-transient)/T),varargin{:});
makePlot = varargAssign('makePlot',false,varargin{:});
doFullRun = varargAssign('doFullRun',true,varargin{:});

% if tLength was specified but T not specified, divide into equal windows
if varargAssign('tLength',[],varargin{:}) && ~varargAssign('T',[],varargin{:})
    T = floor((numel(0:ts:endtime)-transient)/tLength);
end

% run kuramoto sims
[~,A_ens,C_ens] = ksims_ens(sims,genFun,kappa_,sigma_,ts,endtime);
n = size(A_ens,2);

% initialize vars
nR = n-numel(toRemove);
if doFullRun
    F = zeros(p,n,T,sims);
    R = zeros(p,nR,T,sims);
else
    F = []; R = [];
end
Fcons = zeros(p,n,T,sims);
Rcons = zeros(p,nR,T,sims);

% do all experiments
for s=1:sims
  if doFullRun
  % run initial full CD experiment    
    A = cell(T,1);
    for t=1:T
        A{t} = squeeze(mean(A_ens(1+transient+(t-1)*tLength:t*tLength+transient,...
                     :,:,s),1));
    end
    C = genlouvainREPs(A,p,gamma_,omega_);
    Ccons = zeros(p,n*T);
    for i=1:T
        Ccons(:,:,i) = consensus_comm_GL2(C(:,:,i));
    end
    %Ccons = reshape(consensus_comm_GL2(Ccons),[p,n,T]);
    F(:,:,:,s) = C;
    Fcons(:,:,:,s) = Ccons;
  end
  
  % run removal CD experiment
  A = cell(T,1);
  for t=1:T
      A{t} = squeeze(mean(A_ens(1+transient+(t-1)*tLength:t*tLength+transient,...
                   removeval(1:n,toRemove),removeval(1:n,toRemove),s),1));
  end
  C_r = genlouvainREPs(A,p,gamma_,omega_);
  Ccons_r = zeros(p,nR,T);
  for i=1:T
      Ccons_r(:,:,i) = consensus_comm_GL2(C_r(:,:,i));
  end
  R(:,:,:,s) = C_r;
  Rcons(:,:,:,s) = Ccons_r;
end

% compute statistical success metrics...
if doFullRun
    Fassoc = zeros(n,n,T);
else
    Fassoc = [];
end
Rassoc = zeros(nR,nR,T);
for t=1:T
  if doFullRun
    Fassoc(:,:,t) = mod_allegiance(squeeze(mode(Fcons(:,:,t,:),1))',0);
  end
    Rassoc(:,:,t) = mod_allegiance(squeeze(mode(Rcons(:,:,t,:),1))',0);
end

if makePlot
    
  % blue colormap
  bluemap = @(n_) flipud( ...
                [linspace(0,1,n_)',...
                linspace(0.31,1,n_)',...
                linspace(0.93,1,n_)'] ...
               );
  b100 = bluemap(100);         
  
  % check of consensus
  for s=1:3
    figure;
    if doFullRun
      ncommsF = numel(unique(F(:,:,:,1)));
    end
      ncommsR = numel(unique(R(:,:,:,1)));
        for i=1:T
          if doFullRun  
            subplot(T,2,i*2-1);
            bcolor(F(:,:,i,s)'); 
            colormap(lines(ncommsF)); caxis([1 ncommsF]); colorbar;
            subplot(T,2,i*2);
            bcolor(Fcons(:,:,i,s)'); 
            colormap(lines(ncommsF)); caxis([1 ncommsF]); colorbar;
          else
            subplot(T,2,i*2-1);
            bcolor(R(:,:,i,s)'); 
            colormap(lines(ncommsR)); caxis([1 ncommsR]); colorbar;
            subplot(T,2,i*2);
            bcolor(Rcons(:,:,i,s)'); 
            colormap(lines(ncommsR)); caxis([1 ncommsR]); colorbar;  
          end
        end
  end
  
  % look over time windows... for consensus communities
  for s=1:5
      % create nx(p*T) matrix
    if doFullRun  
      pFmat = zeros(n,p*T);
    end
      pRmat = zeros(nR,p*T);
      for t=1:T
        if doFullRun  
          pFmat(:,(t-1)*p+1:t*p) = Fcons(:,:,t,s)';
        end
          pRmat(:,(t-1)*p+1:t*p) = Rcons(:,:,t,s)';
      end
    if doFullRun  
      ncommsF = numel(unique(pFmat));
      figure; bcolor(pFmat); 
              colormap(lines(ncommsF)); caxis([1 ncommsF]); colorbar;
    end
      ncommsR = numel(unique(pRmat));
      figure; bcolor(pRmat); 
              colormap(lines(ncommsR)); caxis([1 ncommsR]); colorbar;        
  end
  
  % figure for full network
  if doFullRun
    %figF = nextFig;
    %figure(figF);
    figure;
      subplot(1,T+1,1);
        bcolor(mean(C_ens,3)); colormap(b100); caxis([0 1]); hold on;
      for t=1:T
        subplot(1,T+1,t+1);
          bcolor(Fassoc(:,:,t)); colormap(b100); caxis([0 1]);
      end
      hold off;
  end
      
  % figure for removed nodes
    %figR = nextFig;
    %figure(figR);
    figure;
      subplot(1,T+1,1);
        bcolor(mean(C_ens(removeval(1:n,toRemove),removeval(1:n,toRemove),:),3));...
        colormap(b100); caxis([0 1]); hold on;
      for t=1:T
        subplot(1,T+1,t+1);
          bcolor(Rassoc(:,:,t)); colormap(b100); caxis([0 1]);
      end
      hold off;  
      
end
      
end





