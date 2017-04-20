function adj = modcoupler(N,M,m,pbase,pin,pout,varargin)
% creates adjacency matrix (binary) for a network with an arbitrary number
% of communities and diffeering community sizes and coupling probabilities.
% N = total number of nodes
% M = vector of length nS, containing the nS unique community sizes
% m = vector of length nS, containing the no. of communities of each size
%       nC = sum(m)  (note: m must be the same size as M OR scalar)
% pbase = scalar containing base coupling probability for non-comm nodes
% pin = vector containing the in-comm coupling probability for each comm.
%       size (if numel(pin)==nS) or for each community (if numel(pin)==nC)
% pout = vector containing the out-comm coupling probability for each comm.
%       size (if numel(pout)==nS) or for each community (if numel(pout)==nC)
% Note: pin and/or pout may be scalar (same for all nodes)
% If optional argument 'EnsureConnect' is included, followed by a max number 
% of iterations, matrices will be generated until no disconnected in-
% community nodes are created. Default number of iterations = 1000. 

M = M(:);
m = m(:);
pin = pin(:);
pout = pout(:);
nS = numel(M);
if isscalar(m); m = repmat(m,nS,1);
else assert(nS == numel(m)); end;
if isscalar(pin); pin = repmat(pin,nS,1);
else assert(nS == numel(pin)); end;
if isscalar(pout); pout = repmat(pout,nS,1);
else assert(nS == numel(pout)); 
     warning('non-matching pout values; first entered is default'); 
end;
%nC = sum(m);

adj = zeros(N,N);
piconn = M.*(M-1)./2;
niconn = floor(pin.*piconn);
urows = [0;cumsum(M.*m)];
for i = 1:nS
  for c = 1:m(i)
    % select within-community connections
      cblock = zeros(M(i));
      cnx = zeros(piconn(i),1);
      cnx(randperm(piconn(i),niconn(i))) = 1;
      cblock(~~triu(ones(M(i)),1)) = cnx;
    % select between-community connections
      srow = urows(i)+(c-1)*M(i)+1;
      erow = urows(i)+c*M(i);
      poconn = (N-erow)*M(i);
      cnx = zeros(poconn,1);
      cnx(randperm(poconn,floor(poconn*pout(i)))) = 1;
      oblock = reshape(cnx,[M(i),N-erow]);
    % fill in adj
      adj(srow:erow,srow:erow) = cblock;
      adj(srow:erow,erow+1:end) = oblock;
  end  % end loop over comms within class
end  % end loop over size classes

% connect singletons @ rate pbase
singles = N-urows(end);
sblock = zeros(singles);
psconn = singles*(singles-1)/2;
cnx = zeros(psconn,1);
cnx(randperm(psconn,floor(psconn*pbase))) = 1;
sblock(~~triu(ones(singles),1)) = cnx;
adj(urows(end)+1:end,urows(end)+1:end) = sblock;

% symmetrize matrix, keeping all connections
adj = adj + adj';

% ensure no disconnected in-community nodes (verbose)
if(any(ismemvar(varargin,'EnsureConnect')))
    try
      maxiter = varargin{find(ismemvar(varargin,'EnsureConnect'),1,'first')+1};
      assert(isnumeric(maxiter) && isscalar(maxiter));
    catch
      maxiter = 1000;
    end
    
    i=0;
    while any(sum(adj(1:M'*m,1:M'*m))==0) && i<maxiter
        disp('ensuring no disconnected in-community nodes...');
        disp(i); i=i+1;
        adj = modcoupler(N,M,m,pbase,pin,pout);  % NxN binary coupling matrix
    end
end

end