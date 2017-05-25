function [theta,A] = ksims(C,kappa_,sigma_,ts,endtime)
%%%%%%% run simulation %%%%%%%%

% read in structure matrix
try
    assert(ismatrix(C));
    [N,N1] = size(C);
    assert(N1==N);
catch
    error('Input structure matrix must be 2D and square.');
end

% initialize simulation time
tspan = 0:ts:endtime;
nsteps = numel(tspan);

% if kappa_ is kappa_ij (i.e. size(C)), use it as-is.
if all([N,N1]==size(kappa_))
    tDepKappa = false;
else % else, create time-dependent kappa_
    kappa_ = create_tdep_kappa(kappa_,nsteps);
    tDepKappa = true;
end

% create initial condition vectors
w = randn(N,1).*sigma_;                % frequencies (normally distributed)

% initialize state variables
theta0 = unifrnd(0,2*pi,[N,1]);
A = zeros(nsteps,N,N);

% solve discretised equation & keep results in theta
theta = zeros(nsteps,numel(theta0));
theta(1,:) = theta0;
for t=2:nsteps
  
  % evolve oscillators forward in time
  sinji = sin(repmat(theta(t-1,:),N,1)-repmat(theta(t-1,:)',1,N));
  if tDepKappa
      theta(t,:) = theta(t-1,:)' + ts.*w + kappa_(t-1).*sum(C.*sinji,2);      
  else
      theta(t,:) = theta(t-1,:)' + ts.*w + sum(kappa_.*C.*sinji,2);    
  end
   
  % compute synchronization matrix
  cosij = cos(repmat(theta(t,:)',1,N)-repmat(theta(t,:),N,1));
  %A(t,:,:) = A(t,:,:) + shiftdim(abs(cosij)./sims,-1);
  A(t,:,:) = shiftdim(abs(cosij),-1);
  
end  % end solve over tspan

end

function outKappa = create_tdep_kappa(inKappa,nsteps)
% takes an input kappa_ (numeric scalar, vector, or nD matrix) and converts
% it to an output of length (nsteps-1). this it does by truncating (if 
% numel(kappa) > nsteps-1) or by expanding. example: if there are N
% elements in kappa_, kappa_(1) covers the first Nth of steps, kappa_(2)
% the 2nd Nth, ...., and kappa_(end) covers the last full Nth division of
% steps PLUS the remaining steps.

if isscalar(inKappa)
    outKappa = repmat(inKappa,nsteps-1,1);
else
    inKappa = inKappa(:);
    nK = numel(inKappa);
    if nK >= nsteps-1
        outKappa = inKappa(1:nsteps-1);
    else
        kdivs = floor((nsteps-1)/nK);
        krem = mod(nsteps-1,nK);
        outKappa = zeros(nsteps-1,1);
        for i=1:nK
            outKappa((i-1)*kdivs+1:i*kdivs) = inKappa(i);
        end
        if ~~krem  % if remainder
            outKappa(nK*kdivs+1:end) = inKappa(end);
        end
    end
end
assert(numel(outKappa) == nsteps-1);

end
