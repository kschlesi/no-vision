function [tPosRate,fPosRate] = calculate_single_acc(Fcons,Cens)
% calculates true positive rate and false positive rate for each sim and tw

% Fcons is p x n x T x sims
% Cens is n x n x sims
[~,n,T,sims] = size(Fcons);
[n1,~,sims1] = size(Cens);
assert(n==n1);
assert(sims==sims1);

tPosRate = zeros(T,sims);
fPosRate = zeros(T,sims);
for s=1:sims
  for t=1:T  
    % find association of nodes in one sim/tw
    Cass = mod_allegiance(squeeze(mode(Fcons(:,:,t,s),1)),0);
    Cstr = Cens(:,:,s);
    Cpos = Cass.*~~Cstr;
    Cneg = Cass.*~Cstr;
    tPosRate(t,s) = sum(sum(Cpos))/sum(sum(~~Cstr));
    fPosRate(t,s) = sum(sum(Cneg))/sum(sum(~Cstr));
  end
end

end