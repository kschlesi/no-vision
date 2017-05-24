function [tPosMat,fPosMat] = calculate_thresh_acc(Fass,C,sims)
% calculates tpos rate and fpos rate for DEEP structure at all
% accuracy thresholds (for plotting ROC curves)

% Fass is n x n x T
% C is n x n
[~,n,T] = size(Fass);
[n1,~] = size(C);
Ctime = repmat(C,[1,1,T]);
assert(n==n1);

tPosMat = zeros(sims+1,T);
fPosMat = zeros(sims+1,T);
for thr=0:sims
    Fthresh = (Fass>=(thr/sims));
    Fpos = Fthresh.*~~Ctime;
    Fneg = Fthresh.*~Ctime;
    tPosRate = sum(sum(Fpos,1),2)./sum(sum(~~Ctime,1),2);
    fPosRate = sum(sum(Fneg,1),2)./sum(sum(~Ctime,1),2);
    tPosMat(thr+1,:) = squeeze(tPosRate)';
    fPosMat(thr+1,:) = squeeze(fPosRate)';
end