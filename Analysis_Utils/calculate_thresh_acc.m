function [tPosMat,fPosMat,tNegMat,fNegMat,precMat] = calculate_thresh_acc(Fass,C,sims)
% calculates tpos rate and fpos rate for DEEP structure at all
% accuracy thresholds (for plotting ROC curves)

% Fass is n x n x T (count of detections)
% C is n x n        (ground truth condition)
[~,n,T] = size(Fass);
[n1,~] = size(C);
Ctime = repmat(C,[1,1,T]);
assert(n==n1);

tPosMat = zeros(sims+1,T);
fPosMat = zeros(sims+1,T);
fNegMat = zeros(sims+1,T);
tNegMat = zeros(sims+1,T);
precMat = zeros(sims+1,T);
for thr=0:sims
    Fthresh = (Fass>=(thr/sims));
    Tpos = Fthresh.*~~Ctime;  % find true positives
    Fpos = Fthresh.*~Ctime;   % find false positives
    Fneg = ~Fthresh.*~~Ctime; % find false negatives
    Tneg = ~Fthresh.*~Ctime;  % find true negatives
    tPosRate = sum(sum(Tpos,1),2)./sum(sum(~~Ctime,1),2); % sensitivity/recall (# true pos / # true)
    fPosRate = sum(sum(Fpos,1),2)./sum(sum(~Ctime,1),2);  % fallout (# false pos / # false)
    fNegRate = sum(sum(Fneg,1),2)./sum(sum(~~Ctime,1),2);  % FNR (# false neg / # true)
    tNegRate = sum(sum(Tneg,1),2)./sum(sum(~Ctime,1),2);  % specificity (# true neg / # false)
    precision = sum(sum(Tpos,1),2)./sum(sum(Fthresh,1),2); % precision
    tPosMat(thr+1,:) = squeeze(tPosRate)';
    fPosMat(thr+1,:) = squeeze(fPosRate)';
    fNegMat(thr+1,:) = squeeze(fNegRate)';
    tNegMat(thr+1,:) = squeeze(tNegRate)';
    precMat(thr+1,:) = squeeze(precision)';
end