function maskmat = mask_remove(inmat,toRemove)
% zeros out the entries 'toRemove' along first two dimensions of inmat

maskmat = inmat;
maskmat(toRemove,:) = NaN;
maskmat(:,toRemove) = NaN;

end