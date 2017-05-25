function kappaMat = make_kappa_mat(kappa_,N,M,m)
% takes a list of kappa_ values, one for each of the size(M) communities,
% and returns kappaMat, an NxN matrix with kappa_ values in the block
% diagonal entries representing the true underlying communities.

assert(all(size(kappa_)==size(M)));
kappaMat = zeros(N,N);
for sz = 1:numel(M) % size
    for sb = 1:m(sz) % sub-block
        startix = cumsum([0;M.*m]);
        kappaMat(startix(sz)+(sb-1)*M(sz)+1:startix(sz)+sb*M(sz),...
                 startix(sz)+(sb-1)*M(sz)+1:startix(sz)+sb*M(sz)) = kappa_(sz);
    end
end

end