function augmat = addback_remove(inmat,N,toRemove,d)
% adds back zeroed filler rows to inmat, assuming added rows were
% originally removed from toRemove indices of the N-element dth dimension 
% of an original matrix

try assert(numel(toRemove)==numel(unique(toRemove)))
catch
    error('toRemove must have no duplicate values');
end

if d==1 % initial support
    augsize = size(inmat);
    augsize(1) = N;
    augmat = zeros(augsize).*NaN;
    notRemoved = removeval(1:N,toRemove);
    augmat(notRemoved,:,:,:) = inmat;
    %for i=1:N
    %    rem_ix = find(notRemoved==i);
    %    if numel(rem_ix)==1
    %        augmat(i,:,:,:) = inmat(rem_ix,:,:,:);
    %    end
    %end
end
      
if d==2 % initial support
    augsize = size(inmat);
    augsize(2) = N;
    augmat = zeros(augsize).*NaN;
    notRemoved = removeval(1:N,toRemove);
    augmat(:,notRemoved,:,:) = inmat;
%     for i=1:N
%         rem_ix = find(notRemoved==i);
%         if numel(rem_ix)==1
%             augmat(:,i,:,:) = inmat(:,rem_ix,:,:);
%         end
%     end
end

%mask_remove(augmat,toRemove);

end