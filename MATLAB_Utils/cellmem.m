function comp = cellmem(A,B)
% this takes two cell arrays, A and B. returns a logical vector C the same
% size as A, with C(i) = 1 if the contents of A(i) are identical to the
% contents of ANY cell in B.

comp = zeros(size(A));
for a=1:numel(A)
    xinb = zeros(size(B));
    for b=1:numel(B)
        % if either cell (A or B) CONTAINS a char, use strcmp;
        % else use all(equals)
        if ischar(A{a}) || ischar(B{b})
%             try
              xinb(b) = all(all(strcmp(A{a},B{b}))) && numel(A{a})==numel(B{b});
%             catch oops
%               if oops
%                 xinb(b) = all(reshape(strcmp(A{a},B{b}),numel(strcmp(A{a},B{b}),1))) ...
%                       && numel(A{a})==numel(B{b});  
%               end
%             end
        else
            try AeqB = A{a}==B{b};
            catch oops
                if oops; AeqB = 0; end;
            end
            xinb(b) = (all(AeqB(:)) && numel(A{a})==numel(B{b}));
        end
    end
    comp(a) = any(xinb(:));
end
        
end