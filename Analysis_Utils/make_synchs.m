function A = make_synchs(A_ens,varargin)
% takes A_ens & makes synch matrices with (optional) T, tLength, transient

% set size variables
[nts,n,n1,sims] = size(A_ens);
assert(n==n1);

% set input variables
T = varargAssign('T',1,varargin{:});
transient = varargAssign('transient',0,varargin{:});
tLength = varargAssign('tLength',floor((nts-transient)/T),varargin{:});
if varargCheck('tLength',varargin{:}) && ~varargCheck('T',varargin{:})
    T = floor((numel(0:ts:endtime)-transient)/tLength);
end

% make synchs
A = cell(T,sims);
for s=1:sims
    for t=1:T
        A{t,s} = squeeze(...
                    mean(...
                      A_ens(1+transient+(t-1)*tLength:t*tLength+transient,...
                            :,:,s)...
                    ,1)...
                 );
    end
end

end