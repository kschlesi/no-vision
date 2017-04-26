function var = varargAssign(varString,varDefault,varargin)
% VARARGASSIGN checks the varargin for occurrence of varString. if exists,
% it returns the contents of the varargin cell immediately following the
% cell containing varString. if it does not exist, returns varDefault.
% if varDefault is empty (i.e. [] is passed), returns boolean indicating
% whether varString is in varargin.
% NOTE: requires ismemvar.m, cellmem.m

firstInstance = find(ismemvar(varargin,varString),1,'first');
if isempty(varDefault) % if only checking
    var = any(firstInstance);
else                   % if assigning
    if any(firstInstance)
        var = varargin{firstInstance+1};
    else
        var = varDefault;
    end
end

end