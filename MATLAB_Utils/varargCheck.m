function var = varargCheck(varString,varargin)
% checks whether varString is contained in a cell in varargin

var = varargAssign(varString,[],varargin{:});

end