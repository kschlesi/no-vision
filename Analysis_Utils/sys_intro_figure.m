function [] = sys_intro_figure(Cens,Rass,varargin)
% does things for sys_intro fig. 

if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
end

if varargCheck('title',varargin{:})
    title(varargAssign('title',' ',varargin{:}));
end

% first subplot
subplot(1,2,1);
structmat_figure(Cens,'ExtFigureCmd',varargin{:});

% second subplot
%subplot(1,3,2);
%synchmat_figure(,'ExtFigureCmd',varargin{:});

% third subplot
subplot(1,2,2);
det_prob_tAvg_figure(Rass,'ExtFigureCmd','colorbar',varargin{:});