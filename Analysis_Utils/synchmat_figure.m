function [] = synchmat_figure(Ain,varargin)
% synchronization matrix plot.

  if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
  end
  
  plotData = ~varargCheck('AccOnly',varargin{:});
  accessorize = ~varargCheck('DataOnly',varargin{:});
  
if plotData
  bcolor(Ain); colormap parula; caxis([0 1]); hold on;
  if any(ismemvar(varargin,'colorbar'))
      colorbar;
  end
end

if accessorize
  if varargCheck('title',varargin{:})
    title('example synchronization matrix');
  end
end