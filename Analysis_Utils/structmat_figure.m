function structmat_figure(Cens,varargin)
% structure subplot.

  if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
  end
  
  plotData = ~varargCheck('AccOnly',varargin{:});
  accessorize = ~varargCheck('DataOnly',varargin{:});
  
  % default blue colormap
  bluemap = @(n_) flipud( ...
                [linspace(0,1,n_)',...
                linspace(0.31,1,n_)',...
                linspace(0.93,1,n_)'] ...
               );
  b100 = bluemap(100); 
  
  %
if plotData
  bcolor(mean(Cens,3)); colormap(b100); caxis([0 1]); hold on;
  if any(ismemvar(varargin,'colorbar'))
      colorbar;
  end
end

if accessorize
  if varargCheck('title',varargin{:})
    title('mean underlying influences');
  end
end