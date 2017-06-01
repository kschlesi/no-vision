function [] = tpfp_figure(Rg,Cens,varargin)
% creates true pos / false pos figure and accessories from input data

  if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
  end
  
  plotData = ~varargCheck('AccOnly',varargin{:});
  accessorize = ~varargCheck('DataOnly',varargin{:});

  if plotData
      [tPos,fPos] = calculate_single_acc(Rg,Cens);
      if varargCheck('xplot',varargin{:})
        scatter(mean(fPos),mean(tPos),'x','MarkerEdgeColor',last_color);
      else
        scatter(mean(fPos),mean(tPos),'o','MarkerFaceColor',next_color);
      end
  end
  
  if accessorize
      xlabel('false positive rate');
      ylabel('true positive rate');
    if varargCheck('legend',varargin{:})  
      legend(varargAssign('legend',varargin{:}),'location','southeast');
    end
    if varargCheck('title',varargin{:})
      title('CD performance by run');
    end
  end
  
end