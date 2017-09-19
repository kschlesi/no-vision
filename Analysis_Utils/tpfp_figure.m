function [tps,fps] = tpfp_figure(Rg,Cens,varargin)
% creates true pos / false pos figure and accessories from input data

  if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
  end
  
  addEllipse = varargAssign('addEllipse',false,varargin{:});
  
  plotData = ~varargCheck('AccOnly',varargin{:});
  accessorize = ~varargCheck('DataOnly',varargin{:});

  if plotData
      [tPos,fPos] = calculate_single_acc(Rg,Cens);
      plotargs = varargAssign('plotargs',{'o'},varargin{:});
      tps = mean(tPos);
      fps = mean(fPos);
      if addEllipse
          ellipseArgs = varargAssign('ellipseArgs',{'NoLine'},varargin{:});
          filled_ellipse(std(fps),std(tps),0,mean(fps),mean(tps),ellipseArgs{:});
      end
      plot(fps,tps,plotargs{:});
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