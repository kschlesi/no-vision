function [] = ROC_figure(Rass,N,M,m,toRemove,sims,varargin)
% creates true pos / false pos figure and accessories from input data

  if ~varargCheck('ExtFigureCmd',varargin{:})
    figure(next_fig);
  end
  
  plotData = ~varargCheck('AccOnly',varargin{:});
  accessorize = ~varargCheck('DataOnly',varargin{:});

  if plotData
        Cdeepfun = @()modcoupler(N,M,m,0,1,0);
        Cdeep = Cdeepfun()+eye(N);
        rix = removeval(1:N,toRemove);
        [tPosRates,fPosRates] = calculate_thresh_acc(Rass,Cdeep(rix,rix,:),sims);
        plot(mean(fPosRates,2),mean(tPosRates,2),'-o','MarkerFaceColor',next_color);
  end
  
  if accessorize
      xlabel('false positive rate');
      ylabel('true positive rate');
      axis([-0.01 1 0 1.01]);
    if varargCheck('legend',varargin{:})  
      legend(varargAssign('legend',varargin{:}),'location','southeast');
    end
    if varargCheck('title',varargin{:})
      title('ROC');
    end
  end
  
end