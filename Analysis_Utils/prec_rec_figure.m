function [] = prec_rec_figure(Rass,N,M,m,toRemove,sims,varargin)
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
      [recall,~,~,~,precision] = calculate_thresh_acc(Rass,Cdeep(rix,rix,:),sims);
      if varargCheck(':plot',varargin{:})
        plot(mean(recall,2),mean(precision,2),':o','Color',last_color,...
                                              'MarkerFaceColor',last_color);  
      else
        plot(mean(recall,2),mean(precision,2),'-o','MarkerFaceColor',next_color);
      end
  end
  
  if accessorize
      xlabel('recall (# true pos / all true)');
        ylabel('precision (# true pos / all pos)');
        axis([0 1.05 0 1.05]);
    if varargCheck('legend',varargin{:})  
      legend(varargAssign('legend',varargin{:}),'location','southwest');
    end
    if varargCheck('title',varargin{:})
      title('precision-recall curves');
    end
  end
  
end