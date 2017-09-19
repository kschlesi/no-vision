function h = bar_error(model_series,model_error,varargin)
% create barplot with error bars (model_series and model_error)

if ~varargCheck('ExtFigureCmd',varargin{:})
    figure;
end

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

legInput = varargAssign('legend',...
            cellfun(@num2str,num2cell(1:numbars),'UniformOutput',false),...
            varargin{:});
xtl = varargAssign('XTickLabel',1:numgroups,varargin{:});
gridstyle = varargAssign('GridLineStyle','-',varargin{:});
ebc = varargAssign('ErrBarColor','k',varargin{:});
barArgs = varargAssign('BarArgs',false,varargin{:});

h = bar(model_series);
set(h,'BarWidth',1);    % The bars will now touch each other
% set bar arguments
if iscell(barArgs) || barArgs~=false
  assert(size(barArgs,1)==numbars+1);
  for b=1:numbars
    set(h(b),barArgs{[1,b+1],:});
  end
end
set(gca,'YGrid','on')
set(gca,'GridLineStyle',gridstyle)
set(gca,'XTicklabel',xtl)
set(get(gca,'XLabel'),'String','Functional System')
set(get(gca,'YLabel'),'String','Mean System Recruitment')
lh = legend(legInput);
set(lh,'Location','BestOutside','Orientation','horizontal')
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, model_series(:,i), model_error(:,i), ebc, 'linestyle', 'none');
end

end