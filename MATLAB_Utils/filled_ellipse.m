function h = filled_ellipse(ra,rb,varargin)
% wrapper for ELLIPSE.M, written by KJ Schlesinger, June 2017

patchColor = varargAssign('ellipseColor',[1 1 1],varargin{:});
patchAlpha = varargAssign('ellipseAlpha',1,varargin{:});

h = ellipse(ra,rb,varargin{:});
x = get(h,'Xdata');
y = get(h,'Ydata');
hold on;
patch(x,y,patchColor,'LineStyle','None','FaceAlpha',patchAlpha);

end