function color = last_color(n)

if nargin==0
    n = 1;
end
% returns last color used in current axes
        colors = get(gca,'ColorOrder');
        index  = get(gca,'ColorOrderIndex');
        if index<=n
            %disp([num2str(index) ', checking 1']);
            index = 1;
        else
            %disp([num2str(index) ', checking ' num2str(index-1)])
            index = index-n;
        end
        color = colors(index,:);
end