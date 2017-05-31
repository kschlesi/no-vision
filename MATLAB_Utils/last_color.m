function color = last_color()
% returns last color used in current axes
        colors = get(gca,'ColorOrder');
        index  = get(gca,'ColorOrderIndex');
        if index<=1
            %disp([num2str(index) ', checking 1']);
            index = 1;
        else
            %disp([num2str(index) ', checking ' num2str(index-1)])
            index = index-1;
        end
        color = colors(index,:);
end