function color = next_color()
% returns color marked for use in current axes
        colors = get(gca,'ColorOrder');
        index  = get(gca,'ColorOrderIndex');
        color = colors(index,:);
end