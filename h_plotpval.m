function h_plotpval(pm)

    drawnow
    ylim = get(gca, 'YLim');
    if ylim(2)<0.5
        set(gca, 'YLim', [ylim(1) 0.5]);
        ylim = get(gca, 'YLim')
    end
    dy = ylim(2)-ylim(1)
    ypos = ylim(2)-0.1*dy;
    if pm<0.001
        pval = sprintf('%1.1e', pm)
    else
        pval = sprintf('%1.3f', pm)
    end
    plot([1 3], ypos*[1 1], 'k')
    text(2,ypos+0.02*dy,['pval = ' pval], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline')
    %plot([0.8 1 1.2], ypos*[1 1 1]-[1 0 1]*0.05, 'k')
    plot([-1 1]*0.25+1, ypos*[1 1]-0.025*dy, 'k')
    plot([-1 1]*0.25+3, ypos*[1 1]-0.025*dy, 'k')
    plot([1 1], ypos - [1 0]*0.025*dy, 'k')
    plot([3 3], ypos - [1 0]*0.025*dy, 'k')
    %plot([-1 -1]*0.2+1, ypos*[1-0.05 1-0.025], 'k')