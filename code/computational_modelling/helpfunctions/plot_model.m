function plot_model(ax, x, mu, se, fillCol, lwModel, hasBL)
    
    if hasBL
        [hl, hp] = boundedline(ax, x, mu, se, 'alpha', 'cmap', [0 0 0]);
        
        try, set(hp, 'FaceColor', fillCol, 'FaceAlpha', 0.35, 'EdgeColor', 'none'); end
        if isgraphics(hl)
            set(hl, 'Color', 'k', 'LineWidth', lwModel);
        elseif isstruct(hl) && isfield(hl, 'mainLine')
            set(hl.mainLine, 'Color', 'k', 'LineWidth', lwModel);
            if isfield(hl, 'patch')
                set(hl.patch, 'FaceColor', fillCol, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            end
        end
    else
        
        px = [x, fliplr(x)];
        py = [mu + se, fliplr(mu - se)];
        patch('Parent', ax, 'XData', px, 'YData', py, ...
              'FaceColor', fillCol, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        plot(ax, x, mu, 'k-', 'LineWidth', lwModel);
    end
end
