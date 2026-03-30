function DrawHC_WAIC_Set1()
%% ==============================================================================
%  Draw HC WAIC Comparison Only (Set1) - Baseline Marker Added
%  Data: Set1_data.csv

%
%  Style:
%   - Arial 9pt, FontWeight 'normal'
%   - LineWidth 0.6 (thin axis lines)
%   - Canvas: 5.5 cm x 6.0 cm
%   - Axes Height: fixed at 4.0 cm
%   - Margins: Left 1.15 cm, Right 0.35 cm
%   - Feature: M3c is marked with a thicker vertical line
% ==============================================================================
    clc; close all;
    
    % -----------------------------
    % 0) Global font settings: Arial 9pt, force normal weight
    % -----------------------------
    set(groot, 'defaultAxesFontName', 'Arial');
    set(groot, 'defaultTextFontName', 'Arial');
    set(groot, 'defaultAxesFontSize', 9);
    set(groot, 'defaultTextFontSize', 9);
    set(groot, 'defaultAxesFontWeight', 'normal');
    set(groot, 'defaultTextFontWeight', 'normal');

    % ---------- PATHS ----------
    waic_csv = 'F:\PIT\model\EEGmodelHC\HC\results\stan\plot\Set1_data.csv';
    save_dir = 'F:\PIT\model\EEGmodelHC\HC\results\stan\plot\';
    
    % ---------- SETTINGS ----------
    figW_cm  = 5.5;   % 55 mm
    figH_cm  = 6.0;   % 60 mm
    fontSz   = 9;
    fontName = 'Arial';
    
    hc_c_best  = [229, 157, 66]/255;
    hc_c_other = [246, 228, 198]/255;
    
    % Line settings
    axisLW    = 0.6;    % axis line width
    zeroLW    = 0.5;    % thin zero line across the panel
    baseMkLW  = 2.0;    % thick baseline marker for M3c
    barEdgeLW = 0.5;    % bar edge width

    % ---------- FIGURE ----------
    f = figure('Units','centimeters', ...
        'Position',[10 10 figW_cm figH_cm], ...
        'Color','w', ...
        'PaperPositionMode','auto', ...
        'InvertHardcopy','off');
    
    set(f, 'PaperUnits','centimeters');
    set(f, 'PaperPosition',[0 0 figW_cm figH_cm]);
    set(f, 'PaperSize',[figW_cm figH_cm]);
    set(f, 'Renderer','painters');

    % ============================================================
    % Size control: fixed height = 40 mm
    % ============================================================
    axH_cm = 4.0;     % fixed height = 40 mm
    left_cm   = 1.15; % left margin
    right_cm  = 0.35; % right margin
    bottom_cm = 1.35; % bottom margin
    axW_cm = figW_cm - left_cm - right_cm;

    ax1 = axes('Parent', f, 'Units','centimeters', ...
        'Position',[left_cm, bottom_cm, axW_cm, axH_cm]);
    hold(ax1,'on');

    % Axis style
    set(ax1, 'XColor','k', 'YColor','k', ...
        'LineWidth', axisLW, ...
        'FontSize', fontSz, ...
        'FontName', fontName, ...
        'FontWeight', 'normal', ...
        'TickDir', 'out', ...
        'Box', 'off');

    if exist(waic_csv,'file')
        T_waic = readtable(waic_csv);

        % 1) Model label mapping
        old_names = {'M9',  'M10', 'M12', 'M13', 'M14'};
        new_names = {'M3c', 'M4a', 'M4b', 'M4c', 'M4d'};
        waic_map  = containers.Map(old_names, new_names);
        
        T_waic.PaperLabel = cell(height(T_waic),1);
        for ii = 1:height(T_waic)
            mName = T_waic.model{ii};
            if isKey(waic_map, mName)
                T_waic.PaperLabel{ii} = waic_map(mName);
            else
                T_waic.PaperLabel{ii} = mName;
            end
        end

        % 2) Compute Delta
        base_model_raw = 'M9';
        base_idx = find(strcmp(T_waic.model, base_model_raw));
        if isempty(base_idx)
            warning('Ref model %s not found. Use max WAIC as baseline.', base_model_raw);
            [~, base_idx] = max(T_waic.WAIC);
        end
        base_val      = T_waic.WAIC(base_idx);
        T_waic.Delta  = T_waic.WAIC - base_val;
        T_waic.IsBest = (T_waic.WAIC == min(T_waic.WAIC));

        % 3) Reorder display
        disp_ord = {'M3c', 'M4a', 'M4b', 'M4c', 'M4d'};
        T_plot = table();
        for k = 1:numel(disp_ord)
            row = T_waic(strcmp(T_waic.PaperLabel, disp_ord{k}), :);
            if ~isempty(row)
                T_plot = [T_plot; row]; %#ok<AGROW>
            end
        end
        T_plot = flipud(T_plot);
        nb = height(T_plot);

        % 4) Draw bars
        for k = 1:nb
            val = T_plot.Delta(k);
            fc  = T_plot.IsBest(k)*hc_c_best + (~T_plot.IsBest(k))*hc_c_other;
            % LineWidth here controls the bar border
            barh(ax1, k, val, 0.65, 'FaceColor', fc, 'EdgeColor','k', 'LineWidth', barEdgeLW);
        end

        % 5) Set Y axis
        yL = [0.4 nb+0.6];
        set(ax1, 'YLim', yL, 'YTick', 1:nb, 'YTickLabel', T_plot.PaperLabel);

        % 6) Set X axis (reversed)
        min_v = min(T_plot.Delta);
        max_v = max(T_plot.Delta);
        limit_low  = min(min_v, 0);
        limit_high = max(max_v, 0);
        step = 50;
        grid_low  = floor(limit_low/step)*step;
        grid_high = ceil(limit_high/step)*step;
        if grid_high == grid_low
            grid_high = grid_low + step;
        end
        set(ax1, 'XLim', [grid_low grid_high], 'XDir','reverse');

        % 7) Draw zero line and baseline marker
        
        % A. Thin zero line across the full panel
        line(ax1, [0 0], yL, 'Color', 'k', 'LineWidth', zeroLW);
        
        % B. Thick short vertical line at M3c as the baseline marker
        idx_M3c = find(strcmp(T_plot.PaperLabel, 'M3c'));
        if ~isempty(idx_M3c)
            % [idx_M3c-0.35, idx_M3c+0.35] matches the bar height range
            hThick = line(ax1, [0 0], [idx_M3c-0.35, idx_M3c+0.35], ...
                'Color', 'k', 'LineWidth', baseMkLW);
            uistack(hThick, 'top'); % keep the thick line on top
        end

        % X-axis label
        xlabel(ax1, '{\Delta}log model evidence (WAIC)', ...
            'FontSize', fontSz, 'FontName', fontName, 'FontWeight', 'normal');
    else
        text(ax1, 0.5, 0.5, 'CSV Not Found', 'HorizontalAlignment','center', ...
            'FontName', fontName, 'FontSize', fontSz, 'FontWeight', 'normal');
    end

    % Final safeguard: enforce text properties
    set(findall(f, '-property', 'FontWeight'), 'FontWeight', 'normal');
    set(findall(f, '-property', 'FontName'), 'FontName', fontName);
    set(findall(f, '-property', 'FontSize'), 'FontSize', fontSz);

    % ---------- SAVE ----------
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    outBase = fullfile(save_dir, 'HC_Set1_WAIC_BaselineBold');
    
    fprintf('Saving figure to: %s\n', outBase);
    print(f, [outBase '_lineart1000.tif'], '-dtiff', '-r1000');
    print(f, [outBase '.svg'], '-dsvg', '-r1000');
    fprintf('Done.\n');
end