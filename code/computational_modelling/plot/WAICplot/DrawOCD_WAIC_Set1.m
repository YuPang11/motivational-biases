function DrawOCD_WAIC_Set1()
%% ==============================================================================
%  Draw OCD WAIC Comparison (Set1) - Final Version
%  Data: Set1_dataOCD.csv
%  Mapping: M5->M3c (Baseline), M10->M4a, M12->M4b, M13->M4c, M14->M4d
%  Ref: M5 (displayed as M3c) is set as baseline (0).
%
%  Style:
%   - OCD color theme (blue)
%   - Arial 9pt, FontWeight 'normal'
%   - Canvas: 5.5 cm x 6.0 cm
%   - Axes height: fixed at 4.0 cm
%   - Margins: Left 1.15 cm, Right 0.35 cm
%   - Feature: explicitly draw the M3c baseline marker (thick line)
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
    % Load OCD data
    waic_csv = 'F:\PIT\model\EEGmodelHC\HC\results\stan\plot\Set1_dataOCD.csv';
    save_dir = 'F:\PIT\model\EEGmodelHC\HC\results\stan\plot\';
    
    % ---------- SETTINGS ----------
    figW_cm  = 5.5;   % 55 mm
    figH_cm  = 6.0;   % 60 mm
    fontSz   = 9;
    fontName = 'Arial';
    
    % OCD color theme (blue)
    ocd_c_best  = [0.25 0.60 0.90];      % dark blue (best model)
    ocd_c_other = [0.70 0.85 1.0];       % light blue (other models)
    
    % Line settings
    axisLW    = 0.6;    % axis line width
    zeroLW    = 0.5;    % thin zero line across the panel
    baseMkLW  = 2.0;    % M3c baseline marker width (bold)
    barEdgeLW = 0.5;    % bar border width

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
        
        % 1) Label mapping (key change: M5 -> M3c)
        % Assume the CSV uses model labels such as M5, M10, M12...
        old_names = {'M5',  'M10', 'M12', 'M13', 'M14'};
        new_names = {'M3c', 'M4a', 'M4b', 'M4c', 'M4d'};
        waic_map  = containers.Map(old_names, new_names);
        
        T_waic.PaperLabel = cell(height(T_waic),1);
        for ii = 1:height(T_waic)
            mName = T_waic.model{ii};
            if isKey(waic_map, mName)
                T_waic.PaperLabel{ii} = waic_map(mName);
            else
                T_waic.PaperLabel{ii} = mName; % fallback
            end
        end

        % 2) Compute Delta (baseline model = M5/M3c)
        base_model_raw = 'M5'; 
        base_idx = find(strcmp(T_waic.model, base_model_raw));
        
        if isempty(base_idx)
            % If M5 is not found in the CSV, try M3c or use max WAIC as fallback
            alt_base = find(strcmp(T_waic.PaperLabel, 'M3c'));
            if ~isempty(alt_base)
                base_idx = alt_base;
            else
                warning('Ref model %s not found. Using max WAIC as baseline.', base_model_raw);
                [~, base_idx] = max(T_waic.WAIC);
            end
        end
        
        base_val      = T_waic.WAIC(base_idx);
        T_waic.Delta  = T_waic.WAIC - base_val;
        T_waic.IsBest = (T_waic.WAIC == min(T_waic.WAIC));

        % 3) Reorder rows (make sure M3c is included)
        disp_ord = {'M3c', 'M4a', 'M4b', 'M4c', 'M4d'};
        T_plot = table();
        for k = 1:numel(disp_ord)
            row = T_waic(strcmp(T_waic.PaperLabel, disp_ord{k}), :);
            if ~isempty(row)
                T_plot = [T_plot; row]; %#ok<AGROW>
            end
        end
        T_plot = flipud(T_plot); % reverse order so the first item is shown on top
        nb = height(T_plot);

        % 4) Draw horizontal bars
        for k = 1:nb
            val = T_plot.Delta(k);
            % Use OCD blue color theme
            fc  = T_plot.IsBest(k)*ocd_c_best + (~T_plot.IsBest(k))*ocd_c_other;
            % Draw bar
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
        
        % A. Thin zero line across the full panel (background reference)
        line(ax1, [0 0], yL, 'Color', 'k', 'LineWidth', zeroLW);
        
        % B. Explicitly draw the M3c baseline marker (thick black line)
        idx_M3c = find(strcmp(T_plot.PaperLabel, 'M3c'));
        
        if ~isempty(idx_M3c)
            % Draw a short thick vertical line at the Y position of M3c
            % The line spans the height of the bar (0.7 units)
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
    % Output filename for OCD
    outBase = fullfile(save_dir, 'OCD_Set1_WAIC_WithM3c');
    
    fprintf('Saving figure to: %s\n', outBase);
    print(f, [outBase '_lineart1000OCD.tif'], '-dtiff', '-r1000');
    print(f, [outBase '.svg'], '-dsvg', '-r1000');
    fprintf('Done.\n');
end