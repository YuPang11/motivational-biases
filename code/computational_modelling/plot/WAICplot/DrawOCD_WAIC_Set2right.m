function DrawOCD_WAIC_Set2right()
%% ==============================================================================
%  Draw OCD WAIC Comparison Only (Set2) - Baseline Marker Added
%  Data: Set2_dataOCD.csv
%  Mapping: M5->M3c (Baseline), M10->M4a, M15->M5a, M16->M5b
%  Ref: M5 (displayed as M3c) is set as baseline (0).
%
%  Style:
%   - OCD Color Theme (Blue)
%   - Arial 9pt, FontWeight 'normal'
%   - Canvas: 5.5 cm x 6.0 cm
%   - Axes Height: Fixed 4.0 cm
%   - Feature: Explicitly draw M3c baseline marker (thick line)
%
%  Current requirements:
%   - Y axis on the RIGHT
%   - X axis shown as 0 (right) -> negative (left)
%   - Reduce left whitespace and increase right whitespace
% ==============================================================================

    clc; close all;

    % -----------------------------
    % 0) Global font settings: Arial 9pt, normal weight
    % -----------------------------
    set(groot, 'defaultAxesFontName', 'Arial');
    set(groot, 'defaultTextFontName', 'Arial');
    set(groot, 'defaultAxesFontSize', 9);
    set(groot, 'defaultTextFontSize', 9);
    set(groot, 'defaultAxesFontWeight', 'normal');
    set(groot, 'defaultTextFontWeight', 'normal');

    % ---------- PATHS ----------
    waic_csv = 'F:\PIT\model\EEGmodelHC\OCD\results\stan\plot\Set2_dataOCD.csv';
    save_dir = 'F:\PIT\ISPSplot\';

    % ---------- SETTINGS ----------
    figW_cm  = 5.5;   % 55 mm
    figH_cm  = 6.0;   % 60 mm
    fontSz   = 9;
    fontName = 'Arial';

    % OCD blue theme
    ocd_c_best  = [0.25 0.60 0.90];
    ocd_c_other = [0.70 0.85 1.00];

    % Line settings
    axisLW    = 0.6;
    zeroLW    = 0.5;
    baseMkLW  = 2.0;   % thick marker for baseline M3c
    barEdgeLW = 0.5;

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
    % Fixed axes height = 40 mm
    % ============================================================
    axH_cm = 4.0;

    % Shift plot left and leave more space on the right
    % Original: left_cm = 1.15, right_cm = 0.35
    % Current:  left_cm = 0.45, right_cm = 1.05
    left_cm   = 0.45;
    right_cm  = 1.05;

    bottom_cm = 1.35;

    % Axes width is recalculated automatically
    axW_cm = figW_cm - left_cm - right_cm;

    ax1 = axes('Parent', f, 'Units','centimeters', ...
        'Position',[left_cm, bottom_cm, axW_cm, axH_cm]);
    hold(ax1,'on');

    % Put Y axis on the right
    set(ax1, 'XColor','k', 'YColor','k', ...
        'LineWidth', axisLW, ...
        'FontSize', fontSz, ...
        'FontName', fontName, ...
        'FontWeight', 'normal', ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'YAxisLocation', 'right');

    if exist(waic_csv,'file')
        T_waic = readtable(waic_csv);

        % ---------- Ensure required columns exist ----------
        if ~ismember('model', T_waic.Properties.VariableNames)
            error('CSV must contain a column named "model".');
        end
        if ~ismember('WAIC', T_waic.Properties.VariableNames)
            error('CSV must contain a column named "WAIC".');
        end

        % Ensure model is cellstr
        if isstring(T_waic.model)
            model_raw = cellstr(T_waic.model);
        elseif iscell(T_waic.model)
            model_raw = T_waic.model;
        else
            model_raw = cellstr(string(T_waic.model));
        end

        % ---------- 1) Mapping (Set2 OCD) ----------
        old_names = {'M5',  'M10', 'M15', 'M16'};
        new_names = {'M3c', 'M4a', 'M5a', 'M5b'};
        waic_map  = containers.Map(old_names, new_names);

        T_waic.PaperLabel = cell(height(T_waic),1);
        for ii = 1:height(T_waic)
            mName = model_raw{ii};
            if isKey(waic_map, mName)
                T_waic.PaperLabel{ii} = waic_map(mName);
            else
                T_waic.PaperLabel{ii} = mName;
            end
        end

        % ---------- 2) Delta (baseline = M5 -> M3c) ----------
        base_model_raw = 'M5';
        base_idx = find(strcmp(model_raw, base_model_raw), 1);

        if isempty(base_idx)
            warning('Ref model %s not found. Using max WAIC as baseline.', base_model_raw);
            [~, base_idx] = max(T_waic.WAIC);
        end

        base_val      = T_waic.WAIC(base_idx);
        T_waic.Delta  = T_waic.WAIC - base_val;
        T_waic.IsBest = (T_waic.WAIC == min(T_waic.WAIC));

        % ---------- 3) Display order ----------
        disp_ord = {'M3c','M4a','M5a','M5b'};
        T_plot = table();
        for k = 1:numel(disp_ord)
            row = T_waic(strcmp(T_waic.PaperLabel, disp_ord{k}), :);
            if ~isempty(row)
                T_plot = [T_plot; row]; %#ok<AGROW>
            end
        end

        % Keep original style: first item displayed at the bottom
        T_plot = flipud(T_plot);
        nb = height(T_plot);

        % ---------- 4) Draw bars ----------
        for k = 1:nb
            val = T_plot.Delta(k);
            fc  = T_plot.IsBest(k)*ocd_c_best + (~T_plot.IsBest(k))*ocd_c_other;
            barh(ax1, k, val, 0.65, 'FaceColor', fc, 'EdgeColor','k', 'LineWidth', barEdgeLW);
        end

        % ---------- 5) Y axis ----------
        yL = [0.4 nb+0.6];
        set(ax1, 'YLim', yL, 'YTick', 1:nb, 'YTickLabel', T_plot.PaperLabel);

        % ---------- 6) X axis: 0 on right, negative on left ----------
        min_v = min(T_plot.Delta);
        limit_low = min(min_v, 0);

        step = 50;
        grid_low = floor(limit_low/step)*step;
        if grid_low == 0
            grid_low = -step;
        end

        set(ax1, 'XLim', [grid_low 0], 'XDir','normal');

        % ---------- 7) Zero line and baseline marker ----------
        line(ax1, [0 0], yL, 'Color', 'k', 'LineWidth', zeroLW);

        idx_M3c = find(strcmp(T_plot.PaperLabel, 'M3c'), 1);
        if ~isempty(idx_M3c)
            hThick = line(ax1, [0 0], [idx_M3c-0.35, idx_M3c+0.35], ...
                'Color', 'k', 'LineWidth', baseMkLW);
            uistack(hThick, 'top');
        end

        xlabel(ax1, '{\Delta}log model evidence (WAIC)', ...
            'FontSize', fontSz, 'FontName', fontName, 'FontWeight', 'normal');

    else
        text(ax1, 0.5, 0.5, 'CSV Not Found', 'HorizontalAlignment','center', ...
            'FontName', fontName, 'FontSize', fontSz, 'FontWeight', 'normal');
    end

    % ---------- Final enforcement ----------
    set(findall(f, '-property', 'FontWeight'), 'FontWeight', 'normal');
    set(findall(f, '-property', 'FontName'), 'FontName', fontName);
    set(findall(f, '-property', 'FontSize'), 'FontSize', fontSz);

    % ---------- SAVE ----------
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    outBase = fullfile(save_dir, 'OCD_Set2_WAIC_BaselineBold');

    fprintf('Saving figure to: %s\n', outBase);
    print(f, [outBase '_lineart1000OCD.tif'], '-dtiff', '-r1000');
    print(f, [outBase '.svg'], '-dsvg', '-r1000');
    fprintf('Done.\n');
end