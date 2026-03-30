function DrawHC_WAIC_Set3()
%% ==============================================================================
%  Draw HC WAIC Comparison Only (Set3) - Baseline Marker Added
%  Data: Set3_data.csv
%  Mapping: M9->M3c, M12->M4b, M19->M5c, M21->M5d
%  Ref: M3c (raw M9) is set as baseline (0).
%
%  Style:
%   - Arial 9pt, FontWeight 'normal'
%   - LineWidth 0.6 (thinner lines for axes)
%   - Canvas: 5.5cm x 6.0cm
%   - Axes Height: Fixed 4.0cm
%   - Margins: Left 1.15cm, Right 0.35cm
%   - ★ Feature: M3c marked with a thicker vertical line
% ==============================================================================
    clc; close all;

    % -----------------------------
    % 0) Global font: Arial 9pt, normal (no bold)
    % -----------------------------
    set(groot, 'defaultAxesFontName', 'Arial');
    set(groot, 'defaultTextFontName', 'Arial');
    set(groot, 'defaultAxesFontSize', 9);
    set(groot, 'defaultTextFontSize', 9);
    set(groot, 'defaultAxesFontWeight', 'normal');
    set(groot, 'defaultTextFontWeight', 'normal');

    % ---------- PATHS (UPDATED) ----------
    waic_csv = 'F:\PIT\model\EEGmodelHC\HC\results\stan\plot\Set3_data.csv';
    save_dir = 'F:\PIT\ISPSplot\';

    % ---------- SETTINGS ----------
    figW_cm  = 5.5;   % 55 mm
    figH_cm  = 6.0;   % 60 mm
    fontSz   = 9;
    fontName = 'Arial';

    hc_c_best  = [229, 157, 66]/255;
    hc_c_other = [246, 228, 198]/255;

    axisLW    = 0.6;    % axis linewidth
    zeroLW    = 0.5;    % thin 0-line
    baseMkLW  = 2.0;    % thick baseline marker at M3c
    barEdgeLW = 0.5;    % bar edge

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
    % Axes size control: fixed height 40mm
    % ============================================================
    axH_cm = 4.0;
    left_cm   = 1.15;
    right_cm  = 0.35;
    bottom_cm = 1.35;
    axW_cm = figW_cm - left_cm - right_cm;

    ax1 = axes('Parent', f, 'Units','centimeters', ...
        'Position',[left_cm, bottom_cm, axW_cm, axH_cm]);
    hold(ax1,'on');

    set(ax1, 'XColor','k', 'YColor','k', ...
        'LineWidth', axisLW, ...
        'FontSize', fontSz, ...
        'FontName', fontName, ...
        'FontWeight', 'normal', ...
        'TickDir', 'out', ...
        'Box', 'off');

    if exist(waic_csv,'file')
        T_waic = readtable(waic_csv);

        % ---------- 1) Mapping (UPDATED) ----------
        old_names = {'M9',  'M12', 'M19', 'M21'};
        new_names = {'M3c', 'M4b', 'M5c', 'M5d'};
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

        % ---------- 2) Delta (baseline = M3c i.e., raw M9) ----------
        base_model_raw = 'M9';  % M3c
        base_idx = find(strcmp(T_waic.model, base_model_raw));
        if isempty(base_idx)
            warning('Ref model %s not found. Use max WAIC as baseline.', base_model_raw);
            [~, base_idx] = max(T_waic.WAIC);
        end
        base_val     = T_waic.WAIC(base_idx);
        T_waic.Delta = T_waic.WAIC - base_val;
        T_waic.IsBest = (T_waic.WAIC == min(T_waic.WAIC));

        % ---------- 3) Display order (UPDATED) ----------
        disp_ord = {'M3c','M4b','M5c','M5d'};
        T_plot = table();
        for k = 1:numel(disp_ord)
            row = T_waic(strcmp(T_waic.PaperLabel, disp_ord{k}), :);
            if ~isempty(row)
                T_plot = [T_plot; row]; %#ok<AGROW>
            end
        end


        T_plot = flipud(T_plot);
        nb = height(T_plot);

        % ---------- 4) Plot bars ----------
        for k = 1:nb
            val = T_plot.Delta(k);
            fc  = T_plot.IsBest(k)*hc_c_best + (~T_plot.IsBest(k))*hc_c_other;
            barh(ax1, k, val, 0.65, 'FaceColor', fc, 'EdgeColor','k', 'LineWidth', barEdgeLW);
        end

        % ---------- 5) Y axis labels ----------
        yL = [0.4 nb+0.6];
        set(ax1, 'YLim', yL, 'YTick', 1:nb, 'YTickLabel', T_plot.PaperLabel);

        % ---------- 6) X axis (reverse) ----------
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

        % ---------- 7) 0-line + baseline marker at M3c ----------
        line(ax1, [0 0], yL, 'Color', 'k', 'LineWidth', zeroLW);

        idx_M3c = find(strcmp(T_plot.PaperLabel, 'M3c'));
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

    % ---------- final font enforcement ----------
    set(findall(f, '-property', 'FontWeight'), 'FontWeight', 'normal');
    set(findall(f, '-property', 'FontName'), 'FontName', fontName);
    set(findall(f, '-property', 'FontSize'), 'FontSize', fontSz);

    % ---------- SAVE (UPDATED DIR) ----------
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    outBase = fullfile(save_dir, 'HC_Set3_WAIC_BaselineBold');

    fprintf('Saving figure to: %s\n', outBase);
    print(f, [outBase '_lineart1000.tif'], '-dtiff', '-r1000');
    print(f, [outBase '.svg'], '-dsvg', '-r1000');
    fprintf('Done.\n');
end
