%% =========================================================
% Partial Spearman correlations of ISPSDeg / Radius_z / Proj_impulse_z
% with each scale score (X)
% Controls: BDI, BAI (residual method)
% Data: F:\PIT\model\EEGmodelHC\OCD\OCD_M5b_Q.xlsx
% Output: F:\PIT\ISPSplot\relation  (CSV + one 80x70 mm SVG per combination)
%
% Plot style:
%  - Canvas: 80×70 mm
%  - Main plot: Arial 9 pt, normal
%  - Legend: top, Arial 8 pt, normal
%  - Black data points (smaller size)
%  - Blue regression line + blue 95% CI shading
%  - Fixed physical axis lengths: X = 4.0 cm, Y = 4.1 cm
%
% Revision notes:
%  - Keep Excel column name as ThetaDeg (do not rename data column)
%  - Display y-axis label as φ_reg (reg as subscript)
%  - Replace ThetaDeg with phi_reg in output filenames as well
%% =========================================================
clc; clear; close all;

%% ---- 0) Paths ------------------------------------------------------
xlsx_path = 'F:\PIT\model\EEGmodelHC\OCD\OCD_M5b_Q.xlsx';
out_dir   = 'F:\PIT\ISPSplot\relation';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ---- 1) Global font fallback: Arial + normal -----------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultLegendFontName', 'Arial');
set(groot, 'defaultLegendFontWeight', 'normal');

%% ---- 2) Read table -------------------------------------------------
if ~isfile(xlsx_path)
    error('Data file not found: %s', xlsx_path);
end
T = readtable(xlsx_path, 'VariableNamingRule','preserve');

%% ---- 3) Variable settings -----------------------------------------
x_vars = {'YBOCSO','YBOCSC','YBOCS','OCR', ...
          'BIS','Responsibility','Perfectionism','Control'};

% Actual Excel column name remains ThetaDeg
y_vars = {'ThetaDeg','Radius_z','Proj_impulse_z'};

covars = {'BDI','BAI'};

need_cols = [x_vars, y_vars, covars];
missing_cols = setdiff(need_cols, T.Properties.VariableNames);
if ~isempty(missing_cols)
    error('The following columns are missing in Excel:\n%s', strjoin(missing_cols, ', '));
end

%% ---- 4) Preallocate output table -----------------------------------
nTotal = numel(x_vars) * numel(y_vars);
pred_col  = cell(nTotal,1);
out_col   = cell(nTotal,1);
n_col     = nan(nTotal,1);
r_col     = nan(nTotal,1);
p_col     = nan(nTotal,1);
plot_col  = cell(nTotal,1);
rowk = 0;

%% ---- 5) Unified plotting parameters --------------------------------
% Canvas size: 80×70 mm
figW_cm = 8.0;   % 80 mm
figH_cm = 7.0;   % 70 mm

% Fixed physical axis lengths
axW_cm = 4.0;    % X-axis length 4 cm
axH_cm = 4.1;    % Y-axis length 4.1 cm

% Axis position
axLeft_cm   = 1.6;
axBottom_cm = 1.3;

% Legend above the main plot (same width as main plot)
lgdGap_cm   = 0.30;
lgdH_cm     = 0.75;
lgdPos_cm   = [axLeft_cm, axBottom_cm + axH_cm + lgdGap_cm, axW_cm, lgdH_cm];

blueCol = [0.20 0.55 0.95];
ciAlpha = 0.18;

ptSize  = 8;     % Smaller black points
lwLine  = 1.6;
lwAxis  = 0.8;

fontMain = 9;    % Main plot font size
fontLgd  = 8;    % Legend font size

%% ---- 6) Loop: partial Spearman for all Y × X -----------------------
for yi = 1:numel(y_vars)
    y_name = y_vars{yi};
    fprintf('>>> Analyzing Y = %s\n', y_name);

    for xi = 1:numel(x_vars)
        x_name = x_vars{xi};
        rowk = rowk + 1;

        pred_col{rowk} = x_name;
        out_col{rowk}  = y_name;

        % ---- Extract columns and convert to double
        %      (compatible with numeric/string/cellstr/categorical) ----
        x  = str2double(string(T.(x_name)));
        y  = str2double(string(T.(y_name)));
        c1 = str2double(string(T.(covars{1})));
        c2 = str2double(string(T.(covars{2})));

        ok = isfinite(x) & isfinite(y) & isfinite(c1) & isfinite(c2);
        x = x(ok); y = y(ok); c1 = c1(ok); c2 = c2(ok);

        n = numel(x);
        n_col(rowk) = n;

        if n < 6
            warning('Sample size too small, skipped: X=%s, Y=%s (n=%d)', x_name, y_name, n);
            r_col(rowk) = NaN;
            p_col(rowk) = NaN;
            plot_col{rowk} = '';
            continue;
        end

        % ---- Residual method: control BDI and BAI ----
        Xcov = [ones(n,1), c1, c2];
        bx = Xcov \ x;
        by = Xcov \ y;
        x_resid = x - Xcov * bx;
        y_resid = y - Xcov * by;

        % ---- Spearman correlation between residuals ----
        [r_s, p_s] = corr(x_resid, y_resid, 'Type','Spearman', 'Rows','complete');
        r_col(rowk) = r_s;
        p_col(rowk) = p_s;

        % ---- Linear fit + 95% CI (for visualization only) ----
        mdl = fitlm(x_resid, y_resid);
        xx = linspace(min(x_resid), max(x_resid), 200)';
        [yhat, yCI] = predict(mdl, xx, 'Alpha', 0.05, 'Prediction','curve');

        % ---- Create figure (80×70 mm) ----
        f = figure('Color','w', 'Units','centimeters', ...
            'Position',[5 5 figW_cm figH_cm], ...
            'PaperUnits','centimeters', ...
            'PaperPosition',[0 0 figW_cm figH_cm], ...
            'PaperSize',[figW_cm figH_cm], ...
            'InvertHardcopy','off');

        ax = axes('Parent', f, 'Units','centimeters', ...
            'Position', [axLeft_cm axBottom_cm axW_cm axH_cm]);
        hold(ax, 'on');

        % 95% CI shading
        fill(ax, [xx; flipud(xx)], [yCI(:,1); flipud(yCI(:,2))], blueCol, ...
            'EdgeColor','none', 'FaceAlpha', ciAlpha);

        % Black points
        scatter(ax, x_resid, y_resid, ptSize, 'o', ...
            'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth',0.5);

        % Blue line
        plot(ax, xx, yhat, '-', 'Color', blueCol, 'LineWidth', lwLine);

        % Axis labels
        xlabel(ax, 'Y-BOCS obsessions', 'FontName','Arial', 'FontSize',fontMain, ...
            'FontWeight','normal', 'Interpreter','none');

        % y-axis label: display φ_reg if y is ThetaDeg
        if strcmp(y_name, 'ThetaDeg')
            ylabel(ax, '\phi_{reg}', 'FontName','Arial', 'FontSize',fontMain, ...
                'FontWeight','normal', 'Interpreter','tex');
        else
            ylabel(ax, y_name, 'FontName','Arial', 'FontSize',fontMain, ...
                'FontWeight','normal', 'Interpreter','none');
        end

        set(ax, 'Box','off', 'TickDir','out', 'LineWidth', lwAxis, ...
            'FontName','Arial', 'FontSize',fontMain, 'FontWeight','normal');

        % Padding (data range)
        xpad = 0.05 * range(x_resid); if xpad==0, xpad = 0.01; end
        ypad = 0.05 * range(y_resid); if ypad==0, ypad = 0.01; end
        xlim(ax, [min(x_resid)-xpad, max(x_resid)+xpad]);
        ylim(ax, [min(y_resid)-ypad, max(y_resid)+ypad]);

        % ---- Top legend (Arial 8 pt) ----
        hCI_leg   = patch(NaN, NaN, blueCol, 'FaceAlpha', ciAlpha, 'EdgeColor','none');
        hData_leg = plot(NaN, NaN, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4, 'LineWidth',0.6);
        hLine_leg = plot(NaN, NaN, '-', 'Color', blueCol, 'LineWidth', lwLine);

        rp_txt = sprintf('r=%.3f   p=%.3f', r_s, p_s);

        lgd = legend(ax, [hCI_leg, hData_leg, hLine_leg], {'95% CI','Data', rp_txt}, ...
            'Location','northoutside', 'Orientation','horizontal', 'Box','on');
        lgd.FontName   = 'Arial';
        lgd.FontSize   = fontLgd;
        lgd.FontWeight = 'normal';
        lgd.Interpreter = 'none';
        lgd.ItemTokenSize = [10 6];

        lgd.Units = 'centimeters';
        lgd.Position = lgdPos_cm;

        % ---- Save SVG (vector) ----
        % Output filename: use phi_reg if y is ThetaDeg
        if strcmp(y_name, 'ThetaDeg')
            y_tag = 'phi_reg';
        else
            y_tag = y_name;
        end

        base_name = sprintf('%s_vs_%s_partialSpearman_BDI_BAI_80x70mm', y_tag, x_name);
        svg_file  = fullfile(out_dir, [base_name '.svg']);

        try
            exportgraphics(f, svg_file, 'ContentType','vector');
        catch
            print(f, svg_file, '-dsvg');
        end

        close(f);
        plot_col{rowk} = svg_file;
    end
end

%% ---- 7) Save result table (CSV) ------------------------------------
results = table(pred_col, out_col, n_col, r_col, p_col, plot_col, ...
    'VariableNames', {'predictor','outcome','n','spearman_r','spearman_p','plot_file'});

results_export = results(:, {'predictor','outcome','n','spearman_r','spearman_p'});

% CSV filename can also use phi_reg for consistency
out_csv = fullfile(out_dir, 'phi_reg_Radius_ProjImpulse_partial_Spearman_BDI_BAI.csv');
writetable(results_export, out_csv);

fprintf('\n✅ Done!\n- Result table: %s\n- Figure directory: %s\n', out_csv, out_dir);
disp(results_export);




%% =========================================================
% Partial Spearman correlations of ISPSDeg / Radius_z / Proj_impulse_z
% with each scale score (X)
% Controls: BDI, BAI, BIS (residual method)
% Data: F:\PIT\model\EEGmodelHC\OCD\OCD_M5b_Q.xlsx
% Output: F:\PIT\ISPSplot\relation  (CSV + one 80x70 mm SVG per combination)
%
% Plot style:
%  - Canvas: 80×70 mm
%  - Main plot: Arial 9 pt, normal
%  - Legend: top, Arial 8 pt, normal
%  - Black data points (smaller size)
%  - Blue regression line + blue 95% CI shading
%  - Fixed physical axis lengths: X = 4.0 cm, Y = 4.1 cm
%
% Revision notes:
%  - Keep Excel column name as ThetaDeg (do not rename data column)
%  - Display y-axis label as φ_reg (reg as subscript)
%  - Replace ThetaDeg with phi_reg in output filenames as well
%  - Add BIS as a covariate (now covars = {BDI, BAI, BIS})
%% =========================================================
clc; clear; close all;

%% ---- 0) Paths ------------------------------------------------------
xlsx_path = 'F:\PIT\model\EEGmodelHC\OCD\OCD_M5b_Q.xlsx';
out_dir   = 'F:\PIT\ISPSplot\relation';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ---- 1) Global font fallback: Arial + normal -----------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultLegendFontName', 'Arial');
set(groot, 'defaultLegendFontWeight', 'normal');

%% ---- 2) Read table -------------------------------------------------
if ~isfile(xlsx_path)
    error('Data file not found: %s', xlsx_path);
end
T = readtable(xlsx_path, 'VariableNamingRule','preserve');

%% ---- 3) Variable settings -----------------------------------------
x_vars = {'YBOCSO','YBOCSC','YBOCS','OCR', ...
          'BIS','Responsibility','Perfectionism','Control'};

% Actual Excel column name remains ThetaDeg
y_vars = {'ThetaDeg','Radius_z','Proj_impulse_z'};

% Covariates: BDI, BAI, BIS
covars = {'BDI','BAI','BIS'};

need_cols = [x_vars, y_vars, covars];
missing_cols = setdiff(need_cols, T.Properties.VariableNames);
if ~isempty(missing_cols)
    error('The following columns are missing in Excel:\n%s', strjoin(missing_cols, ', '));
end

%% ---- 4) Preallocate output table -----------------------------------
nTotal = numel(x_vars) * numel(y_vars);
pred_col  = cell(nTotal,1);
out_col   = cell(nTotal,1);
n_col     = nan(nTotal,1);
r_col     = nan(nTotal,1);
p_col     = nan(nTotal,1);
plot_col  = cell(nTotal,1);
rowk = 0;

%% ---- 5) Unified plotting parameters --------------------------------
figW_cm = 8.0;   % 80 mm
figH_cm = 7.0;   % 70 mm

axW_cm = 4.0;    % X-axis length 4 cm
axH_cm = 4.1;    % Y-axis length 4.1 cm

axLeft_cm   = 1.6;
axBottom_cm = 1.3;

lgdGap_cm   = 0.30;
lgdH_cm     = 0.75;
lgdPos_cm   = [axLeft_cm, axBottom_cm + axH_cm + lgdGap_cm, axW_cm, lgdH_cm];

blueCol = [0.20 0.55 0.95];
ciAlpha = 0.18;

ptSize  = 8;
lwLine  = 1.6;
lwAxis  = 0.8;

fontMain = 9;
fontLgd  = 8;

%% ---- 6) Loop: partial Spearman for all Y × X -----------------------
for yi = 1:numel(y_vars)
    y_name = y_vars{yi};
    fprintf('>>> Analyzing Y = %s\n', y_name);

    for xi = 1:numel(x_vars)
        x_name = x_vars{xi};
        rowk = rowk + 1;

        pred_col{rowk} = x_name;
        out_col{rowk}  = y_name;

        % ---- Extract X/Y and convert to double ----
        x = str2double(string(T.(x_name)));
        y = str2double(string(T.(y_name)));

        % ---- Extract covariate matrix C (one covariate per column) ----
        C = nan(height(T), numel(covars));
        for ci = 1:numel(covars)
            C(:,ci) = str2double(string(T.(covars{ci})));
        end

        % ---- Remove missing values: x, y, and all covariates must be valid ----
        ok = isfinite(x) & isfinite(y) & all(isfinite(C), 2);
        x = x(ok); y = y(ok); C = C(ok,:);

        n = numel(x);
        n_col(rowk) = n;

        if n < 6
            warning('Sample size too small, skipped: X=%s, Y=%s (n=%d)', x_name, y_name, n);
            r_col(rowk) = NaN;
            p_col(rowk) = NaN;
            plot_col{rowk} = '';
            continue;
        end

        % ---- Residual method: control BDI, BAI, and BIS ----
        Xcov = [ones(n,1), C];
        bx = Xcov \ x;
        by = Xcov \ y;
        x_resid = x - Xcov * bx;
        y_resid = y - Xcov * by;

        % ---- Spearman correlation between residuals ----
        [r_s, p_s] = corr(x_resid, y_resid, 'Type','Spearman', 'Rows','complete');
        r_col(rowk) = r_s;
        p_col(rowk) = p_s;

        % ---- Linear fit + 95% CI (for visualization only) ----
        mdl = fitlm(x_resid, y_resid);
        xx = linspace(min(x_resid), max(x_resid), 200)';
        [yhat, yCI] = predict(mdl, xx, 'Alpha', 0.05, 'Prediction','curve');

        % ---- Create figure (80×70 mm) ----
        f = figure('Color','w', 'Units','centimeters', ...
            'Position',[5 5 figW_cm figH_cm], ...
            'PaperUnits','centimeters', ...
            'PaperPosition',[0 0 figW_cm figH_cm], ...
            'PaperSize',[figW_cm figH_cm], ...
            'InvertHardcopy','off');

        ax = axes('Parent', f, 'Units','centimeters', ...
            'Position', [axLeft_cm axBottom_cm axW_cm axH_cm]);
        hold(ax, 'on');

        % 95% CI shading
        fill(ax, [xx; flipud(xx)], [yCI(:,1); flipud(yCI(:,2))], blueCol, ...
            'EdgeColor','none', 'FaceAlpha', ciAlpha);

        % Black points
        scatter(ax, x_resid, y_resid, ptSize, 'o', ...
            'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth',0.5);

        % Blue line
        plot(ax, xx, yhat, '-', 'Color', blueCol, 'LineWidth', lwLine);

        % Axis labels
        xlabel(ax, 'Y-BOCS obsessions', 'FontName','Arial', 'FontSize',fontMain, ...
            'FontWeight','normal', 'Interpreter','none');

        % y-axis label: ThetaDeg -> φ_reg
        if strcmp(y_name, 'ThetaDeg')
            ylabel(ax, '\phi_{reg}', 'FontName','Arial', 'FontSize',fontMain, ...
                'FontWeight','normal', 'Interpreter','tex');
        else
            ylabel(ax, y_name, 'FontName','Arial', 'FontSize',fontMain, ...
                'FontWeight','normal', 'Interpreter','none');
        end

        set(ax, 'Box','off', 'TickDir','out', 'LineWidth', lwAxis, ...
            'FontName','Arial', 'FontSize',fontMain, 'FontWeight','normal');

        % Padding (data range)
        xpad = 0.05 * range(x_resid); if xpad==0, xpad = 0.01; end
        ypad = 0.05 * range(y_resid); if ypad==0, ypad = 0.01; end
        xlim(ax, [min(x_resid)-xpad, max(x_resid)+xpad]);
        ylim(ax, [min(y_resid)-ypad, max(y_resid)+ypad]);

        % ---- Top legend (Arial 8 pt) ----
        hCI_leg   = patch(NaN, NaN, blueCol, 'FaceAlpha', ciAlpha, 'EdgeColor','none');
        hData_leg = plot(NaN, NaN, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4, 'LineWidth',0.6);
        hLine_leg = plot(NaN, NaN, '-', 'Color', blueCol, 'LineWidth', lwLine);

        rp_txt = sprintf('r=%.3f   p=%.3f', r_s, p_s);

        lgd = legend(ax, [hCI_leg, hData_leg, hLine_leg], {'95% CI','Data', rp_txt}, ...
            'Location','northoutside', 'Orientation','horizontal', 'Box','on');
        lgd.FontName   = 'Arial';
        lgd.FontSize   = fontLgd;
        lgd.FontWeight = 'normal';
        lgd.Interpreter = 'none';
        lgd.ItemTokenSize = [10 6];

        lgd.Units = 'centimeters';
        lgd.Position = lgdPos_cm;

        % ---- Save SVG (vector) ----
        if strcmp(y_name, 'ThetaDeg')
            y_tag = 'phi_reg';
        else
            y_tag = y_name;
        end

        base_name = sprintf('%s_vs_%s_partialSpearman_CTRL_BDI_BAI_BIS_80x70mm', y_tag, x_name);
        svg_file  = fullfile(out_dir, [base_name '.svg']);

        try
            exportgraphics(f, svg_file, 'ContentType','vector');
        catch
            print(f, svg_file, '-dsvg');
        end

        close(f);
        plot_col{rowk} = svg_file;
    end
end

%% ---- 7) Save result table (CSV) ------------------------------------
results = table(pred_col, out_col, n_col, r_col, p_col, plot_col, ...
    'VariableNames', {'predictor','outcome','n','spearman_r','spearman_p','plot_file'});

results_export = results(:, {'predictor','outcome','n','spearman_r','spearman_p'});

out_csv = fullfile(out_dir, 'phi_reg_Radius_ProjImpulse_partial_Spearman_CTRL_BDI_BAI_BIS.csv');
writetable(results_export, out_csv);

fprintf('\n✅ Done!\n- Result table: %s\n- Figure directory: %s\n', out_csv, out_dir);
disp(results_export);