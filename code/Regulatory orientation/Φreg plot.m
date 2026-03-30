%% ================================================================
%  HC: M16 beta + M19 beta (MATLAB Version) - PUBLICATION READY
%  Metric calculation & plotting (W 60 mm × H 70 mm)
%  Output: F:/PIT/ISPSplot/
%
%  Revisions:
%   1) Figure size: W=60 mm (6.0 cm), H=70 mm (7.0 cm)
%   2) Larger white triangle for Mean
%   3) Mean legend at lower right
%   4) Force Arial 9 pt + normal for all objects (including legend)
%   5) Put axes/ticks on top to avoid contourf overlap
%  ================================================================
clc; clear; close all;

%% ---- 0) Paths and settings ----------------------------------------
data_dir = 'F:/PIT/model/EEGmodelHC/HC/results/stan/';
out_dir  = 'F:/PIT/ISPSplot/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Load CSV data exported from R (same as your runnable version)
try
    draws16  = readmatrix(fullfile(data_dir, 'M16_draws.csv'));
    drawsAlt = readmatrix(fullfile(data_dir, 'ALT_draws.csv'));
    fid = fopen(fullfile(data_dir, 'model_type_flag.txt'), 'r');
    flag_str = fgetl(fid);
    fclose(fid);
    is_M19 = strcmpi(strtrim(flag_str), 'TRUE');
catch
    error('Please export the CSV data in R first.');
end

nSub = size(draws16, 2);

%% ---- Global font fallback: Arial 9 pt + normal --------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultLegendFontName', 'Arial');
set(groot, 'defaultLegendFontSize', 9);
set(groot, 'defaultLegendFontWeight', 'normal');

%% ---- 1) Metric calculation ----------------------------------------
calc_metrics = @(d) struct(...
    'mean',  mean(d, 1)', ...
    'sd',    std(d, 0, 1)', ...
    'p_neg', mean(d < 0, 1)', ...
    'p_pos', mean(d > 0, 1)');

stats16  = calc_metrics(draws16);
statsAlt = calc_metrics(drawsAlt);

% Build table
T = table();
T.sub        = (1:nSub)';
T.beta16     = stats16.mean;
T.sd_beta16  = stats16.sd;
T.p16_neg    = stats16.p_neg;
T.p16_pos    = stats16.p_pos;

T.betaALT    = statsAlt.mean;
T.sd_betaALT = statsAlt.sd;
T.pALT_neg   = statsAlt.p_neg;
T.pALT_pos   = statsAlt.p_pos;

%% ---- 2) Composite metrics -----------------------------------------
z = @(x) (x - mean(x)) ./ std(x);

T.PCI = -T.beta16;
if is_M19
    T.QCI = -T.betaALT;
else
    T.QCI =  T.betaALT;
end

T.PCI_z     = z(T.PCI);
T.QCI_z     = z(T.QCI);
T.H_NGtW    = T.PCI_z + T.QCI_z;
T.H_GtA     = T.PCI_z - T.QCI_z;
T.Balance   = T.PCI_z - T.QCI_z;
T.Opponency = -T.beta16 .* T.betaALT;
T.R_mag     = sqrt(T.PCI_z.^2 + T.QCI_z.^2);
T.ThetaDeg  = atan2(T.QCI_z, T.PCI_z) * 180 / pi;
T.w         = 1 ./ (T.sd_beta16.^2 + T.sd_betaALT.^2);

writetable(T, fullfile(out_dir, 'subject_metrics_M16_M19_matlab.csv'));

%% ---- 3) Common plotting settings ----------------------------------
% Figure size: W 60 mm, H 70 mm
figW_cm = 6.0;
figH_cm = 7.0;

fontName = 'Arial';
fontSize = 9;

% Color settings
hexColors = {'FFF3E0', 'FFE0B2', 'FFCC80', 'FFB74D', 'FFA726', 'FB8C00', 'EF6C00'};
cmap = zeros(length(hexColors), 3);
for i = 1:length(hexColors)
    cmap(i, :) = sscanf(hexColors{i}, '%2x%2x%2x', [1 3]) / 255;
end

brks = [0.10, 0.25, 0.40, 0.60, 0.80, 0.95];

% Axis labels (TeX subscripts)
x_lab = [char(946) '_{ISPS-motor} (Pavlovian bias)'];
y_lab = [char(946) '_{ISPS-lPFC} (Instrumental action values)'];

y_data   = T.betaALT;
g_beta16 = mean(T.beta16);
g_y      = mean(y_data);

%% ---- 4) Plot KDE scatter ------------------------------------------
f1 = figure('Units', 'centimeters', 'Position', [5 5 figW_cm figH_cm], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], 'PaperSize', [figW_cm figH_cm], ...
    'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

ax1 = axes('Parent', f1);
hold(ax1, 'on');

% --- 4.1 Compute 2D KDE ---
gridPts = 200;
x_lims = [-3.5, 1.5];
y_lims = [-2, 1];

[Xgrid, Ygrid] = meshgrid(linspace(x_lims(1), x_lims(2), gridPts), ...
                          linspace(y_lims(1), y_lims(2), gridPts));
data = [T.beta16, y_data];
[f, ~] = ksdensity(data, [Xgrid(:), Ygrid(:)]);
F = reshape(f, gridPts, gridPts);
F_norm = F / max(F(:));

% --- 4.2 Filled contour plot ---
contourf(ax1, Xgrid, Ygrid, F_norm, brks, 'LineStyle', 'none');
colormap(ax1, cmap);
caxis(ax1, [0 1]);

% --- 4.3 Reference lines ---
xline(ax1, 0, '--k', 'LineWidth', 0.5);
yline(ax1, 0, '--k', 'LineWidth', 0.5);

% --- 4.4 Scatter points (subjects) ---
scatter(ax1, T.beta16, y_data, 6, 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.3 0.3 0.3], ...
    'MarkerFaceAlpha', 0.9, ...
    'LineWidth', 0.5);

% --- 4.5 Mean point (white triangle) ---
% Revision (2): larger white triangle
meanMarkerSize = 7;   % originally 4; use 8/9 if you want it larger
p_mean = plot(ax1, g_beta16, g_y, '^', ...
    'MarkerSize', meanMarkerSize, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 0.8);

% --- 4.6 Axis styling ---
xlabel(ax1, x_lab, 'Interpreter','tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');
ylabel(ax1, y_lab, 'Interpreter','tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

xlim(ax1, x_lims);
ylim(ax1, y_lims);
xticks(ax1, -3:1:1);
yticks(ax1, -2:1:1);

set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, ...
    'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

axis(ax1, 'square');

% Key fix: put axes/ticks on top so contourf does not cover them
set(ax1, 'Layer', 'top');

% --- 4.7 Legend: lower right ---
% Revision (3): legend at lower right + 9 pt + normal
lgd = legend(ax1, p_mean, {'Mean'}, 'Location', 'southeast', 'Box', 'off');
% ---- Manually move legend 2 mm to the right ----
lgd.Units = 'centimeters';
pos = lgd.Position;      % [x y w h], unit = cm (relative to figure)
pos(1) = pos(1) + 0.6;   % move right by 0.2 cm = 2 mm
lgd.Position = pos;

lgd.FontName   = 'Arial';
lgd.FontSize   = 9;
lgd.FontWeight = 'normal';
lgd.ItemTokenSize = [10, 10];

% --- 4.8 Final fallback: force Arial 9 pt + normal for all objects in f1 ---
set(findall(f1, '-property','FontName'),   'FontName','Arial');
set(findall(f1, '-property','FontSize'),   'FontSize',9);
set(findall(f1, '-property','FontWeight'), 'FontWeight','normal');

% --- 4.9 Export ---
out_base1 = fullfile(out_dir, 'HC_beta_scatter_60x70mm');
print(f1, [out_base1 '_halftone.tif'], '-dtiff', '-r500');
print(f1, [out_base1 '_lineart.tif'],  '-dtiff', '-r1000');

%% ---- 5) Plot Theta histogram --------------------------------------
f2 = figure('Units', 'centimeters', 'Position', [12 5 figW_cm figH_cm], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], 'PaperSize', [figW_cm figH_cm], ...
    'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

ax2 = axes('Parent', f2);
histogram(ax2, T.ThetaDeg, 18, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'w');
hold(ax2, 'on');

xline(ax2, 0,  '--k', 'LineWidth', 0.6);
xline(ax2, 90, '--k', 'LineWidth', 0.6);

xlabel(ax2, 'ThetaDeg (\circ)', 'Interpreter', 'tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');
ylabel(ax2, 'Count', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, ...
    'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

% --- 5.1 Final fallback: force Arial 9 pt + normal for all objects in f2 ---
set(findall(f2, '-property','FontName'),   'FontName','Arial');
set(findall(f2, '-property','FontSize'),   'FontSize',9);
set(findall(f2, '-property','FontWeight'), 'FontWeight','normal');

% --- 5.2 Export ---
out_base2 = fullfile(out_dir, 'HC_theta_hist_60x70mm');
print(f2, [out_base2 '_halftone.tif'], '-dtiff', '-r500');
print(f2, [out_base2 '_lineart.tif'],  '-dtiff', '-r1000');

fprintf('MATLAB processing finished.\nOutput path: %s\n', out_dir);
fprintf('- Size: W 60 mm × H 70 mm\n- Font: Arial 9 pt (including legend) + normal\n- Larger Mean triangle + legend at lower right\n- Quality: 500 dpi & 1000 dpi TIFF\n');





%% ================================================================
%  OCD: M16 beta + M21 beta (MATLAB Version) - PUBLICATION READY
%  KDE contours + Scatter (blue palette)
%  Size: W 60 mm × H 70 mm
%  Output: F:/PIT/model/EEGmodelHC/OCD/plots/betaR_matlab
%
%  Key revisions (same as HC version):
%   1) Figure size: W=60 mm (6.0 cm), H=70 mm (7.0 cm)
%   2) Larger white triangle for Mean
%   3) Mean legend at lower right (optional: move 2 mm right)
%   4) Force Arial 9 pt + normal for all objects (including legend)
%   5) Put axes/ticks on top to avoid contourf overlap
%   6) Change y-axis label to: β_{ISPS-motor} (instrumental action value)
%  ================================================================
clc; clear; close all;

%% ---- 0) Paths and settings ----------------------------------------
data_dir = 'F:/PIT/model/EEGmodelHC/OCD/results/stan/';
out_dir  = 'F:/PIT/model/EEGmodelHC/OCD/plots/betaR_matlab/';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Load OCD CSV data
try
    draws16 = readmatrix(fullfile(data_dir, 'M16_draws_OCD.csv'));
    draws21 = readmatrix(fullfile(data_dir, 'M21_draws_OCD.csv'));
catch
    error('Please run the R script first to export M16_draws_OCD.csv and M21_draws_OCD.csv.');
end

% Check whether the subject counts match
if size(draws16, 2) ~= size(draws21, 2)
    error('The subject counts in M16 and M21 do not match!');
end
nSub = size(draws16, 2);

%% ---- Global font fallback: Arial 9 pt + normal --------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultLegendFontName', 'Arial');
set(groot, 'defaultLegendFontSize', 9);
set(groot, 'defaultLegendFontWeight', 'normal');

%% ---- 1) Metric calculation ----------------------------------------
calc_metrics = @(d) struct(...
    'mean',  mean(d, 1)', ...
    'sd',    std(d, 0, 1)', ...
    'p_neg', mean(d < 0, 1)', ...
    'p_pos', mean(d > 0, 1)');

stats16 = calc_metrics(draws16);
stats21 = calc_metrics(draws21);

T = table();
T.sub       = (1:nSub)';
T.beta16    = stats16.mean;
T.sd_beta16 = stats16.sd;
T.p16_neg   = stats16.p_neg;
T.p16_pos   = stats16.p_pos;

T.beta21    = stats21.mean;
T.sd_beta21 = stats21.sd;
T.p21_neg   = stats21.p_neg;
T.p21_pos   = stats21.p_pos;

%% ---- 2) Composite metrics -----------------------------------------
z = @(x) (x - mean(x)) ./ std(x);

% Direction definition (same as in your original OCD script)
T.PCI = -T.beta16;
T.QCI = -T.beta21;

T.PCI_z     = z(T.PCI);
T.QCI_z     = z(T.QCI);
T.H_NGtW    = T.PCI_z + T.QCI_z;
T.H_GtA     = T.PCI_z - T.QCI_z;
T.Balance   = T.PCI_z - T.QCI_z;
T.Opponency = -T.beta16 .* T.beta21;
T.R_mag     = sqrt(T.PCI_z.^2 + T.QCI_z.^2);
T.ThetaDeg  = atan2(T.QCI_z, T.PCI_z) * 180 / pi;

T.Type1_prob     = T.p16_neg .* T.p21_neg;
T.Type2_prob     = T.p16_pos .* T.p21_pos;
T.TypeContrast   = T.Type1_prob - T.Type2_prob;

writetable(T, fullfile(out_dir, 'subject_metrics_M16_M21_matlab.csv'));

%% ---- 3) Common plotting settings ----------------------------------
% Figure size: W 60 mm, H 70 mm
figW_cm = 6.0;
figH_cm = 7.0;

fontName = 'Arial';
fontSize = 9;

% Blue palette (same as your original script)
hexColors = {'EFF3FB','DCE9F5','C6DBEF','9ECAE1','6BAED6','3182BD','08519C'};
cmap = zeros(numel(hexColors), 3);
for i = 1:numel(hexColors)
    cmap(i, :) = sscanf(hexColors{i}, '%2x%2x%2x', [1 3]) / 255;
end

% KDE contour levels (to match HC style, use 0.10~0.95)
brks = [0.10, 0.25, 0.40, 0.60, 0.80, 0.95];

% Axis labels (TeX subscripts; y-axis updated as requested)
x_lab = [char(946) '_{ISPS-motor} (Pavlovian bias)'];
y_lab = [char(946) '_{ISPS-motor} (Instrumental action values)'];

% Data to plot (x=beta16, y=beta21)
x_data = T.beta16;
y_data = T.beta21;

g_x = mean(x_data);
g_y = mean(y_data);

% Axis limits with a small margin to avoid clipping
pad_range = @(v, mult) [min(v) - range(v)*mult, max(v) + range(v)*mult];
x_lims = pad_range(x_data, 0.08);
y_lims = pad_range(y_data, 0.08);

%% ---- 4) Plot KDE scatter ------------------------------------------
f1 = figure('Units', 'centimeters', 'Position', [5 5 figW_cm figH_cm], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], 'PaperSize', [figW_cm figH_cm], ...
    'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

ax1 = axes('Parent', f1);
hold(ax1, 'on');

% --- 4.1 2D KDE ---
gridPts = 200;
[Xgrid, Ygrid] = meshgrid(linspace(x_lims(1), x_lims(2), gridPts), ...
                          linspace(y_lims(1), y_lims(2), gridPts));
data2 = [x_data, y_data];
[f, ~] = ksdensity(data2, [Xgrid(:), Ygrid(:)]);
F = reshape(f, gridPts, gridPts);
F_norm = F / max(F(:));

% --- 4.2 Filled contours ---
contourf(ax1, Xgrid, Ygrid, F_norm, brks, 'LineStyle', 'none');
colormap(ax1, cmap);
caxis(ax1, [0 1]);

% --- 4.3 Reference lines ---
xline(ax1, 0, '--k', 'LineWidth', 0.5);
yline(ax1, 0, '--k', 'LineWidth', 0.5);

% --- 4.4 Scatter points (small points for compact figure) ---
scatter(ax1, x_data, y_data, 6, 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.3 0.3 0.3], ...
    'MarkerFaceAlpha', 0.9, ...
    'LineWidth', 0.5);

% --- 4.5 Mean point (larger white triangle) ---
meanMarkerSize = 7;
p_mean = plot(ax1, g_x, g_y, '^', ...
    'MarkerSize', meanMarkerSize, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 0.8);

% --- 4.6 Axis styling ---
xlabel(ax1, x_lab, 'Interpreter','tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');
ylabel(ax1, y_lab, 'Interpreter','tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

xlim(ax1, x_lims);
ylim(ax1, y_lims);

set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, ...
    'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

axis(ax1, 'square');
set(ax1, 'Layer', 'top'); % Put axes/ticks on top

% --- 4.7 Legend: lower right ---
lgd = legend(ax1, p_mean, {'Mean'}, 'Location', 'southeast', 'Box', 'off');
lgd.FontName   = 'Arial';
lgd.FontSize   = 9;
lgd.FontWeight = 'normal';
lgd.ItemTokenSize = [10, 10];

% (Optional) move legend 2 mm to the right: uncomment if needed
% lgd.Units = 'centimeters';
% pos = lgd.Position;
% pos(1) = pos(1) + 0.2;   % 0.2 cm = 2 mm
% lgd.Position = pos;

% --- 4.8 Final fallback: force Arial 9 pt + normal for all objects in f1 ---
set(findall(f1, '-property','FontName'),   'FontName','Arial');
set(findall(f1, '-property','FontSize'),   'FontSize',9);
set(findall(f1, '-property','FontWeight'), 'FontWeight','normal');

% --- 4.9 Export (500/1000 dpi TIFF) ---
out_base1 = fullfile(out_dir, 'OCD_beta16_vs_beta21_60x70mm');
print(f1, [out_base1 '_halftone.tif'], '-dtiff', '-r500');
print(f1, [out_base1 '_lineart.tif'],  '-dtiff', '-r1000');

%% ---- 5) Theta histogram (same 60×70 mm) ---------------------------
f2 = figure('Units', 'centimeters', 'Position', [12 5 figW_cm figH_cm], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], 'PaperSize', [figW_cm figH_cm], ...
    'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

ax2 = axes('Parent', f2);
histogram(ax2, T.ThetaDeg, 18, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'w');
hold(ax2, 'on');

xline(ax2, 0,  '--k', 'LineWidth', 0.6);
xline(ax2, 90, '--k', 'LineWidth', 0.6);

xlabel(ax2, 'ThetaDeg (\circ)', 'Interpreter', 'tex', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');
ylabel(ax2, 'Count', 'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, ...
    'FontName', fontName, 'FontSize', fontSize, 'FontWeight','normal');

set(findall(f2, '-property','FontName'),   'FontName','Arial');
set(findall(f2, '-property','FontSize'),   'FontSize',9);
set(findall(f2, '-property','FontWeight'), 'FontWeight','normal');

out_base2 = fullfile(out_dir, 'OCD_theta_hist_60x70mm');
print(f2, [out_base2 '_halftone.tif'], '-dtiff', '-r500');
print(f2, [out_base2 '_lineart.tif'],  '-dtiff', '-r1000');

fprintf('OCD MATLAB processing finished. Results saved to: %s\n', out_dir);
fprintf('- Size: W 60 mm × H 70 mm\n- Font: Arial 9 pt (including legend) + normal\n');