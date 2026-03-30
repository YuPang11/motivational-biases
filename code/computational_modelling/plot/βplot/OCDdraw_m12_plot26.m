function OCDdraw_m12_plot26()
%% =========================================
%  M12: beta distribution plot (OCD data, blue theme, 70 x 35 mm)
%  Features:
%  1. Data: read OCD group data (M12_f2_drawsocd.csv)
%  2. Color: switch to OCD blue theme
%  3. Size: total figure width = 70 mm, height = 35 mm
%  4. Layout: shortened y-axis, left-right arrangement
%  5. Key change: right-panel scatter points are constrained
%     within the blue violin shape
% =========================================
close all; clc;

% ---- Edit your folder path here ----
dir_hc1 = 'H:\eegplot\M4b';
% ------------------------------------

% Define OCD file paths
file_f2_csv = fullfile(dir_hc1, 'M12_f2_drawsocd.csv');
file_indiv_csv = '';

% Output file name with identifier
out_tiff = fullfile(dir_hc1, 'M12_beta_distribution_OCD_70x35mm_Constrained2.tif');

% Check file existence
if ~isfile(file_f2_csv)
    error("File not found: %s\nPlease export the CSV file first using the R code.", file_f2_csv);
end

%% ---- 1. Read data ----
fprintf('Reading OCD data...\n');
try
    opts = detectImportOptions(file_f2_csv);
    opts.VariableNamingRule = 'preserve';
    T = readtable(file_f2_csv, opts);
catch
    T = readtable(file_f2_csv, 'PreserveVariableNames', true);
end
coln = T.Properties.VariableNames;

% Get data
beta_group = get_group_draws_from_f2(T, coln);
subj_tbl   = get_subject_medians(file_indiv_csv, T, coln);

%% ---- 2. Global plot settings ----
% Font settings
baseFont = 9;
fontNameStr = "Arial";

% Greek-letter label
xLabelStr = [char(946), '_{', char(952), '-power}'];

% Color settings (OCD blue theme)
faceCol  = [0.25 0.60 0.90];   % OCD blue
lineCol  = [0 0 0];
lw       = 0.8;                % Line width

% --------------------------------------------------------
% Size and margins (keep 70 x 35 mm settings)
% --------------------------------------------------------
figW_cm = 7.0;  % Total width = 70 mm
figH_cm = 3.5;  % Total height = 35 mm

% Margin settings
marginL = 0.90; % Left margin
marginR = 0.10; % Right margin
gap     = 0.90; % Gap between panels

% Bottom and top margins
marginB = 1.25;
marginT = 0.45;

% Compute panel width and height
plotW = (figW_cm - marginL - marginR - gap) / 2;
plotH = figH_cm - marginB - marginT;

% Create figure canvas
fig = figure("Color","w","Units","centimeters","Position",[5 5 figW_cm figH_cm], ...
    "PaperPositionMode","auto","InvertHardcopy","off");
set(fig, "PaperUnits","centimeters", ...
    "PaperPosition",[0 0 figW_cm figH_cm], ...
    "PaperSize",[figW_cm figH_cm]);

%% =========================
% Panel 1: Left - group level
% =========================
ax1 = axes('Parent', fig, 'Units', 'centimeters', ...
    'Position', [marginL, marginB, plotW, plotH]);
hold(ax1, "on");

[fg, xg] = ksdensity(beta_group(:), "Function","pdf");

patch(ax1, [xg fliplr(xg)], [fg zeros(size(fg))], faceCol, ...
    "EdgeColor","none", "FaceAlpha",1);
plot(ax1, xg, fg, "Color", lineCol, "LineWidth", lw);

title(ax1, "Group-level", "FontWeight","normal", "FontName",fontNameStr, "FontSize", baseFont);
xlabel(ax1, xLabelStr, "FontName",fontNameStr, "FontSize", baseFont);
ylabel(ax1, "Density", "FontName",fontNameStr, "FontSize", baseFont);

set(ax1, "FontName",fontNameStr, "FontSize", baseFont, ...
    "Box","off", "TickDir","out", "LineWidth", 0.6, ...
    "XColor","k", "YColor","k", ...
    "XAxisLocation", "bottom", "YAxisLocation", "left");

%% =========================
% Panel 2: Right - subject level
% =========================
leftPos2 = marginL + plotW + gap;
ax2 = axes('Parent', fig, 'Units', 'centimeters', ...
    'Position', [leftPos2, marginB, plotW, plotH]);
hold(ax2, "on");

y = subj_tbl.beta_median;
y = y(isfinite(y));

% A. Draw vertical violin (compute shape boundary)
[fy, yy] = ksdensity(y, "Function","pdf");
fy = fy / max(fy);
maxHalfWidth = 0.35;
fy = fy * maxHalfWidth; % fy is now the half-width at each yy level
x0 = 1;

patch(ax2, [x0 - fy, fliplr(x0 + fy)], [yy, fliplr(yy)], faceCol, ...
    "EdgeColor", lineCol, "LineWidth", lw, "FaceAlpha",1);

% B. Draw scatter points (key change: constrain points inside the shape)
rng(1);
markerSize = 14;

% 1. Interpolate the allowed half-width at each y value
allowed_half_width = interp1(yy, fy, y, 'linear', 0);

% 2. Generate random jitter ratios between -1 and 1
raw_jitter_ratio = (rand(size(y)) - 0.5) * 2;

% 3. Set scaling factor for aesthetics
scale_factor = 0.85;

% 4. Compute the final constrained jitter
constrained_jit = raw_jitter_ratio .* allowed_half_width * scale_factor;

scatter(ax2, x0 + constrained_jit, y, markerSize, "o", ...
    "MarkerFaceColor","w", "MarkerEdgeColor","k", "LineWidth",0.5);

title(ax2, "Subject-level", "FontWeight","normal", "FontName",fontNameStr, "FontSize", baseFont);
xlabel(ax2, xLabelStr, "FontName",fontNameStr, "FontSize", baseFont);
ylabel(ax2, "Median estimate", "FontName",fontNameStr, "FontSize", baseFont);

set(ax2, "FontName",fontNameStr, "FontSize", baseFont, ...
    "Box","off", "TickDir","out", "LineWidth", 0.6, ...
    "XColor","k", "YColor","k", ...
    "XLim", [x0-0.5, x0+0.5], "XTick", [], ...
    "YTickMode", "auto", "YTickLabelMode", "auto", ...
    "XAxisLocation", "bottom", "YAxisLocation", "left");

%% ---- 3. Export figure (TIFF, 1000 dpi) ----
print(fig, out_tiff, "-dtiff", "-r1000");
fprintf("Plot completed!\n");
fprintf("TIFF file saved (1000 dpi): %s\n", out_tiff);

end

%% ==========================================================
% Subfunctions
% ==========================================================
function beta = get_group_draws_from_f2(T, coln)
    candidates = {'X[6]','X_6','X6','X.6.','X(6)'};
    idx = find(ismember(coln, candidates), 1);
    if isempty(idx)
        idx = find(~cellfun(@isempty, regexp(coln, '^(X\[\s*6\s*\]|X[_\.]?6)$', 'once')), 1);
    end
    if isempty(idx)
        idx = find(contains(lower(string(coln)), "beta"), 1);
    end
    if isempty(idx)
        error("Cannot find group-level beta in f2 column names.");
    end
    beta = T{:, idx};
    beta = beta(:);
end

function subj_tbl = get_subject_medians(file_indiv_csv, T, coln)
    if ~isempty(file_indiv_csv) && isfile(file_indiv_csv)
        try
            opts = detectImportOptions(file_indiv_csv);
            opts.VariableNamingRule = 'preserve';
            df = readtable(file_indiv_csv, opts);
        catch
            df = readtable(file_indiv_csv, 'PreserveVariableNames', true);
        end
    end

    x6mask = ~cellfun(@isempty, regexp(coln, '^x_?6\[\d+\]$', 'once'));
    x6_cols = coln(x6mask);
    if isempty(x6_cols)
         error("Cannot find subject-level beta data in f2.");
    end

    idxnum = cellfun(@(s) str2double(regexprep(s,'^.*\[(\d+)\].*$','$1')), x6_cols);
    [~, ord] = sort(idxnum);
    x6_cols = x6_cols(ord);

    med = nan(numel(x6_cols),1);
    for i = 1:numel(x6_cols)
        col_data = T{:, strcmp(coln, x6_cols{i})};
        med(i) = median(col_data, "omitnan");
    end

    med = med(:);
    subj_tbl = table((1:numel(med))', med, 'VariableNames', {'subject', 'beta_median'});
    fprintf("Extracted medians for %d subjects.\n", height(subj_tbl));
end