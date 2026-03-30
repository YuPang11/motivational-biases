function HCdraw_m12_plot26()
%% =========================================
%  M12: beta distribution plot (70 x 26 mm)
%  M12-specific version adapted from M16 script:
%   1) Input file changed to M12_f2_drawshc.csv
%   2) Output file changed to M12_beta_distribution_70x26mm_Compact.tif
%   3) X-axis label changed to \beta_{\theta}
%   4) Original layout and size settings kept unchanged
%% =========================================
close all; clc;

% -----------------------------
% Paths
% -----------------------------
dir_out     = 'F:\PIT\ISPSplot';
file_f2_csv = fullfile(dir_out, 'M12_f2_drawshc.csv');
out_tiff    = fullfile(dir_out, 'M12_beta_distribution_70x26mm_Compact.tif');

file_indiv_csv = '';
BETA_INDEX = 6;   % In M12, beta is the 6th parameter

if ~isfile(file_f2_csv)
    error("File not found: %s\nPlease check the path or filename.", file_f2_csv);
end

%% ---- 1. Load data ----
fprintf('Reading data...\n');
try
    opts = detectImportOptions(file_f2_csv);
    opts.VariableNamingRule = 'preserve';
    T = readtable(file_f2_csv, opts);
catch
    T = readtable(file_f2_csv, 'PreserveVariableNames', true);
end
coln = T.Properties.VariableNames;

beta_group = get_group_draws_from_f2(T, coln, BETA_INDEX);
subj_tbl   = get_subject_medians(file_indiv_csv, T, coln, BETA_INDEX);

%% ---- 2. Global plot settings ----
baseFont    = 9;
fontNameStr = 'Arial';
xLabelStr   = '\beta_{\theta}';

faceCol  = [242 183 108]/255;
lineCol  = [0 0 0];
lw       = 0.8;

% -----------------------------
% Canvas size (70 mm x 26 mm)
% -----------------------------
figW_cm = 7.0;
figH_cm = 2.6;

% Left and right margins
marginL = 1.15;
marginR = 0.05;
gap     = 0.90;

% Top and bottom margins
marginB = 0.72;
marginT = 0.26;

plotW = (figW_cm - marginL - marginR - gap) / 2;

% ---- Axis height compression ----
rawH_cm   = figH_cm - marginB - marginT;
vertComp  = 0.78;
plotH     = rawH_cm * vertComp;
axY       = marginB + (rawH_cm - plotH)/2;

% -----------------------------
% Global font settings
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');

fig = figure("Color","w","Units","centimeters", ...
    "Position",[5 5 figW_cm figH_cm], ...
    "PaperPositionMode","auto","InvertHardcopy","off");
set(fig, "PaperUnits","centimeters", ...
    "PaperPosition",[0 0 figW_cm figH_cm], ...
    "PaperSize",[figW_cm figH_cm]);

%% =========================
% Panel 1: Left - group level
% =========================
ax1 = axes('Parent', fig, 'Units', 'centimeters', ...
    'Position', [marginL, axY, plotW, plotH]);
hold(ax1, "on");

[fg, xg] = ksdensity(beta_group(:), "Function","pdf");
patch(ax1, [xg fliplr(xg)], [fg zeros(size(fg))], faceCol, ...
    "EdgeColor","none", "FaceAlpha",1);
plot(ax1, xg, fg, "Color", lineCol, "LineWidth", lw);

t1 = title(ax1, "group-level", ...
    "FontWeight","normal", "FontName",fontNameStr, "FontSize", baseFont);
set(t1,'Units','normalized');
t1.Position(2) = 0.98;

xl1 = xlabel(ax1, xLabelStr, ...
    "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal", "Interpreter","tex");
ylabel(ax1, "density", ...
    "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal");

set(ax1, "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal", "Box","off", "TickDir","out", "LineWidth", 0.6, ...
    "XColor","k", "YColor","k", "XAxisLocation", "bottom", "YAxisLocation", "left");

ylim(ax1, [0, max(fg)*1.01]);

% Make left-panel ticks more compact
try
    ax1.XAxis.TickLabelGapOffset = -2;
catch
    try
        ax1.TickLabelGapOffset = -2;
    catch
    end
end

% Move left-panel xlabel downward
set(xl1,'Units','normalized');
xl1.Position(2) = -0.14;

%% =========================
% Panel 2: Right - subject level
% =========================
leftPos2 = marginL + plotW + gap;
ax2 = axes('Parent', fig, 'Units', 'centimeters', ...
    'Position', [leftPos2, axY, plotW, plotH]);
hold(ax2, "on");

y = subj_tbl.beta_median;
y = y(isfinite(y));

[fy, yy] = ksdensity(y, "Function","pdf");
fy = fy / max(fy);
maxHalfWidth = 0.35;
fy = fy * maxHalfWidth;
x0 = 1;

patch(ax2, [x0 - fy, fliplr(x0 + fy)], [yy, fliplr(yy)], faceCol, ...
    "EdgeColor", lineCol, "LineWidth", lw, "FaceAlpha",1);

rng(1);
jitWidth   = 0.32;
markerSize = 4;
jit = (rand(size(y)) - 0.5) * jitWidth;

scatter(ax2, x0 + jit, y, markerSize, "o", ...
    "MarkerFaceColor","w", "MarkerEdgeColor","k", "LineWidth",0.35);

t2 = title(ax2, "subject-level", ...
    "FontWeight","normal", "FontName",fontNameStr, "FontSize", baseFont);
set(t2,'Units','normalized');
t2.Position(2) = 0.98;

xlabel(ax2, xLabelStr, ...
    "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal", "Interpreter","tex");
yl2 = ylabel(ax2, "median estimate", ...
    "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal");

set(ax2, "FontName",fontNameStr, "FontSize", baseFont, ...
    "FontWeight","normal", "Box","off", "TickDir","out", "LineWidth", 0.6, ...
    "XColor","k", "YColor","k", ...
    "XLim", [x0-0.5, x0+0.5], ...
    "XAxisLocation", "bottom", "YAxisLocation", "left");

% Keep X-axis line but remove ticks and tick labels
set(ax2, 'XTick', [], 'XTickLabel', []);
set(ax2, 'XColor', 'k');

% "median estimate": shift downward by 1 mm and slightly rightward
dy_norm = 0.1 / plotH;    % 0.1 cm = 1 mm
dx_cm   = 0.06;           % 0.06 cm = 0.6 mm
dx_norm = dx_cm / plotW;

yl2.Units = 'normalized';
yl2.Position(2) = yl2.Position(2) - dy_norm;
yl2.Position(1) = yl2.Position(1) + dx_norm;

q   = quantile(y, [0.05 0.95]);
pad = 0.03 * (q(2) - q(1) + eps);
ylim(ax2, [q(1)-pad, q(2)+pad]);

%% ---- Final safeguard: enforce Arial 9pt, non-bold ----
set(findall(fig, '-property','FontName'),   'FontName','Arial');
set(findall(fig, '-property','FontSize'),   'FontSize', 9);
set(findall(fig, '-property','FontWeight'), 'FontWeight','normal');

%% ---- Export TIFF at 1000 dpi ----
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end
print(fig, out_tiff, "-dtiff", "-r1000");

fprintf("Plotting completed!\n");
fprintf("Saved to: %s\n", out_tiff);
end

%% ==========================================================
% Subfunctions
%% ==========================================================
function beta = get_group_draws_from_f2(T, coln, beta_idx)
    candidates = {
        sprintf('X[%d]', beta_idx), sprintf('X[ %d ]', beta_idx), ...
        sprintf('X_%d', beta_idx), sprintf('X%d', beta_idx), ...
        sprintf('X.%d.', beta_idx), sprintf('X(%d)', beta_idx)
    };
    idx = find(ismember(coln, candidates), 1);

    if isempty(idx)
        pat = sprintf('^(X\\[\\s*%d\\s*\\]|X[_\\.]?%d)$', beta_idx, beta_idx);
        idx = find(~cellfun(@isempty, regexp(coln, pat, 'once')), 1);
    end

    if isempty(idx)
        cname = lower(string(coln));
        hit = find(contains(cname,"beta") | contains(cname,"theta"), 1);
        if ~isempty(hit)
            idx = hit;
        end
    end

    if isempty(idx)
        error("Cannot find group-level beta in f2 column names. Please check BETA_INDEX or CSV column names.");
    end

    beta = T{:, idx};
    beta = beta(:);
end

function subj_tbl = get_subject_medians(file_indiv_csv, T, coln, beta_idx)
    %#ok<NASGU>
    if ~isempty(file_indiv_csv) && isfile(file_indiv_csv)
        try
            opts = detectImportOptions(file_indiv_csv);
            opts.VariableNamingRule = 'preserve';
            readtable(file_indiv_csv, opts);
        catch
            readtable(file_indiv_csv, 'PreserveVariableNames', true);
        end
    end

    pat = sprintf('^x_?%d\\[\\d+\\]$', beta_idx);
    xmask  = ~cellfun(@isempty, regexp(lower(coln), pat, 'once'));
    x_cols = coln(xmask);

    if isempty(x_cols)
        error("Cannot find subject-level beta in f2 (e.g., x_%d[i] or x%d[i]). Please check BETA_INDEX or CSV column names.", beta_idx, beta_idx);
    end

    idxnum = cellfun(@(s) str2double(regexprep(s,'^.*\[(\d+)\].*$','$1')), x_cols);
    [~, ord] = sort(idxnum);
    x_cols = x_cols(ord);

    med = nan(numel(x_cols),1);
    for i = 1:numel(x_cols)
        col_data = T{:, strcmp(coln, x_cols{i})};
        med(i) = median(col_data, "omitnan");
    end

    subj_tbl = table((1:numel(med))', med(:), ...
        'VariableNames', {'subject', 'beta_median'});
    fprintf("Extracted medians for %d subjects.\n", height(subj_tbl));
end