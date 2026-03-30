%% cue
%% Two-condition plot: permutation results (Arial 9pt + 50 mm + no legend + 1000 dpi)
% Required variables (already available):
% TF1, TF2, MFchans, iFreq, freq.time, iTime, MFtime

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Compute time course for each subject
% -----------------------------
nSub = numel(TF1);
nT   = numel(freq.time);

plot_congruent   = nan(nSub, nT);
plot_incongruent = nan(nSub, nT);

for iSub = 1:nSub
    plot_congruent(iSub, :)   = squeeze(mean(TF1{iSub}.powspctrm(MFchans, iFreq, :), [1, 2]));
    plot_incongruent(iSub, :) = squeeze(mean(TF2{iSub}.powspctrm(MFchans, iFreq, :), [1, 2]));

    % ---- Fallback (for older MATLAB versions, use this instead of the two lines above) ----
    % tmp1 = TF1{iSub}.powspctrm(MFchans, iFreq, :);
    % tmp2 = TF2{iSub}.powspctrm(MFchans, iFreq, :);
    % plot_congruent(iSub, :)   = squeeze(mean(mean(tmp1,1,'omitnan'),2,'omitnan'));
    % plot_incongruent(iSub, :) = squeeze(mean(mean(tmp2,1,'omitnan'),2,'omitnan'));
end

m_cong = mean(plot_congruent,   1, 'omitnan');
m_inco = mean(plot_incongruent, 1, 'omitnan');

% Significant time-window indices (compatible with both logical mask and index vector)
if islogical(iTime)
    idxSig = find(iTime);
else
    idxSig = iTime;
end

% -----------------------------
% 2) Create figure: width = 50 mm (5 cm)
% -----------------------------
figW_cm = 5.0;   % 50 mm
figH_cm = 4.0;   % Adjust height if needed for journal layout

fig = figure('Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep exported physical size stable
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

ax = axes('Parent', fig);
hold(ax, 'on');

% -----------------------------
% 3) Unified line width and colors
% -----------------------------
colCong = [0 0 1];   % Blue
colInco = [1 0 0];   % Red

lwMain   = 1.2;      % Main curves
lwSigWin = 2.0;      % Significant time window: moderately thicker
lwRef    = 0.6;      % Thinner vertical dashed line
lwSigBar = 1.0;      % Significant interval bar: not bold

% -----------------------------
% 4) Main curves
% -----------------------------
plot(ax, freq.time, m_cong(:), 'Color', colCong, 'LineWidth', lwMain);
plot(ax, freq.time, m_inco(:), 'Color', colInco, 'LineWidth', lwMain);

% -----------------------------
% 5) Significant time window: moderately thicker (overlay only on significant interval)
% -----------------------------
if ~isempty(idxSig)
    plot(ax, freq.time(idxSig), m_cong(idxSig), 'Color', colCong, 'LineWidth', lwSigWin);
    plot(ax, freq.time(idxSig), m_inco(idxSig), 'Color', colInco, 'LineWidth', lwSigWin);
end

% -----------------------------
% 6) Axes and labels
% -----------------------------
ylabel(ax, 'Power (dB)', 'FontName','Arial','FontSize',9);
xlabel(ax, 'Time (s)',   'FontName','Arial','FontSize',9);   % Added x-axis title

box(ax, 'off');

set(ax, 'XTick', [0 0.5 1], ...
        'XTickLabel', {'0','0.5','1'}, ...
        'XLim', [-0.25 1.3], ...
        'YLim', [-0.25 3], ...
        'FontName', 'Arial', ...
        'FontSize', 9, ...
        'TickDir', 'out');

% Vertical dashed line (thinner)
plot(ax, [0.700 0.700], get(ax,'YLim'), ':k', 'LineWidth', lwRef);%HC plot(ax, [0.679 0.679], get(ax,'YLim'), ':k', 'LineWidth', lwRef);

% Significant interval bar (not bold)
ySig = 2.3;
plot(ax, MFtime, [ySig ySig], 'k', 'LineWidth', lwSigBar);

% Significance stars
text(ax, mean(MFtime), ySig + 0.2, '**', ...
    'FontName','Arial', 'FontSize',9, 'FontWeight','bold', ...
    'HorizontalAlignment','center');

% -----------------------------
% 7) Final safeguard: ensure all text is Arial 9pt
% -----------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

% -----------------------------
% 8) Export (1000 dpi recommended for line plots)
% -----------------------------
outBase = 'theta_timecourse_50mm';
print(fig, [outBase '.tif'], '-dtiff', '-r1000');

% For halftone standard (300–500 dpi), if needed:
% print(fig, [outBase '_300dpi.tif'], '-dtiff', '-r300');
% print(fig, [outBase '_500dpi.tif'], '-dtiff', '-r500');







%% Four-condition plot: line plot of separate conditions (Arial 9pt + 50 mm width + NO legend + DPI export)

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Data (keep original calculation)
% -----------------------------
chan2plot = [6 11 127];

plot_gtw  = squeeze(mean(mean(mean(gtw_fin(:,chan2plot,iFreq,:),1),2),3));
plot_gta  = squeeze(mean(mean(mean(gta_fin(:,chan2plot,iFreq,:),1),2),3));
plot_ngtw = squeeze(mean(mean(mean(ngtw_fin(:,chan2plot,iFreq,:),1),2),3));
plot_ngta = squeeze(mean(mean(mean(ngta_fin(:,chan2plot,iFreq,:),1),2),3));

% -----------------------------
% 2) Create figure: width = 50 mm (5 cm)
% -----------------------------
figW_cm = 5.0;   % 50 mm
figH_cm = 4.0;   % Adjustable

fig = figure('Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep physical size stable on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

ax = axes('Parent', fig);
hold(ax, 'on');

% -----------------------------
% 3) Line styles (keep original colors; significant bar not bold)
% -----------------------------
lwMain = 1.0;   % Four main curves
lwRef  = 0.6;   % Vertical dashed line
lwSig  = 1.0;   % Significant interval bar (not bold: recommended 0.8~1.2)

colG = [0 .5 0];
colR = [.8 0 0];

plot(ax, freq.time, plot_gtw,  '-',  'LineWidth', lwMain, 'Color', colG);
plot(ax, freq.time, plot_gta,  '-',  'LineWidth', lwMain, 'Color', colR);
plot(ax, freq.time, plot_ngtw, '--', 'LineWidth', lwMain, 'Color', colG);
plot(ax, freq.time, plot_ngta, '--', 'LineWidth', lwMain, 'Color', colR);

% -----------------------------
% 4) Axes and labels (Arial 9pt)
% -----------------------------
ylim(ax, [-1 3]);
set(ax, 'XTick', [0 .5 1], ...
        'XLim',  [-.25 1.3], ...
        'YTick', 0:3, ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'TickDir','out');

% Vertical dashed line (more stable after ylim is set)
plot(ax, [.679 .679], get(ax,'YLim'), ':k', 'LineWidth', lwRef);

% Significant interval bar + stars (bar not bold)
plot(ax, MFtime, [2.5 2.5], 'k', 'LineWidth', lwSig);
text(ax, mean(MFtime), 2.6, '**', ...
    'FontName','Arial', 'FontSize', 9, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

ylabel(ax, 'Power (dB)', 'FontName','Arial', 'FontSize',9);
xlabel(ax, 'Time (s)',         'FontName','Arial', 'FontSize',9);

box(ax, 'off');

% -----------------------------
% 5) Remove legend (as requested)
% -----------------------------
% legend(ax, ... )   % No legend
% legend boxoff

% -----------------------------
% 6) Final safeguard: ensure all text is Arial 9pt
% -----------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

% -----------------------------
% 7) Export (journal dpi requirements)
% -----------------------------
outBase = 'theta_timecourse_50mm_conditions';

% Line plot: 1000 dpi (recommended for this figure)
print(fig, [outBase '_lineart1000OCD.tif'], '-dtiff', '-r1000');

% Backup: halftone 300–500 dpi
print(fig, [outBase '_halftone300OCD.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500OCD.tif'], '-dtiff', '-r500');





%% TF map: difference TF map (Arial 9pt + 50 mm width + 40 mm height + DPI export)

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Data
% -----------------------------
plot_differ = squeeze(mean(fre_differ(chan2plot,:,:), 1));  % freq x time

% -----------------------------
% 2) Create figure: width = 50 mm (5 cm) × height = 40 mm (4 cm)
% -----------------------------
figW_cm = 5.0;   % 50 mm
figH_cm = 4.0;   % 40 mm

fig = figure('Name','midfrontal TF power difference', ...
    'Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep physical size stable on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

ax = axes('Parent', fig);
hold(ax, 'on');

% -----------------------------
% 3) Main plot: contourf
% -----------------------------
colormap(ax, jet);
contourf(ax, freq.time, freq.freq, plot_differ, 36, 'LineStyle', 'none');

% -----------------------------
% 4) Axes settings (with black axis lines)
% -----------------------------
set(ax, 'YTick', [2 4 8 16 32], ...
        'YScale','log', ...
        'XTick', [0 .5 1], ...
        'YLim',  [1.5 50], ...
        'XLim',  [-.25 1.3], ...
        'CLim',  [-0.5 0.5], ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'TickDir','out', ...
        'XColor','k', ...
        'YColor','k', ...
        'LineWidth',0.9);

box(ax, 'off');   % Keep only left/bottom axes (common in papers)

% -----------------------------
% 5) Analysis window box
% -----------------------------
lwBox = 0.8;
plot(ax, [MFtime(2) MFtime(2)], MFfreq, 'k', 'LineWidth', lwBox);
plot(ax, [MFtime(1) MFtime(1)], MFfreq, 'k', 'LineWidth', lwBox);
plot(ax, MFtime, [MFfreq(1) MFfreq(1)], 'k', 'LineWidth', lwBox);
plot(ax, MFtime, [MFfreq(2) MFfreq(2)], 'k', 'LineWidth', lwBox);

% -----------------------------
% 6) Colorbar: ticks + dB units (-0.5 dB / 0 dB / 0.5 dB)
% -----------------------------
cb = colorbar(ax);
set(cb, 'FontName','Arial', 'FontSize', 9);

cb.Ticks = [-0.5 0 0.5];
cb.TickLabels = {'-0.5 dB','0 dB','0.5 dB'};

% Optional: colorbar label
cb.Label.String   = 'dB';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9;

% -----------------------------
% 7) Axis labels
% -----------------------------
ylabel(ax, 'Frequency (Hz)', 'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');
xlabel(ax, 'Time (s)',      'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');

% -----------------------------
% 8) Final safeguard: ensure all text is Arial 9pt
% -----------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

% -----------------------------
% 9) Export (dpi)
% -----------------------------
outBase = 'TF_difference_50x40mm';

print(fig, [outBase '_lineart1000OCD.tif'], '-dtiff', '-r1000');
print(fig, [outBase '_halftone300OCD.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500OCD.tif'], '-dtiff', '-r500');




%% Topoplot: Arial 9pt + 30 mm × 30 mm + NO colorbar + small black dots + smaller white circles + thin outline/ears/nose + dpi export

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Colormap (same as TF map)
% -----------------------------
cmap = jet(64);

% -----------------------------
% 2) Build freq_differ structure (keep original logic)
% -----------------------------
freq_differ = struct();
freq_differ.powspctrm = fre_differ;
freq_differ.time      = freq.time;
freq_differ.freq      = freq.freq;
freq_differ.label     = freq.label;
freq_differ.elec      = freq.elec;

% Manually adjust electrode coordinates (keep original logic)
elec = freq.elec;
elec.chanpos(:,1) = -freq.elec.chanpos(:,2);
elec.chanpos(:,2) =  freq.elec.chanpos(:,1);
freq_differ.elec  = elec;

% -----------------------------
% 3) Create figure: width = 30 mm (3 cm) × height = 30 mm (3 cm)
% -----------------------------
figW_cm = 3.0;   % 30 mm
figH_cm = 3.0;   % 30 mm

fig = figure('Name','midfrontal TF power Valence x Required action', ...
    'Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep physical size stable on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

% -----------------------------
% 4) FieldTrip topo settings
% -----------------------------
cfg = [];
cfg.ylim     = MFfreq;
cfg.xlim     = MFtime;
cfg.zlim     = [-0.5 0.5];

cfg.comment  = 'no';
cfg.style    = 'straight';
cfg.elec     = elec;
cfg.colormap = cmap;

% Other electrodes: small black dots
cfg.marker       = 'on';
cfg.markersymbol = '.';
cfg.markercolor  = [0 0 0];
cfg.markersize   = 3;     % Recommended: 3~5

% Highlight electrodes: smaller white hollow circles
cfg.highlight        = 'on';
cfg.highlightsymbol  = 'o';
cfg.highlightcolor   = [1 1 1];   % White edge
cfg.highlightsize    = 2;         % Circle size

% Key: ensure MFchans is highlighted correctly (index -> label)
if exist('MFchans','var') && isnumeric(MFchans)
    cfg.highlightchannel = freq_differ.label(MFchans);
elseif exist('MFchans','var')
    cfg.highlightchannel = MFchans; % Already label
else
    cfg.highlightchannel = freq_differ.label([6 11 127]);
end

% -----------------------------
% 5) Draw topo
% -----------------------------
ft_topoplotTFR(cfg, freq_differ);
ax = gca;

colormap(ax, cmap);
caxis(ax, [-0.5 0.5]);

% Remove colorbar completely
cb = findall(fig,'Type','ColorBar');
if ~isempty(cb), delete(cb); end

% -----------------------------
% 6) Make outline/ears/nose thin
% -----------------------------
outlineLW = 0.5;  % Recommended: 0.5~0.8

hLines = findobj(ax, 'Type','line');
for i = 1:numel(hLines)
    mk = get(hLines(i),'Marker');
    % Only modify lines without markers (usually scalp outline/nose/ears)
    if (ischar(mk) && strcmp(mk,'none')) || isempty(mk)
        set(hLines(i), 'LineWidth', outlineLW, 'Color','k');
    end
end

% -----------------------------
% 7) Make white hollow circles more balanced
% -----------------------------
hlLW = 0.6;  % White circle edge width
hHL = findobj(ax,'Type','line','Marker','o');
for i = 1:numel(hHL)
    set(hHL(i), 'MarkerFaceColor','none', ...
        'LineWidth', hlLW, ...
        'MarkerEdgeColor', [1 1 1]);
end

% -----------------------------
% 8) Final safeguard: ensure Arial 9pt
% -----------------------------
set(findall(fig, '-property','FontName'), 'FontName','Arial');
set(findall(fig, '-property','FontSize'), 'FontSize', 9);

% -----------------------------
% 9) Export (dpi)
% -----------------------------
outBase = 'topo_diff_30x30mm';

% Line art: 1000 dpi
print(fig, [outBase '_lineart1000.tif'], '-dtiff', '-r1000');

% Halftone: 300 / 500 dpi
print(fig, [outBase '_halftone300.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500.tif'], '-dtiff', '-r500');





%% resp2: permutation results plot (Arial 9pt + 50×40 mm + no title/legend + DPI export)
% Required variables: TF1, TF2, MFchans, iFreq, freq.time

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Compute time course for each subject
% -----------------------------
nSub = 36;
nT   = numel(freq.time);

congruent   = nan(nSub, nT);
incongruent = nan(nSub, nT);

for iSub = 1:nSub
    tmp1 = TF1{iSub}.powspctrm(MFchans, iFreq, :);
    tmp2 = TF2{iSub}.powspctrm(MFchans, iFreq, :);

    congruent(iSub,:)   = squeeze(mean(mean(tmp1,1,'omitnan'),2,'omitnan'));
    incongruent(iSub,:) = squeeze(mean(mean(tmp2,1,'omitnan'),2,'omitnan'));
end

m_cong = mean(congruent,   1, 'omitnan');
m_inco = mean(incongruent, 1, 'omitnan');

% Significant time window (based on the specified range)
% iTime = find(freq.time >= -0.676 & freq.time <= -0.250);
iTime = find(freq.time >= -0.550 & freq.time <= -0.126);

% -----------------------------
% 2) Create figure: 50×40 mm (5×4 cm)
% -----------------------------
figW_cm = 5.0;   % 50 mm
figH_cm = 4.0;   % 40 mm

fig = figure('Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Fix physical size on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

ax = axes('Parent', fig);
hold(ax, 'on');

% -----------------------------
% 3) Line styles
% -----------------------------
colCong = [0 0.5 0];    % Green
colInco = [0.8 0 0];    % Red (same as four-condition plot)

lwMain   = 1.2;   % Main curve width
lwSigBar = 1.0;   % Significant bar not bold
lwRef    = 0.6;   % Thin line if needed

% -----------------------------
% 4) Main curves (do NOT bold the significant interval)
% -----------------------------
plot(ax, freq.time, m_cong, 'Color', colCong, 'LineWidth', lwMain);
plot(ax, freq.time, m_inco, 'Color', colInco, 'LineWidth', lwMain);

% -----------------------------
% 5) Significant interval bar + label (keep, but do not bold curves)
% -----------------------------
yBar = 1.85;
if ~isempty(iTime)
    plot(ax, [freq.time(iTime(1)) freq.time(iTime(end))], [yBar yBar], 'k', 'LineWidth', lwSigBar);

    text(ax, mean(freq.time(iTime)), yBar+0.20, '**', ...
        'FontName','Arial','FontSize',9,'FontWeight','normal', ...
        'HorizontalAlignment','center');
end

% -----------------------------
% 6) Axes settings (Arial 9pt)
% -----------------------------
ylabel(ax, 'Power (dB)', 'FontName','Arial', 'FontSize',9, 'FontWeight','normal');

% Added x-axis title
xlabel(ax, 'Time (s)', 'FontName','Arial', 'FontSize',9, 'FontWeight','normal');

box(ax, 'off');
set(ax, 'XTick', [-1 0 0.5], ...
        'YTick', [-0.5 0 1 2], ...
        'YLim',  [-1 2.3], ...
        'XLim',  [-1 0.5], ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'TickDir','out');

% Remove title: no title
% Remove legend: no legend

% -----------------------------
% 7) Final safeguard: ensure all text is Arial 9pt
% -----------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

% -----------------------------
% 8) Export (dpi)
% -----------------------------
outBase = 'perm_lineplot_50x40mm';

% Line art: 1000 dpi
print(fig, [outBase '_lineart1000resp.tif'], '-dtiff', '-r1000');

% Halftone: 300 / 500 dpi
print(fig, [outBase '_halftone300resp.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500resp.tif'], '-dtiff', '-r500');





%% resp TF map: difference TF map (Arial 9pt + 50 mm width + 40 mm height + DPI export)
% X-axis changed to: -1 ~ 0.5

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Data
% -----------------------------
plot_differ = squeeze(mean(fre_differ(chan2plot,:,:), 1));  % freq x time

% -----------------------------
% 2) Create figure: width = 50 mm (5 cm) × height = 40 mm (4 cm)
% -----------------------------
figW_cm = 5.0;   % 50 mm
figH_cm = 4.0;   % 40 mm

fig = figure('Name','midfrontal TF power difference', ...
    'Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep physical size stable on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

ax = axes('Parent', fig);
hold(ax, 'on');

% -----------------------------
% 3) Main plot: contourf
% -----------------------------
colormap(ax, jet);
contourf(ax, freq.time, freq.freq, plot_differ, 36, 'LineStyle', 'none');

% -----------------------------
% 4) Axes settings (x-axis: -1 ~ 0.5)
% -----------------------------
set(ax, 'YTick', [4 8 16 32], ...
        'YScale','log', ...
        'XTick', [-1 -0.5 0 0.5], ...   % You can also change this to [-1 0 0.5] if needed
        'YLim',  [3 50], ...
        'XLim',  [-1 0.5], ...          % Key change
        'CLim',  [-0.5 0.5], ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'TickDir','out', ...
        'XColor','k', ...
        'YColor','k', ...
        'LineWidth',0.9);

box(ax, 'off');   % Keep only left/bottom axes

% -----------------------------
% 5) Analysis window box
% -----------------------------
lwBox = 0.8;
plot(ax, [MFtime(2) MFtime(2)], MFfreq, 'k', 'LineWidth', lwBox);
plot(ax, [MFtime(1) MFtime(1)], MFfreq, 'k', 'LineWidth', lwBox);
plot(ax, MFtime, [MFfreq(1) MFfreq(1)], 'k', 'LineWidth', lwBox);
plot(ax, MFtime, [MFfreq(2) MFfreq(2)], 'k', 'LineWidth', lwBox);

% -----------------------------
% 6) Colorbar: ticks + dB units (-0.5 dB / 0 dB / 0.5 dB)
% -----------------------------
cb = colorbar(ax);
set(cb, 'FontName','Arial', 'FontSize', 9);

cb.Ticks = [-0.5 0 0.5];
cb.TickLabels = {'-0.5 dB','0 dB','0.5 dB'};

% Optional: colorbar label
cb.Label.String   = 'dB';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9;

% -----------------------------
% 7) Axis labels
% -----------------------------
ylabel(ax, 'Frequency (Hz)', 'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');
xlabel(ax, 'Time (s)',      'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');

% -----------------------------
% 8) Final safeguard: ensure all text is Arial 9pt
% -----------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

% -----------------------------
% 9) Export (dpi)
% -----------------------------
outBase = 'TF_difference_50x40mm';

print(fig, [outBase '_lineart1000respOCD.tif'], '-dtiff', '-r1000');
print(fig, [outBase '_halftone300respOCD.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500respOCD.tif'], '-dtiff', '-r500');





%% resp topoplot: Arial 9pt + 30 mm × 30 mm + NO colorbar + small black dots + smaller white circles + thin outline/ears/nose + dpi export

% -----------------------------
% 0) Global font settings (recommended)
% -----------------------------
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultTextFontSize', 9);

% -----------------------------
% 1) Colormap (same as TF map)
% -----------------------------
cmap = jet(64);

% -----------------------------
% 2) Build freq_differ structure (keep original logic)
% -----------------------------
freq_differ = struct();
freq_differ.powspctrm = fre_differ;
freq_differ.time      = freq.time;
freq_differ.freq      = freq.freq;
freq_differ.label     = freq.label;
freq_differ.elec      = freq.elec;

% Manually adjust electrode coordinates (keep original logic)
elec = freq.elec;
elec.chanpos(:,1) = -freq.elec.chanpos(:,2);
elec.chanpos(:,2) =  freq.elec.chanpos(:,1);
freq_differ.elec  = elec;

% -----------------------------
% 3) Create figure: width = 30 mm (3 cm) × height = 30 mm (3 cm)
% -----------------------------
figW_cm = 3.0;   % 30 mm
figH_cm = 3.0;   % 30 mm

fig = figure('Name','midfrontal TF power Valence x Required action', ...
    'Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off');

% Keep physical size stable on export
set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

% -----------------------------
% 4) FieldTrip topo settings
% -----------------------------
cfg = [];
cfg.ylim     = MFfreq;
cfg.xlim     = MFtime;
cfg.zlim     = [-0.5 0.5];

cfg.comment  = 'no';
cfg.style    = 'straight';
cfg.elec     = elec;
cfg.colormap = cmap;

% Other electrodes: small black dots
cfg.marker       = 'on';
cfg.markersymbol = '.';
cfg.markercolor  = [0 0 0];
cfg.markersize   = 3;     % Recommended: 3~5

% Highlight electrodes: smaller white hollow circles
cfg.highlight        = 'on';
cfg.highlightsymbol  = 'o';
cfg.highlightcolor   = [1 1 1];   % White edge
cfg.highlightsize    = 2;         % Circle size

% Key: ensure MFchans is highlighted correctly (index -> label)
if exist('MFchans','var') && isnumeric(MFchans)
    cfg.highlightchannel = freq_differ.label(MFchans);
elseif exist('MFchans','var')
    cfg.highlightchannel = MFchans; % Already label
else
    cfg.highlightchannel = freq_differ.label([6 11 127]);
end

% -----------------------------
% 5) Draw topo
% -----------------------------
ft_topoplotTFR(cfg, freq_differ);
ax = gca;

colormap(ax, cmap);
caxis(ax, [-0.5 0.5]);

% Remove colorbar completely
cb = findall(fig,'Type','ColorBar');
if ~isempty(cb), delete(cb); end

% -----------------------------
% 6) Make outline/ears/nose thin
% -----------------------------
outlineLW = 0.6;  % Recommended: 0.5~0.8

hLines = findobj(ax, 'Type','line');
for i = 1:numel(hLines)
    mk = get(hLines(i),'Marker');
    % Only modify lines without markers (usually scalp outline/nose/ears)
    if (ischar(mk) && strcmp(mk,'none')) || isempty(mk)
        set(hLines(i), 'LineWidth', outlineLW, 'Color','k');
    end
end

% -----------------------------
% 7) Make white hollow circles more balanced
% -----------------------------
hlLW = 0.6;  % White circle edge width
hHL = findobj(ax,'Type','line','Marker','o');
for i = 1:numel(hHL)
    set(hHL(i), 'MarkerFaceColor','none', ...
        'LineWidth', hlLW, ...
        'MarkerEdgeColor', [1 1 1]);
end

% -----------------------------
% 8) Final safeguard: ensure Arial 9pt
% -----------------------------
set(findall(fig, '-property','FontName'), 'FontName','Arial');
set(findall(fig, '-property','FontSize'), 'FontSize', 9);

% -----------------------------
% 9) Export (dpi)
% -----------------------------
outBase = 'topo_diff_30x30mm';

% Line art: 1000 dpi
print(fig, [outBase '_lineart1000respOCD.tif'], '-dtiff', '-r1000');

% Halftone: 300 / 500 dpi
print(fig, [outBase '_halftone300respOCD.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500respOCD.tif'], '-dtiff', '-r500');
