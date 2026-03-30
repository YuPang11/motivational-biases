%% Relation midfrontal theta power and Pavlovian-Instrumental conflict (FINAL)
clear; clc; close all;

% =============================
% 0) Global style (Arial 9pt, no bold)
% =============================
set(groot, 'defaultAxesFontName',   'Arial');
set(groot, 'defaultTextFontName',   'Arial');
set(groot, 'defaultAxesFontSize',   9);
set(groot, 'defaultTextFontSize',   9);
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');

% =============================
% 1) Load data
% =============================
thetapow = csvread('F:\PIT\model\EEGmodelHC\HC\results\stan\theta&conflict\thetapower.csv',1,0);
thetapow(thetapow==0) = NaN;

conflicttrials = csvread('F:\PIT\model\EEGmodelHC\HC\results\stan\theta&conflict\conflicttrials.csv',1,0);
conflicttrials = reshape(conflicttrials,37,480);

conflict = csvread('F:\PIT\model\EEGmodelHC\HC\results\stan\theta&conflict\M9_osap_conflict.csv',1,0);

nSub   = size(thetapow,1);
x2plot = -3:0.1:3;

% =============================
% 2) Correlation across all trials (optional print)
% =============================
rho = nan(1,nSub);
for iSub = 1:nSub
    idx = ~isnan(thetapow(iSub,:));
    rho(iSub) = corr(thetapow(iSub,idx)', conflict(iSub,idx)');
end
[~,p,~,stats] = ttest(rho); %#ok<ASGLU>
disp([p stats.tstat])

% =============================
% 3) Color mapping for correlation (R -> color)
% =============================
cols = jet(50);                 % 50 levels
rMin = -0.25; rMax = 0.25;      % clamp range shown in colorbar

% =============================
% 4) Figure size (cm) + leave bottom space so 3-line xlabel won't clip
%    (You asked earlier for axes height 40mm, so we enforce that below)
% =============================
figW_cm = 5.5;      % ~55 mm width
figH_cm = 6.8;      % increase height to fit 3-line xlabel

fig = figure('Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off', ...
    'Renderer','painters');

set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

% =============================
% 5) Axes + plot lines
% =============================
ax = axes('Parent',fig);
hold(ax,'on');
colormap(ax, cols);

P1 = nan(1,nSub);
P2 = nan(1,nSub);

for iSub = 1:nSub
    idx = ~isnan(thetapow(iSub,:)) & logical(conflicttrials(iSub,:));
    if nnz(idx) < 5
        continue
    end

    x  = zscore(conflict(iSub,idx));
    y  = thetapow(iSub,idx);

    P  = polyfit(x, y, 1);
    P1(iSub) = P(1);
    P2(iSub) = P(2);

    yfit = P(1)*x2plot + P(2);

    coef = corr(thetapow(iSub,idx)', conflict(iSub,idx)');
    coef = max(min(coef, rMax), rMin);

    cIdx = round((coef - rMin) / (rMax - rMin) * (size(cols,1)-1)) + 1;
    cIdx = max(min(cIdx, size(cols,1)), 1);

    plot(ax, x2plot, yfit, '-', 'Color', cols(cIdx,:), 'LineWidth', 0.8);
end

% Mean line (black)
mP1 = mean(P1,'omitnan');
mP2 = mean(P2,'omitnan');
plot(ax, x2plot, mP1*x2plot + mP2, 'k', 'LineWidth', 2);

% =============================
% 6) Axes formatting
% =============================
set(ax, 'XLim', [-3 3], ...
        'TickDir','out', ...
        'Box','off', ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'FontWeight','normal', ...
        'XColor','k', ...
        'YColor','k', ...
        'LineWidth',0.8);

ylabel(ax, 'Midfrontal theta power', 'FontName','Arial','FontSize',9,'FontWeight','normal');

% =============================
% 7) Three-line xlabel (with arrows + weak/strong)
%    Use spaces to shift weak/strong if needed.
% =============================
arrowL = char(8592);  % ←
arrowR = char(8594);  % →

weakShift   = '      ';           % add spaces to move "weak conflict" right
gapSpaces   = 14;                 % controls spacing between weak and strong
gapShift    = repmat(' ', 1, gapSpaces);

xlabel(ax, { ...
    'Pavlovian-instrumental conflict (z)', ...
    [arrowL '                             ' arrowR], ...
    [weakShift 'Weak conflict' gapShift 'Strong conflict'] ...
    }, 'FontName','Arial','FontSize',9,'FontWeight','normal','Interpreter','none');

% =============================
% 8) Colorbar (R), keep Arial 9pt normal
% =============================
caxis(ax, [rMin rMax]);
cb = colorbar(ax, 'Ticks', -0.2:0.1:0.2);
cb.Label.String = 'R';
set(cb, 'FontName','Arial', 'FontSize', 9);
set(cb.Label, 'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');

% =============================
% 9) Fix axes physical height (Y-axis length) to 40mm = 4.0cm
%    and ensure bottom space for 3-line xlabel.
% =============================
desiredAxH_cm = 4.0;     % 40 mm
bottom_cm     = 2.3;     % bottom space for three-line xlabel (adjust if needed)

ax.Units = 'centimeters';
axPos = ax.Position;
axPos(2) = bottom_cm;
axPos(4) = desiredAxH_cm;
ax.Position = axPos;

% Sync colorbar height with axes height (so colorbar matches)
cb.Units = 'centimeters';
cbPos = cb.Position;
cbPos(2) = axPos(2);
cbPos(4) = axPos(4);
cb.Position = cbPos;

% =============================
% 10) Final enforcement: Arial 9pt, no bold everywhere
% =============================
set(findall(fig, '-property','FontName'),   'FontName','Arial');
set(findall(fig, '-property','FontSize'),   'FontSize',9);
set(findall(fig, '-property','FontWeight'), 'FontWeight','normal');

% =============================
% 11) Export
% =============================
outBase = 'theta_conflict_final';

% line art (1000 dpi)
print(fig, [outBase '_lineart1000.tif'], '-dtiff', '-r1000');

% halftone backups (optional)
print(fig, [outBase '_halftone300.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500.tif'], '-dtiff', '-r500');






%% Relation midfrontal theta power and Pavlovian-Instrumental conflict (FINAL, OCD)
clear; clc; close all;

% =============================
% 0) Global style (Arial 9pt, no bold)
% =============================
set(groot, 'defaultAxesFontName',   'Arial');
set(groot, 'defaultTextFontName',   'Arial');
set(groot, 'defaultAxesFontSize',   9);
set(groot, 'defaultTextFontSize',   9);
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');

% =============================
% 1) Load data (OCD paths)
% =============================
thetapow = csvread('F:\PIT\model\EEGmodelHC\OCD\results\stan\theta&conflict\thetapower.csv', 1, 0);
thetapow(thetapow==0) = NaN;

conflicttrials = csvread('F:\PIT\model\EEGmodelHC\OCD\results\stan\theta&conflict\conflicttrials.csv', 1, 0);

conflict = csvread('F:\PIT\model\EEGmodelHC\OCD\results\stan\theta&conflict\M5_osap_conflict.csv', 1, 0);

nSub = size(thetapow, 1);

conflicttrials = reshape(conflicttrials, nSub, []);

x2plot = -3:0.1:3;

% =============================
% 2) Correlation across all trials (optional print)
% =============================
rho = nan(1, nSub);
for iSub = 1:nSub
    idx = ~isnan(thetapow(iSub,:));
    if nnz(idx) < 5, continue; end
    rho(iSub) = corr(thetapow(iSub,idx)', conflict(iSub,idx)');
end
[~,p,~,stats] = ttest(rho); %#ok<ASGLU>
disp([p stats.tstat])

% =============================
% 3) Color mapping for correlation (R -> color)
% =============================
cols = jet(50);            % 50 levels
rMin = -0.25; rMax = 0.25; % clamp range shown in colorbar

% =============================
% 4) Figure size (cm) + leave bottom space so 3-line xlabel won't clip
%    (axes height fixed to 40mm below)
% =============================
figW_cm = 5.5;   % ~55 mm width
figH_cm = 6.8;   % enough height for 3-line xlabel

fig = figure('Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'Color','w', ...
    'PaperPositionMode','auto', ...
    'InvertHardcopy','off', ...
    'Renderer','painters');

set(fig, 'PaperUnits','centimeters');
set(fig, 'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig, 'PaperSize',[figW_cm figH_cm]);

% =============================
% 5) Axes + plot lines
% =============================
ax = axes('Parent',fig);
hold(ax,'on');
colormap(ax, cols);

P1 = nan(1, nSub);
P2 = nan(1, nSub);

for iSub = 1:nSub
    idx = ~isnan(thetapow(iSub,:)) & logical(conflicttrials(iSub,:));
    if nnz(idx) < 5
        continue
    end

    x  = zscore(conflict(iSub,idx));
    y  = thetapow(iSub,idx);

    P  = polyfit(x, y, 1);
    P1(iSub) = P(1);
    P2(iSub) = P(2);

    yfit = P(1)*x2plot + P(2);

    coef = corr(thetapow(iSub,idx)', conflict(iSub,idx)');
    coef = max(min(coef, rMax), rMin);

    cIdx = round((coef - rMin) / (rMax - rMin) * (size(cols,1)-1)) + 1;
    cIdx = max(min(cIdx, size(cols,1)), 1);

    plot(ax, x2plot, yfit, '-', 'Color', cols(cIdx,:), 'LineWidth', 0.8);
end

% Mean line (black)
mP1 = mean(P1,'omitnan');
mP2 = mean(P2,'omitnan');
plot(ax, x2plot, mP1*x2plot + mP2, 'k', 'LineWidth', 2);

% =============================
% 6) Axes formatting
% =============================
set(ax, 'XLim', [-3 3], ...
        'TickDir','out', ...
        'Box','off', ...
        'FontName','Arial', ...
        'FontSize', 9, ...
        'FontWeight','normal', ...
        'XColor','k', ...
        'YColor','k', ...
        'LineWidth',0.8);

ylabel(ax, 'Midfrontal theta power', 'FontName','Arial','FontSize',9,'FontWeight','normal');

% =============================
% 7) Three-line xlabel (with arrows + weak/strong)
% =============================
arrowL = char(8592);  % ←
arrowR = char(8594);  % →

weakShift   = '      '; % right shift weak
gapSpaces   = 14;       % spacing between weak and strong
gapShift    = repmat(' ', 1, gapSpaces);

xlabel(ax, { ...
    'Pavlovian-instrumental conflict (z)', ...
    [arrowL '                             ' arrowR], ...
    [weakShift 'Weak conflict' gapShift 'Strong conflict'] ...
    }, 'FontName','Arial','FontSize',9,'FontWeight','normal','Interpreter','none');

% =============================
% 8) Colorbar (R), keep Arial 9pt normal
% =============================
caxis(ax, [rMin rMax]);
cb = colorbar(ax, 'Ticks', -0.2:0.1:0.2);
cb.Label.String = 'R';
set(cb, 'FontName','Arial', 'FontSize', 9);
set(cb.Label, 'FontName','Arial', 'FontSize', 9, 'FontWeight','normal');

% =============================
% 9) Fix axes physical height (Y-axis length) to 40mm = 4.0cm
%    and ensure bottom space for 3-line xlabel.
% =============================
desiredAxH_cm = 4.0;  % 40 mm
bottom_cm     = 2.3;  % bottom space for 3-line xlabel

ax.Units = 'centimeters';
axPos = ax.Position;
axPos(2) = bottom_cm;
axPos(4) = desiredAxH_cm;
ax.Position = axPos;

% Sync colorbar height with axes height
cb.Units = 'centimeters';
cbPos = cb.Position;
cbPos(2) = axPos(2);
cbPos(4) = axPos(4);
cb.Position = cbPos;

% =============================
% 10) Final enforcement: Arial 9pt, no bold everywhere
% =============================
set(findall(fig, '-property','FontName'),   'FontName','Arial');
set(findall(fig, '-property','FontSize'),   'FontSize',9);
set(findall(fig, '-property','FontWeight'), 'FontWeight','normal');

% =============================
% 11) Export
% =============================
outBase = 'theta_conflict_final_OCD';

% line art (1000 dpi)
print(fig, [outBase '_lineart1000.tif'], '-dtiff', '-r1000');

% halftone backups (optional)
print(fig, [outBase '_halftone300.tif'], '-dtiff', '-r300');
print(fig, [outBase '_halftone500.tif'], '-dtiff', '-r500');