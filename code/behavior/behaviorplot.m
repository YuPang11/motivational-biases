%% ============================================================
%  Combined Motivational Go/NoGo Task Behaviour (HC & OCD)
%  [Final Adjustment: Vertical Gap Increased by 2mm]
%  - Layout: 50mm width per panel, 2 rows (HC/OCD)
%  - Vertical Gap: Increased from 3.2 to 3.4 cm (+2mm)
%  - Figure Height: Increased from 13.2 to 13.4 cm (+2mm)
%  - Legend: Positioned with 0.2 offset (unchanged relative pos)
%  - Background: White, 1000 DPI TIFF
% ============================================================

clear; clc; close all;

% --- 1. Basic setup ---
try
    addpath('C:\Users\16055\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\boundedline.m\boundedline');
catch
    warning('boundedline path not found. Please check the path.');
end

% Load data
dataFile = 'H:\mt_behavalluse\EEGsingletrial.mat';
if exist(dataFile, 'file')
    load(dataFile);
else
    error('Data file not found: %s', dataFile);
end
nTrialPerBlock = 60;

% --- 2. Subject lists ---
HC_list = [102 105 106 107 108 109 110 113 115 116 117 118 124 125 130 132 133 134 ...
           136 137 138 139 140 141 142 143 145 146 147 148 150 151 152 154 155 156 157];
OCD_list = [201 207 209 210 211 213 214 216 217 218 219 220 221 224 225 227 231 232 ...
            233 234 238 241 243 245 247 249 251 253 254 255 256 258 259 260 261 262];

% --- 3. Colors and fonts ---
clr.win   = [0    0.60 0.25];
clr.avoid = [0.90 0.20 0.10];
lighten   = @(c,a) c + (1-c).*a;
darken    = @(c,a) c.*(1-a);

c.lineWin     = clr.win;
c.lineAvoid   = clr.avoid;
c.barWinL     = lighten(clr.win,  0.35);
c.barWinD     = darken (clr.win,  0.00);
c.barAvdL     = lighten(clr.avoid,0.35);
c.barAvdD     = darken (clr.avoid,0.00);
c.ngWin       = clr.win;
c.ngAvd       = clr.avoid;

fontSz.text  = 9;
fontSz.name  = 'Arial';
lw = 1.0;

% --- 4. Layout parameters (spacing adjusted) ---
figW = 19.0;
% [Adjustment 1] Increase total height to allow a larger row gap (13.2 -> 13.4)
figH = 13.4;

margL = 1.0;
margB = 1.2;
gapX  = 0.7;
pltW = 5.3;
pltH = 4.2;

% [Adjustment 2] Increase vertical gap by 2 mm (3.2 -> 3.4)
gapY = 3.4;

yPos_OCD = margB;
yPos_HC  = margB + pltH + gapY;

% Create figure canvas
f = figure('Units', 'centimeters', 'Position', [5 5 figW figH], 'Color', 'w');

% --- 5. Main plotting loop ---
groups = {'HC', 'OCD'};

for row = 1:2
    groupName = groups{row};
    if row == 1
        yPos = yPos_HC;
        curList = HC_list;
    else
        yPos = yPos_OCD;
        curList = OCD_list;
    end
    
    % === Data summary ===
    idx = find(ismember(par.subjList, curList));
    nSub = numel(idx);
    
    pGo2Win=nan(nSub,nTrialPerBlock); pGo2Avoid=pGo2Win; pNoGo2Win=pGo2Win; pNoGo2Avoid=pGo2Win;
    pCorrectGo2Win=pGo2Win; pCorrectGo2Avoid=pGo2Win;
    
    for ii=1:nSub
        si = idx(ii);
        pGo2Win(ii,:)     = nanmean([data.go(si,data.stim(si,:)==1); data.go(si,data.stim(si,:)==2)]);
        pGo2Avoid(ii,:)   = nanmean([data.go(si,data.stim(si,:)==3); data.go(si,data.stim(si,:)==4)]);
        pNoGo2Win(ii,:)   = nanmean([data.go(si,data.stim(si,:)==5); data.go(si,data.stim(si,:)==6)]);
        pNoGo2Avoid(ii,:) = nanmean([data.go(si,data.stim(si,:)==7); data.go(si,data.stim(si,:)==8)]);
        pCorrectGo2Win(ii,:)   = nanmean([data.acc(si,data.stim(si,:)==1); data.acc(si,data.stim(si,:)==2)]);
        pCorrectGo2Avoid(ii,:) = nanmean([data.acc(si,data.stim(si,:)==3); data.acc(si,data.stim(si,:)==4)]);
    end
    pIncorrectGo2Win = pGo2Win - pCorrectGo2Win;
    pIncorrectGo2Avoid = pGo2Avoid - pCorrectGo2Avoid;
    
    mGoCorr    = mean(mean(pCorrectGo2Win,2));
    mAvoidCorr = mean(mean(pCorrectGo2Avoid,2));
    mGoInc     = mean(mean(pGo2Win,2)) - mGoCorr;
    mAvoidInc  = mean(mean(pGo2Avoid,2)) - mAvoidCorr;
    mNoGoWin   = mean(mean(pNoGo2Win,2));
    mNoGoAvoid = mean(mean(pNoGo2Avoid,2));

    % === Column 1: Trial-by-trial ===
    xPos1 = margL;
    ax1 = axes('Units','centimeters', 'Position', [xPos1 yPos pltW pltH], 'Color', 'none'); hold(ax1, 'on');
    
    boundedline(ax1, 1:nTrialPerBlock, mean(pGo2Win),   std(pGo2Win)/sqrt(nSub),   'cmap', c.lineWin,   'alpha');
    boundedline(ax1, 1:nTrialPerBlock, mean(pGo2Avoid), std(pGo2Avoid)/sqrt(nSub), 'cmap', c.lineAvoid, 'alpha');
    hN1 = boundedline(ax1, 1:nTrialPerBlock, mean(pNoGo2Win),   std(pNoGo2Win)/sqrt(nSub),   'cmap', c.ngWin,   'alpha');
    hN2 = boundedline(ax1, 1:nTrialPerBlock, mean(pNoGo2Avoid), std(pNoGo2Avoid)/sqrt(nSub), 'cmap', c.ngAvd, 'alpha');
    
    if isstruct(hN1) && isfield(hN1,'mainLine')
        set([hN1.mainLine hN2.mainLine], 'LineStyle', '--', 'LineWidth', lw);
    else
        set([hN1 hN2], 'LineStyle', '--', 'LineWidth', lw);
    end
    
    ylabel(ax1, 'P(Go)', 'FontSize', fontSz.text);
    xlabel(ax1, 'Trial', 'FontSize', fontSz.text);
    set(ax1, 'YLim', [0 1], 'XLim', [0 nTrialPerBlock], 'XTick', 0:20:60, 'YTick', 0:0.2:1);

    % === Column 2: Go cues ===
    xPos2 = margL + pltW + gapX;
    ax2 = axes('Units','centimeters', 'Position', [xPos2 yPos pltW pltH], 'Color', 'none'); hold(ax2, 'on');
    
    boundedline(ax2, 1:nTrialPerBlock, mean(pCorrectGo2Win),   std(pCorrectGo2Win)/sqrt(nSub),   'cmap', clr.win,   'alpha');
    boundedline(ax2, 1:nTrialPerBlock, mean(pCorrectGo2Avoid), std(pCorrectGo2Avoid)/sqrt(nSub), 'cmap', clr.avoid, 'alpha');
    boundedline(ax2, 1:nTrialPerBlock, mean(pIncorrectGo2Win),   std(pIncorrectGo2Win)/sqrt(nSub),   'cmap', clr.win,   'alpha');
    boundedline(ax2, 1:nTrialPerBlock, mean(pIncorrectGo2Avoid), std(pIncorrectGo2Avoid)/sqrt(nSub), 'cmap', clr.avoid, 'alpha');
    
    ylabel(ax2, 'P(response)', 'FontSize', fontSz.text);
    xlabel(ax2, 'Trial', 'FontSize', fontSz.text);
    set(ax2, 'YLim', [0 1], 'XLim', [0 nTrialPerBlock], 'XTick', 0:20:60, 'YTick', 0:0.2:1);
    
    text(ax2, 5, 0.92, 'Correct', 'FontSize', fontSz.text, 'Color', 'k');
    text(ax2, 5, 0.10, 'Incorrect', 'FontSize', fontSz.text, 'Color', 'k');

    % === Column 3: Task effects ===
    xPos3 = margL + (pltW + gapX)*2;
    ax3 = axes('Units','centimeters', 'Position', [xPos3 yPos pltW pltH], 'Color', 'none'); hold(ax3, 'on');
    bw = 0.6;
    
    % Draw bars
    bar(ax3, 1, mGoCorr, bw, 'FaceColor', c.barWinL, 'EdgeColor', 'k', 'LineWidth', 1);
    patch(ax3, [1-bw/2 1+bw/2 1+bw/2 1-bw/2], [mGoCorr mGoCorr mGoCorr+mGoInc mGoCorr+mGoInc], c.barWinD, 'EdgeColor', 'k', 'LineWidth', 1);
    bar(ax3, 2, mAvoidCorr, bw, 'FaceColor', c.barAvdL, 'EdgeColor', 'k', 'LineWidth', 1);
    patch(ax3, [2-bw/2 2+bw/2 2+bw/2 2-bw/2], [mAvoidCorr mAvoidCorr mAvoidCorr+mAvoidInc mAvoidCorr+mAvoidInc], c.barAvdD, 'EdgeColor', 'k', 'LineWidth', 1);
    bar(ax3, 3.5, mNoGoWin, bw, 'FaceColor', c.ngWin, 'EdgeColor', 'k', 'LineWidth', 1);
    bar(ax3, 4.5, mNoGoAvoid, bw, 'FaceColor', c.ngAvd, 'EdgeColor', 'k', 'LineWidth', 1);
    
    % Double error bars
    y_go_all   = mean(mean((pGo2Win + pGo2Avoid)./2, 2));
    err_go_all = std(mean(pGo2Win - pGo2Avoid, 2)) / sqrt(nSub);
    errorbar(ax3, 1.5, y_go_all, err_go_all, 'k', 'linestyle', 'none', 'CapSize', 4, 'LineWidth', 1);
    
    y_go_corr   = mean(mean((pCorrectGo2Win + pCorrectGo2Avoid)./2, 2));
    err_go_corr = std(mean(pCorrectGo2Win - pCorrectGo2Avoid, 2)) / sqrt(nSub);
    errorbar(ax3, 1.5, y_go_corr, err_go_corr, 'k', 'linestyle', 'none', 'CapSize', 4, 'LineWidth', 1);
    
    y_nogo_all   = mean(mean((pNoGo2Win + pNoGo2Avoid)./2, 2));
    err_nogo_all = std(mean(pNoGo2Win - pNoGo2Avoid, 2)) / sqrt(nSub);
    errorbar(ax3, 4.0, y_nogo_all, err_nogo_all, 'k', 'linestyle', 'none', 'CapSize', 4, 'LineWidth', 1);
    
    ylabel(ax3, 'P(Go)', 'FontSize', fontSz.text);
    xlabel(ax3, 'Required response', 'FontSize', fontSz.text);
    set(ax3, 'YLim', [0 1], 'XLim', [0.4 5.1], 'XTick', [1.5 4], 'XTickLabel', {'Go','NoGo'}, 'YTick', 0:0.2:1);
end

% --- 6. Draw legends (auto-centered, keep 0.2 offset) ---
legendH = 2.0;
legendW = 8.0;
legendDownOffset = 0.2;

legendY_center = yPos_OCD + pltH + gapY / 2;
legendY_start = legendY_center - legendH / 2 - legendDownOffset;

% >>> 1. Line legend (left)
axL1 = axes('Units','centimeters', 'Position', [margL legendY_start (pltW*2+gapX) legendH], 'Color', 'none');
axis(axL1, 'off'); hold(axL1, 'on');

L_y = 1.0; L_len = 0.5;
plot(axL1, [0 L_len], [L_y L_y], '-', 'Color', c.lineWin, 'LineWidth', lw);
text(axL1, L_len+0.1, L_y, 'Go-to-Win', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');
offset = 2.2;
plot(axL1, [offset offset+L_len], [L_y L_y], '-', 'Color', c.lineAvoid, 'LineWidth', lw);
text(axL1, offset+L_len+0.1, L_y, 'Go-to-Avoid', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');
offset = 4.6;
plot(axL1, [offset offset+L_len], [L_y L_y], '--', 'Color', c.ngWin, 'LineWidth', lw);
text(axL1, offset+L_len+0.1, L_y, 'NoGo-to-Win', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');
offset = 7.2;
plot(axL1, [offset offset+L_len], [L_y L_y], '--', 'Color', c.ngAvd, 'LineWidth', lw);
text(axL1, offset+L_len+0.1, L_y, 'NoGo-to-Avoid', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');
xlim(axL1, [0 10]); ylim(axL1, [0 2.0]);

% >>> 2. Bar legend (right, aligned to column 3 y-axis)
xPosCol3_start = margL + (pltW + gapX)*2;
legendX_start = xPosCol3_start;

axL2 = axes('Units','centimeters', 'Position', [legendX_start legendY_start legendW legendH], 'Color', 'none');
axis(axL2, 'off'); hold(axL2, 'on');

blkW = 0.8;
blkH = 0.45;
pX = 0.0;
pY1 = 1.1;
pY2 = 0.6;

% Win legend
patch(axL2, [pX pX+blkW/2 pX+blkW/2 pX], [pY1 pY1 pY1+blkH pY1+blkH], c.barWinL, 'EdgeColor','k');
patch(axL2, [pX+blkW/2 pX+blkW pX+blkW pX+blkW/2], [pY1 pY1 pY1+blkH pY1+blkH], c.barWinD, 'EdgeColor','k');
text(axL2, pX+blkW+0.3, pY1+blkH/2, 'Win (correct/incorrect)', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');

% Avoid legend
patch(axL2, [pX pX+blkW/2 pX+blkW/2 pX], [pY2 pY2 pY2+blkH pY2+blkH], c.barAvdL, 'EdgeColor','k');
patch(axL2, [pX+blkW/2 pX+blkW pX+blkW pX+blkW/2], [pY2 pY2 pY2+blkH pY2+blkH], c.barAvdD, 'EdgeColor','k');
text(axL2, pX+blkW+0.3, pY2+blkH/2, 'Avoid (correct/incorrect)', 'FontSize', fontSz.text, 'FontName','Arial', 'VerticalAlignment','middle');

xlim(axL2, [0 8]); ylim(axL2, [0 2.0]);

% --- 7. Export (TIFF) ---
set(findall(f, 'type', 'axes'), 'FontName', fontSz.name, 'FontSize', fontSz.text, 'Box', 'off', 'LineWidth', 0.8);

fprintf('Generating 1000 DPI TIFF (Gap +2mm)...\n');
print(f, 'Final_Submission_GapIncreased.tif', '-dtiffn', '-r1000');

fprintf('Done! Output file: Final_Submission_GapIncreased.tif\n');