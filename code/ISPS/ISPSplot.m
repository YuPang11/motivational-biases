%% HC
%% =========================
%  A) Midfrontal–Lateral PFC  [NO exclusion]
%  Go = mean(contra, ipsi);  NoGo = nogo
%  width = 50 mm, Arial 9 pt (no bold), save as SVG
%% =========================

% If results_dir is not defined, use the current folder:
% results_dir = pwd;

% -------- 1) Organize data to plot: TF_pfc2plot(sub, valence, resp[Go/NoGo]) --------
TF_pfc2plot = nan(size(TF_pfc,1), 2, 2);
TF_pfc2plot(:,:,1) = nanmean(TF_pfc(:,:,1:2), 3);  % Go (avg contra & ipsi)
TF_pfc2plot(:,:,2) = TF_pfc(:,:,3);                % NoGo

% -------- 2) Mean + SEM --------
Xwin   = squeeze(TF_pfc2plot(:,1,:));   % sub × 2
Xavoid = squeeze(TF_pfc2plot(:,2,:));   % sub × 2

mean_win   = nanmean(Xwin,   1);
sd_win     = nanstd (Xwin,   0, 1);
n_win      = sum(~isnan(Xwin),1);
sem_win    = sd_win ./ sqrt(n_win);

mean_avoid = nanmean(Xavoid, 1);
sd_avoid   = nanstd (Xavoid, 0, 1);
n_avoid    = sum(~isnan(Xavoid),1);
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% -------- 3) Plot: 50 mm width + Arial 9 pt + no bold --------
figW_mm = 50; figH_mm = 50;
figW_cm = figW_mm/10; figH_cm = figH_mm/10;

fig = figure('Units','centimeters', 'Position',[2 2 figW_cm figH_cm], ...
             'Color','w', 'Renderer','painters');

h1 = bar(1:2,       mean_win,   0.4, 'FaceColor',[0 0.6 0]); hold on; % Win
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'FaceColor',[0.7 0 0]);         % Avoid

errorbar(1:2,       mean_win,   sem_win,   'k','LineStyle','none','LineWidth',0.8);
errorbar(1.45:2.45, mean_avoid, sem_avoid, 'k','LineStyle','none','LineWidth',0.8);

ax = gca;
set(ax,'XTick',1.23:2.23,'XTickLabel',{'Go','NoGo'}, ...
       'XLim',[.7 2.7], 'YLim',[0 25], 'YTick',0:10:20, ...
       'TickDir','out','LineWidth',0.8, ...
       'FontName','Arial','FontSize',9,'FontWeight','normal');
box off;

ylabel('ISPS (% change)', 'FontName','Arial','FontSize',9,'FontWeight','normal');
% lgd = legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');
% set(lgd,'FontName','Arial','FontSize',9,'FontWeight','normal','Box','off');

title('Lateral PFC', 'FontName','Arial','FontSize',9,'FontWeight','bold', 'Color', [0.8 0.4 1]);

% -------- 4) Final pass: force Arial 9 pt + normal across the full figure --------
set(findall(fig,'-property','FontName'),   'FontName','Arial');
set(findall(fig,'-property','FontSize'),   'FontSize',9);
set(findall(fig,'-property','FontWeight'), 'FontWeight','normal');

% -------- 5) Export SVG --------
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig,'PaperSize',[figW_cm figH_cm]);

out_svg = fullfile(results_dir, 'ICPC_Midfrontal_LateralPFC_GoNoGo_NOexcl2.svg');
print(fig, out_svg, '-dsvg', '-painters');

fprintf('Saved SVG: %s\n', out_svg);



%% =========================
%  B) Midfrontal–Motor (exclude contralateral)  [NO exclusion]
%  Plot response 2:3 = Ipsi(Go) & NoGo
%  width = 50 mm, Arial 9 pt (no bold), save as SVG
%% =========================

% If results_dir is not defined, use the current folder:
% results_dir = pwd;

% -------- 1) Mean + SEM (Win vs Avoid; Go=ipsi, NoGo=nogo) --------
Xwin   = squeeze(TF_par(:,1,2:3));   % sub × 2
Xavoid = squeeze(TF_par(:,2,2:3));   % sub × 2

mean_win   = nanmean(Xwin,   1);
sd_win     = nanstd (Xwin,   0, 1);
n_win      = sum(~isnan(Xwin),1);
sem_win    = sd_win ./ sqrt(n_win);

mean_avoid = nanmean(Xavoid, 1);
sd_avoid   = nanstd (Xavoid, 0, 1);
n_avoid    = sum(~isnan(Xavoid),1);
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% -------- 2) Plot: 50 mm width + Arial 9 pt + no bold --------
figW_mm = 50; figH_mm = 50;
figW_cm = figW_mm/10; figH_cm = figH_mm/10;

fig = figure('Units','centimeters', 'Position',[2 2 figW_cm figH_cm], ...
             'Color','w', 'Renderer','painters');

h1 = bar(1:2,       mean_win,   0.4, 'FaceColor',[0 0.6 0]); hold on;
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'FaceColor',[0.7 0 0]);

errorbar(1:2,       mean_win,   sem_win,   'k','LineStyle','none','LineWidth',0.8);
errorbar(1.45:2.45, mean_avoid, sem_avoid, 'k','LineStyle','none','LineWidth',0.8);

ax = gca;
set(ax,'XTick',1.23:2.23,'XTickLabel',{'Go','NoGo'}, ...
       'XLim',[.7 2.7], 'YLim',[0 25], 'YTick',0:10:20, ...
       'TickDir','out','LineWidth',0.8, ...
       'FontName','Arial','FontSize',9,'FontWeight','normal');
box off;

ylabel('ISPS (% change)', 'FontName','Arial','FontSize',9,'FontWeight','normal');
% lgd = legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');
% set(lgd,'FontName','Arial','FontSize',9,'FontWeight','normal','Box','off');

title('Motor', 'FontName','Arial','FontSize',9,'FontWeight','bold', 'Color', [0.4 0.7 1]);

% -------- 3) Final pass: force Arial 9 pt + normal across the full figure --------
set(findall(fig,'-property','FontName'),   'FontName','Arial');
set(findall(fig,'-property','FontSize'),   'FontSize',9);
set(findall(fig,'-property','FontWeight'), 'FontWeight','normal');

% -------- 4) Export SVG --------
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig,'PaperSize',[figW_cm figH_cm]);

out_svg = fullfile(results_dir, 'ICPC_Midfrontal_Motor_GoNoGo_NOexcl2.svg');
print(fig, out_svg, '-dsvg', '-painters');

fprintf('Saved SVG: %s\n', out_svg);





%% OCDdelete6
%% =========================
%  A) Midfrontal–Lateral PFC  (exclude subjects)
%  width = 50 mm, Arial 9 pt (no bold), save as SVG
%% =========================

% -------- 0) Save path (delete this line if results_dir is already defined above) --------
% results_dir = pwd;

% -------- 1) Exclude subjects (by row index) --------
sub_excl = [1 10 18 21 23 31];

TF_pfc_excl = TF_pfc;
TF_pfc_excl(sub_excl,:,:) = NaN;

% -------- 2) Go/NoGo data (Go=mean(contra,ipsi), NoGo=nogo) --------
TF_pfc2plot = nan(size(TF_pfc_excl,1), 2, 2); % sub × val(win/avoid) × resp(Go/NoGo)
TF_pfc2plot(:,:,1) = nanmean(TF_pfc_excl(:,:,1:2), 3);  % Go
TF_pfc2plot(:,:,2) = TF_pfc_excl(:,:,3);                % NoGo

% -------- 3) Mean + SEM (automatically ignore NaN) --------
mean_win   = squeeze(nanmean(TF_pfc2plot(:,1,:), 1));                 % 1×2
sd_win     = squeeze(nanstd (TF_pfc2plot(:,1,:), 0, 1));              % 1×2
n_win      = squeeze(sum(~isnan(TF_pfc2plot(:,1,:)), 1));             % 1×2
sem_win    = sd_win ./ sqrt(n_win);

mean_avoid = squeeze(nanmean(TF_pfc2plot(:,2,:), 1));
sd_avoid   = squeeze(nanstd (TF_pfc2plot(:,2,:), 0, 1));
n_avoid    = squeeze(sum(~isnan(TF_pfc2plot(:,2,:)), 1));
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% -------- 4) Plot: 50 mm width + Arial 9 pt + no bold --------
figW_mm = 50; figH_mm = 50;
figW_cm = figW_mm/10; figH_cm = figH_mm/10;

fig = figure('Units','centimeters', 'Position',[2 2 figW_cm figH_cm], ...
             'Color','w', 'Renderer','painters');

h1 = bar(1:2,       mean_win,   0.4, 'FaceColor',[0 0.6 0]); hold on;
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'FaceColor',[0.7 0 0]);

errorbar(1:2,       mean_win,   sem_win,   'k', 'LineStyle','none', 'LineWidth',0.8);
errorbar(1.45:2.45, mean_avoid, sem_avoid, 'k', 'LineStyle','none', 'LineWidth',0.8);

ax = gca;
set(ax, 'XTick',1.23:2.23, 'XTickLabel',{'Go','NoGo'}, ...
        'XLim',[.7 2.7], 'YLim',[0 25], 'YTick',0:10:20, ...
        'TickDir','out', 'LineWidth',0.8, ...
        'FontName','Arial', 'FontSize',9, 'FontWeight','normal');
box off;

ylabel('ISPS (% change)', 'FontName','Arial', 'FontSize',9, 'FontWeight','normal');
% lgd = legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');  % legend
% set(lgd, 'FontName','Arial', 'FontSize',9, 'FontWeight','normal', 'Box','off'); % legend

title('Lateral PFC', 'FontName','Arial','FontSize',9,'FontWeight','bold', 'Color', [0.8 0.4 1]);
% —— Final pass: force unified font / size / no bold across the full figure —— 
set(findall(fig,'-property','FontName'),   'FontName','Arial');
set(findall(fig,'-property','FontSize'),   'FontSize',9);
set(findall(fig,'-property','FontWeight'), 'FontWeight','normal');

% -------- 5) Export SVG (vector) --------
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig,'PaperSize',[figW_cm figH_cm]);

out_svg = fullfile(results_dir, 'ICPC_Midfrontal_LateralPFC_GoNoGo.svg');
print(fig, out_svg, '-dsvg', '-painters');

% (Optional) If you also need 1000 dpi line-art PNG:
% out_png = fullfile(results_dir, 'ICPC_Midfrontal_LateralPFC_GoNoGo_1000dpi.png');
% print(fig, out_png, '-dpng', '-r1000');

fprintf('Saved SVG: %s\n', out_svg);




%% =========================
%  B) Midfrontal–Motor (exclude contralateral)  (exclude subjects)
%  width = 50 mm, Arial 9 pt (no bold), save as SVG
%% =========================

% -------- 0) Save path (delete this line if results_dir is already defined above) --------
% results_dir = pwd;

% -------- 1) Exclude subjects (by row index) --------
sub_excl = [1 10 18 21 23 31];

TF_par_excl = TF_par;
TF_par_excl(sub_excl,:,:) = NaN;

% -------- 2) Mean + SEM (Go=ipsi; NoGo=nogo, i.e. response 2:3) --------
Xwin   = TF_par_excl(:,1,2:3);   % sub × 2
Xavoid = TF_par_excl(:,2,2:3);

mean_win  = squeeze(nanmean(Xwin,   1));                 % 1×2
sd_win    = squeeze(nanstd (Xwin,   0, 1));              % 1×2
n_win     = squeeze(sum(~isnan(Xwin),1));                % 1×2
sem_win   = sd_win ./ sqrt(n_win);

mean_avoid = squeeze(nanmean(Xavoid, 1));
sd_avoid   = squeeze(nanstd (Xavoid, 0, 1));
n_avoid    = squeeze(sum(~isnan(Xavoid),1));
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% -------- 3) Plot: 50 mm width + Arial 9 pt + no bold --------
figW_mm = 50; figH_mm = 50;
figW_cm = figW_mm/10; figH_cm = figH_mm/10;

fig = figure('Units','centimeters', 'Position',[2 2 figW_cm figH_cm], ...
             'Color','w', 'Renderer','painters');

h1 = bar(1:2,       mean_win,   0.4, 'FaceColor',[0 0.6 0]); hold on;
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'FaceColor',[0.7 0 0]);

errorbar(1:2,       mean_win,   sem_win,   'k', 'LineStyle','none', 'LineWidth',0.8);
errorbar(1.45:2.45, mean_avoid, sem_avoid, 'k', 'LineStyle','none', 'LineWidth',0.8);

ax = gca;
set(ax, 'XTick',1.23:2.23, 'XTickLabel',{'Go','NoGo'}, ...
        'XLim',[.7 2.7], 'YLim',[0 25], 'YTick',0:10:20, ...
        'TickDir','out', 'LineWidth',0.8, ...
        'FontName','Arial', 'FontSize',9, 'FontWeight','normal');
box off;

ylabel('ISPS (% change)', 'FontName','Arial', 'FontSize',9, 'FontWeight','normal');
% lgd = legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');
% set(lgd, 'FontName','Arial', 'FontSize',9, 'FontWeight','normal', 'Box','off');

title('Motor', 'FontName','Arial','FontSize',9,'FontWeight','bold', 'Color', [0.4 0.7 1]);

% —— Final pass: force unified font / size / no bold across the full figure ——
set(findall(fig,'-property','FontName'),   'FontName','Arial');
set(findall(fig,'-property','FontSize'),   'FontSize',9);
set(findall(fig,'-property','FontWeight'), 'FontWeight','normal');

% -------- 4) Export SVG (vector) --------
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',[0 0 figW_cm figH_cm]);
set(fig,'PaperSize',[figW_cm figH_cm]);

out_svg = fullfile(results_dir, 'ICPC_Midfrontal_Motor_GoNoGo.svg');
print(fig, out_svg, '-dsvg', '-painters');

% (Optional) 1000 dpi line-art PNG:
% out_png = fullfile(results_dir, 'ICPC_Midfrontal_Motor_GoNoGo_1000dpi.png');
% print(fig, out_png, '-dpng', '-r1000');

fprintf('Saved SVG: %s\n', out_svg);





%% top
%% HC Step 5: Compute condition effects (final revised version: smaller ROI dots + 50x50 mm + colorbar moved upward)
%% ================================================================
clc; clear; close all;

%% 1) Load data
load('H:\MT_all(afterqujizhi)\cueoriandbu\correct\TOP\HC\CONGRUENT\tICPC_all.mat', 'tICPC','parTF');  
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); 

fprintf('tICPC size = %s (sub × cond × chan × freq × time)\n', mat2str(size(tICPC)));

% Axes
labels    = freq.label(:);          
freq_axis = parTF.freq4tf(:)';      
time_axis = parTF.toi4tf(:)';       

% Analysis time window
ANALYSIS_WIN = [0.45 0.7];

%% 2) EGI → FieldTrip planar rotation & layout
elec = freq.elec;
x0 = elec.chanpos(:,1); y0 = elec.chanpos(:,2);
elec.chanpos(:,1) = -y0; elec.chanpos(:,2) =  x0;
if isfield(elec,'elecpos') && ~isempty(elec.elecpos)
    x0 = elec.elecpos(:,1); y0 = elec.elecpos(:,2);
    elec.elecpos(:,1) = -y0; elec.elecpos(:,2) =  x0;
end
cfgLay      = []; cfgLay.elec = elec;
lay         = ft_prepare_layout(cfgLay);

% Align layout
[~, ia, ib] = intersect(lay.label, labels, 'stable');
lay.label   = lay.label(ia);
lay.pos     = lay.pos(ia,:);
lay.width   = lay.width(ia);
lay.height  = lay.height(ia);
labels      = labels(ib);
tICPC       = tICPC(:,:,ib,:,:); 

%% 3) ROI channels
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2  3 121 122];
lChan    = [36 42 41 47];
rChan    = [92 97 101 102];

%% 4) Data structure & computation
ICPC            = struct();
ICPC.label      = labels;
ICPC.freq       = freq_axis;
ICPC.time       = time_axis;
ICPC.dimord     = 'chan_freq_time';

% Condition definition
CON   = [1 2 6];
INCON = [3 4 5];
LEFT  = [1 3];
RIGHT = [2 4];

% Compute (mask motor)
dat = tICPC;                              
dat(:,LEFT,  rChan,:,:) = NaN;            
dat(:,RIGHT, lChan,:,:) = NaN;            
ICPC.powspctrm = squeeze( nanmean( nanmean( dat(:,INCON,:,:,:) - dat(:,CON,:,:,:), 2), 1) );

%% ================================================================
% Plotting
% ================================================================

% 1. Set canvas (50 mm x 50 mm)
figW_cm = 5.0;
figH_cm = 5.0;

figure('name','ICPC_Congruency', ...
    'Units', 'centimeters', ...
    'Position', [10 10 figW_cm figH_cm], ...
    'Color', 'w', ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], ...
    'PaperSize', [figW_cm figH_cm]);

% 2. Topoplot configuration
cfg = [];
cfg.layout       = lay;
cfg.parameter    = 'powspctrm';
cfg.comment      = 'no';
cfg.style        = 'straight';
cfg.xlim         = ANALYSIS_WIN;
cfg.ylim         = [freq_axis(1) freq_axis(end)];
cfg.zlim         = [-30 30];       
cfg.colormap     = jet(64);

% Tiny black background dots (already set to minimum = 1)
cfg.marker       = 'on';
cfg.markersymbol = '.';
cfg.markercolor  = [0 0 0];
cfg.markersize   = 1;             
cfg.highlight    = 'off';         

% 3. Draw topography
ft_topoplotTFR(cfg, ICPC);

% 4. Adjust axis position
ax = gca;
set(ax, 'Position', [0.02, 0.20, 0.96, 0.80]); 

hold on;

% 5. Manually draw ROI scatter points
% ★★★ Change: reduce size from 15 to 8 ★★★
sz = 8;  
lw = 0.5; 

% Lateral PFC (dark red)
scatter(lay.pos(lpfcChan,1), lay.pos(lpfcChan,2), sz, [0.8 0.4 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);
scatter(lay.pos(rpfcChan,1), lay.pos(rpfcChan,2), sz, [0.8 0.4 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% Motor (dark blue)
scatter(lay.pos(lChan,1), lay.pos(lChan,2), sz, [0.4 0.7 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);
scatter(lay.pos(rChan,1), lay.pos(rChan,2), sz, [0.4 0.7 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% Midfrontal (white)
scatter(lay.pos(MFchan,1), lay.pos(MFchan,2), sz, [1 1 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% 6. Refine outline lines
hLines = findobj(ax, 'Type','line');
for i = 1:numel(hLines)
    mk = get(hLines(i),'Marker');
    if (ischar(mk) && strcmp(mk,'none')) || isempty(mk)
        set(hLines(i), 'LineWidth', 0.5, 'Color','k'); 
    end
end

% 7. Colorbar settings
c = colorbar('Location', 'southoutside'); 
c.LineWidth = 0.5;
c.Ticks = [-30, 0, 30]; 
c.TickLabels = {'-30%', '0%', '30%'}; 

% Colorbar position
cpos = c.Position;
cpos(1) = 0.10;      % Left margin
cpos(2) = 0.12;      % Bottom position (including 2 mm upward shift)
cpos(3) = 0.80;      % Width
cpos(4) = 0.04;      % Thickness
c.Position = cpos;

hold off;

% 8. Force font settings
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');



%% OCD Step 5: Compute condition effects (final revised version: smaller ROI dots + 50x50 mm + colorbar moved upward)
%% ================================================================
clc; clear; close all;

%% 1) Load data
load('H:\MT_all(afterqujizhi)\cueoriandbu\correct\TOP\OCD\congruent\tICPC_all.mat', 'tICPC','parTF');  
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); 

fprintf('tICPC size = %s (sub × cond × chan × freq × time)\n', mat2str(size(tICPC)));

% Axes
labels    = freq.label(:);          
freq_axis = parTF.freq4tf(:)';      
time_axis = parTF.toi4tf(:)';       

% Analysis time window
ANALYSIS_WIN = [0.4 0.8];

%% 2) EGI → FieldTrip planar rotation & layout
elec = freq.elec;
x0 = elec.chanpos(:,1); y0 = elec.chanpos(:,2);
elec.chanpos(:,1) = -y0; elec.chanpos(:,2) =  x0;
if isfield(elec,'elecpos') && ~isempty(elec.elecpos)
    x0 = elec.elecpos(:,1); y0 = elec.elecpos(:,2);
    elec.elecpos(:,1) = -y0; elec.elecpos(:,2) =  x0;
end
cfgLay      = []; cfgLay.elec = elec;
lay         = ft_prepare_layout(cfgLay);

% Align layout
[~, ia, ib] = intersect(lay.label, labels, 'stable');
lay.label   = lay.label(ia);
lay.pos     = lay.pos(ia,:);
lay.width   = lay.width(ia);
lay.height  = lay.height(ia);
labels      = labels(ib);
tICPC       = tICPC(:,:,ib,:,:); 

%% 3) ROI channels
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2  3 121 122];
lChan    = [36 42 41 47];
rChan    = [92 97 101 102];

%% 4) Data structure & computation
ICPC            = struct();
ICPC.label      = labels;
ICPC.freq       = freq_axis;
ICPC.time       = time_axis;
ICPC.dimord     = 'chan_freq_time';

% Condition definition
CON   = [1 2 6];
INCON = [3 4 5];
LEFT  = [1 3];
RIGHT = [2 4];

% Compute (mask motor)
dat = tICPC;                              
dat(:,LEFT,  rChan,:,:) = NaN;            
dat(:,RIGHT, lChan,:,:) = NaN;            
ICPC.powspctrm = squeeze( nanmean( nanmean( dat(:,INCON,:,:,:) - dat(:,CON,:,:,:), 2), 1) );

%% ================================================================
% Plotting
% ================================================================

% 1. Set canvas (50 mm x 50 mm)
figW_cm = 5.0;
figH_cm = 5.0;

figure('name','ICPC_Congruency', ...
    'Units', 'centimeters', ...
    'Position', [10 10 figW_cm figH_cm], ...
    'Color', 'w', ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], ...
    'PaperSize', [figW_cm figH_cm]);

% 2. Topoplot configuration
cfg = [];
cfg.layout       = lay;
cfg.parameter    = 'powspctrm';
cfg.comment      = 'no';
cfg.style        = 'straight';
cfg.xlim         = ANALYSIS_WIN;
cfg.ylim         = [freq_axis(1) freq_axis(end)];
cfg.zlim         = [-30 30];       
cfg.colormap     = jet(64);

% Tiny black background dots (already set to minimum = 1)
cfg.marker       = 'on';
cfg.markersymbol = '.';
cfg.markercolor  = [0 0 0];
cfg.markersize   = 1;             
cfg.highlight    = 'off';         

% 3. Draw topography
ft_topoplotTFR(cfg, ICPC);

% 4. Adjust axis position
ax = gca;
set(ax, 'Position', [0.02, 0.20, 0.96, 0.80]); 

hold on;

% 5. Manually draw ROI scatter points
% ★★★ Change: reduce size from 15 to 8 ★★★
sz = 8;  
lw = 0.5; 

% Lateral PFC (dark red)
scatter(lay.pos(lpfcChan,1), lay.pos(lpfcChan,2), sz, [0.8 0.4 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);
scatter(lay.pos(rpfcChan,1), lay.pos(rpfcChan,2), sz, [0.8 0.4 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% Motor (dark blue)
scatter(lay.pos(lChan,1), lay.pos(lChan,2), sz, [0.4 0.7 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);
scatter(lay.pos(rChan,1), lay.pos(rChan,2), sz, [0.4 0.7 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% Midfrontal (white)
scatter(lay.pos(MFchan,1), lay.pos(MFchan,2), sz, [1 1 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth', lw);

% 6. Refine outline lines
hLines = findobj(ax, 'Type','line');
for i = 1:numel(hLines)
    mk = get(hLines(i),'Marker');
    if (ischar(mk) && strcmp(mk,'none')) || isempty(mk)
        set(hLines(i), 'LineWidth', 0.5, 'Color','k'); 
    end
end

% 7. Colorbar settings
c = colorbar('Location', 'southoutside'); 
c.LineWidth = 0.5;
c.Ticks = [-30, 0, 30]; 
c.TickLabels = {'-30%', '0%', '30%'}; 

% Colorbar position
cpos = c.Position;
cpos(1) = 0.10;      % Left margin
cpos(2) = 0.12;      % Bottom position (including 2 mm upward shift)
cpos(3) = 0.80;      % Width
cpos(4) = 0.04;      % Thickness
c.Position = cpos;

hold off;

% 8. Force font settings
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'normal');

