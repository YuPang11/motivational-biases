%%
%% Step 1: Merge all subject results into one large file
%% Final dimension: ICPCall_all(sub, condi, seed, chan, freq, time)
%% Merge all subject ISPC results
clc; clear; close all;

% === 1️⃣ Set path ===
baseDir = 'F:\PIT\correct_all\HC';
cd(baseDir);

fileList = dir('*_ISPC.mat');
nSub = numel(fileList);
fprintf('Detected %d subject files in total.\n', nSub);

% === 2️⃣ Initialize the large matrix ===
load(fileList(1).name, 'ICPCall', 'parTF');  % Get size info
nChan = parTF.nChan;
nFreq = numel(parTF.freq4tf);
nTime = numel(parTF.toi4tf);

ICPCall_all = NaN(nSub, 6, numel(parTF.MFchan), nChan, nFreq, nTime);
ctSub = 1;

% === 3️⃣ Loop through each subject ===
for i = 1:nSub
    fname = fileList(i).name;
    fprintf('Loading subject %d/%d: %s\n', i, nSub, fname);

    load(fname, 'ICPCall', 'parTF');

    % Skip this condition if it is empty (insufficient trial count)
    for iCondi = 1:6
        if isempty(ICPCall{iCondi})
            continue;
        end

disp(fname)
disp(size(ICPCall{iCondi}))

        ICPCall_all(ctSub, iCondi, :, :, :, :) = ICPCall{iCondi};
    end

    ctSub = ctSub + 1;
end

% === 4️⃣ Save in -v7.3 format ===
save(fullfile(baseDir, 'ICPCall_all.mat'), 'ICPCall_all', 'parTF', '-v7.3');
fprintf('\nSaved as ICPCall_all.mat (v7.3 format)\n');






clc; clear; close all;

% === 1️⃣ Set path ===
baseDir = 'F:\PIT\correct_all\OCD';
cd(baseDir);

fileList = dir('*_ISPC.mat');
nSub = numel(fileList);
fprintf('Detected %d subject files in total.\n', nSub);

% === 2️⃣ Initialize the large matrix ===
load(fileList(1).name, 'ICPCall', 'parTF');  % Get size info
nChan = parTF.nChan;
nFreq = numel(parTF.freq4tf);
nTime = numel(parTF.toi4tf);

ICPCall_all = NaN(nSub, 6, numel(parTF.MFchan), nChan, nFreq, nTime);
ctSub = 1;

% === 3️⃣ Loop through each subject ===
for i = 1:nSub
    fname = fileList(i).name;
    fprintf('Loading subject %d/%d: %s\n', i, nSub, fname);

    load(fname, 'ICPCall', 'parTF');

    % Skip this condition if it is empty (insufficient trial count)
    for iCondi = 1:6
        if isempty(ICPCall{iCondi})
            continue;
        end
        ICPCall_all(ctSub, iCondi, :, :, :, :) = ICPCall{iCondi};
    end

    ctSub = ctSub + 1;
end

% === 4️⃣ Save in -v7.3 format ===
save(fullfile(baseDir, 'ICPCall_all.mat'), 'ICPCall_all', 'parTF', '-v7.3');
fprintf('\nSaved as ICPCall_all.mat (v7.3 format)\n');





%%
%% Step 2: Baseline Correction
clc; clear; close all;

% === 1️⃣ Load the merged file ===
load('F:\PIT\correct_all\HC\ICPCall_all.mat');  % Contains ICPCall_all and parTF
disp('Loaded ICPCall_all.mat');

% === 2️⃣ Extract time indices ===
baseTime = dsearchn(parTF.toi4tf', parTF.baselinetime');  % Find baseline interval indices
baseIdx  = baseTime(1):baseTime(2);
fprintf('Baseline index range: %d ~ %d (%.3f–%.3f s)\n', baseIdx(1), baseIdx(end), ...
    parTF.toi4tf(baseIdx(1)), parTF.toi4tf(baseIdx(end)));

% === 3️⃣ Compute condition-averaged baseline ===
baseline_timeavg = nanmean(ICPCall_all(:,:,:,:,:,baseIdx), 6);% (a) Average across the time dimension (dim 6)

baseline_condavg = nanmean(baseline_timeavg, 2);% (b) Average across the condition dimension (dim 2), giving a condition-averaged baseline

nTime = numel(parTF.toi4tf);
baseline_expand = repmat(baseline_condavg, [1 6 1 1 1 nTime]);% (c) Expand back to 6 conditions × all time points

% === 4️⃣ Compute percentage change (% change) ===
ICPCall_uncorrected = ICPCall_all;  % Backup
ICPCall_all = (ICPCall_all - baseline_expand) ./ baseline_expand * 100;

% === 5️⃣ Save baseline-corrected result ===
save_path = 'F:\PIT\correct_all\HC\ICPCall_all_baseline.mat';
save(save_path, 'ICPCall_all', 'parTF', '-v7.3');
fprintf('Baseline correction completed. Result saved to:\n%s\n', save_path);



clc; clear; close all;

% === 1️⃣ Load the merged file ===
load('F:\PIT\correct_all\OCD\ICPCall_all.mat');  % Contains ICPCall_all and parTF
disp('Loaded ICPCall_all.mat');

% === 2️⃣ Extract time indices ===
baseTime = dsearchn(parTF.toi4tf', parTF.baselinetime');  % Find baseline interval indices
baseIdx  = baseTime(1):baseTime(2);
fprintf('Baseline index range: %d ~ %d (%.3f–%.3f s)\n', baseIdx(1), baseIdx(end), ...
    parTF.toi4tf(baseIdx(1)), parTF.toi4tf(baseIdx(end)));

% === 3️⃣ Compute condition-averaged baseline ===
baseline_timeavg = nanmean(ICPCall_all(:,:,:,:,:,baseIdx), 6);% (a) Average across the time dimension (dim 6)

baseline_condavg = nanmean(baseline_timeavg, 2);% (b) Average across the condition dimension (dim 2), giving a condition-averaged baseline

nTime = numel(parTF.toi4tf);
baseline_expand = repmat(baseline_condavg, [1 6 1 1 1 nTime]);% (c) Expand back to 6 conditions × all time points

% === 4️⃣ Compute percentage change (% change) ===
ICPCall_uncorrected = ICPCall_all;  % Backup
ICPCall_all = (ICPCall_all - baseline_expand) ./ baseline_expand * 100;

% === 5️⃣ Save baseline-corrected result ===
save_path = 'F:\PIT\correct_all\OCD\ICPCall_all_baseline.mat';
save(save_path, 'ICPCall_all', 'parTF', '-v7.3');
fprintf('Baseline correction completed. Result saved to:\n%s\n', save_path);




%% 
%% Step 3: t-weighted Midfrontal ISPC Aggregation
clc; clear; close all;

% === 1️⃣ Load the baseline-corrected file ===
load('F:\PIT\correct_all\HC\ICPCall_all_baseline.mat');
disp('Loaded baseline-corrected ICPCall_all');

% === 2️⃣ Define parameters ===
MFchans = [6 127];  % Midfrontal channels (FCz, Cz)
tvals   = [2.9194, 3.4082];  % Corresponding t values%tvals   = [3.2878, 3.8499];  % Corresponding t values
%tvals   = [5.3689, 5.0774];  % Corresponding t values
weights = tvals ./ sum(tvals);  % Normalized weights
fprintf('Weights = [%.4f, %.4f]\n', weights(1), weights(2));

% === 3️⃣ Compute t-weighted ISPC ===
% ICPCall_all structure: (sub, condi, seed, chan, freq, time)
% Apply weighting on the 3rd dimension (seed channel)
tICPC = squeeze( ...
    weights(1) .* ICPCall_all(:,:,1,:,:,:) + ...
    weights(2) .* ICPCall_all(:,:,2,:,:,:) ...
);

% Check output dimensions
disp('tICPC dimensions: (sub × condi × chan × freq × time)');
disp(size(tICPC));

% === 4️⃣ Save result ===
save_path = 'F:\PIT\correct_all\HC\tICPC_all.mat';
save(save_path, 'tICPC', 'weights', 'MFchans', 'parTF', '-v7.3');

fprintf('Saved weighted result file:\n%s\n', save_path);
fprintf('Included variables: tICPC (weighted ISPC), weights, MFchans, parTF\n');






clc; clear; close all;

% === 1️⃣ Load the baseline-corrected file ===
load('F:\PIT\correct_all\OCD\ICPCall_all_baseline.mat');
disp('Loaded baseline-corrected ICPCall_all');

% === 2️⃣ Define parameters ===
MFchans = [6 127];  % Midfrontal channels (FCz, Cz)
tvals   = [3.4548, 4.4553];  % Corresponding t values
%tvals   = [3.1478, 4.4980];  % Corresponding t values
weights = tvals ./ sum(tvals);  % Normalized weights
fprintf('Weights = [%.4f, %.4f]\n', weights(1), weights(2));

% === 3️⃣ Compute t-weighted ISPC ===
% ICPCall_all structure: (sub, condi, seed, chan, freq, time)
% Apply weighting on the 3rd dimension (seed channel)
tICPC = squeeze( ...
    weights(1) .* ICPCall_all(:,:,1,:,:,:) + ...
    weights(2) .* ICPCall_all(:,:,2,:,:,:) ...
);

% Check output dimensions
disp('tICPC dimensions: (sub × condi × chan × freq × time)');
disp(size(tICPC));

% === 4️⃣ Save result ===
save_path = 'F:\PIT\correct_all\OCD\tICPC_all.mat';
save(save_path, 'tICPC', 'weights', 'MFchans', 'parTF', '-v7.3');

fprintf('Saved weighted result file:\n%s\n', save_path);
fprintf('Included variables: tICPC (weighted ISPC), weights, MFchans, parTF\n');



%% 
%% Step 4: Visualization Check (Topoplot)

clc; clear; close all;

%% 1) Load data
load('H:\MT_all(afterqujizhi)\cueoriandbu\correct\change1\tICPC_all.mat', 'tICPC', 'parTF'); % tICPC, parTF
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); % Use freq.elec/freq.label

fprintf('tICPC size = %s (sub × cond × chan × freq × time)\n', mat2str(size(tICPC)));

%% 2) Group average (subject × condition) → [chan × freq × time]
ICPC_avg = squeeze(nanmean(nanmean(tICPC, 1), 2));   % 127 × 7 × 131

% Axes and labels
freq_axis = parTF.freq4tf(:)';      % 7
time_axis = parTF.toi4tf(:)';       % 131
labels    = freq.label(:);          % 127

%% 3) Theta-band index (4–8 Hz)
theta_band = [4 8];
iFreq = dsearchn(freq_axis', theta_band');
iFreq = iFreq(1):iFreq(2);          % Cover 4–8 Hz

%% 4) Define 10 time windows (same as the original paper)
win = [ ...
   -0.25  -0.05; ...
    0.00   0.20; ...
    0.20   0.40; ...
    0.40   0.60; ...
    0.60   0.80; ...
    0.80   1.00; ...
    1.30   1.50; ...
    2.00   2.20; ...
    2.20   2.40; ...
    2.40   2.60];
titles = {'-250–-50 ms','0–200 ms','200–400 ms','400–600 ms','600–800 ms', ...
          '800–1000 ms','1300–1500 ms','2000–2200 ms','2200–2400 ms','2400–2600 ms'};

%% 5) Prepare electrode coordinates (EGI -> FieldTrip rotation fix)
elec = freq.elec;
x0 = elec.chanpos(:,1); y0 = elec.chanpos(:,2);
elec.chanpos(:,1) = -y0;            % x = -y
elec.chanpos(:,2) =  x0;            % y =  x
if isfield(elec,'elecpos') && ~isempty(elec.elecpos)
    x0 = elec.elecpos(:,1); y0 = elec.elecpos(:,2);
    elec.elecpos(:,1) = -y0;
    elec.elecpos(:,2) =  x0;
end

% Generate layout once
cfgLay      = [];
cfgLay.elec = elec;
lay         = ft_prepare_layout(cfgLay);

% Align layout order with data label order
[~, ia, ib] = intersect(lay.label, labels, 'stable');
lay.label   = lay.label(ia);
lay.pos     = lay.pos(ia,:);
lay.width   = lay.width(ia);
lay.height  = lay.height(ia);
labels      = labels(ib);
ICPC_avg    = ICPC_avg(ib,:,:);      % 127 × 7 × 131

%% 6) Output directory
outDir = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\change1\ICPC_topos_theta';
if ~exist(outDir,'dir'); mkdir(outDir); end

%% 7) For each time window: average over frequency+time → single-channel vector → build ER → plot/save/close
basecfg             = [];
basecfg.layout      = lay;
basecfg.parameter   = 'avg';
basecfg.zlim        = [-50 50];      % Same as the check figure
basecfg.comment     = 'no';
basecfg.style       = 'straight';
basecfg.marker      = 'on';
colormap('jet');

for k = 1:size(win,1)
    % Time indices
    iTime = dsearchn(time_axis', win(k,:)');
    iTime = iTime(1):iTime(2);

    % Average within 4–8 Hz and current time window → 127×1
    topo_k = squeeze(mean(mean(ICPC_avg(:, iFreq, iTime), 3), 2));

    % Build a minimal ER structure (single "time point")
    ER = [];
    ER.label  = labels;
    ER.time   = 1;                   % Single point
    ER.avg    = topo_k;              % 127×1
    ER.dimord = 'chan_time';

    % Create new figure, draw, save, close
    h = figure('Visible','off','Position',[100 100 700 600]); %#ok<NASGU>
    cfg = basecfg;
    cfg.xlim = [1 1];                % Plot the only "time point"
    ft_topoplotER(cfg, ER);
    title(sprintf('All responses (theta 4–8 Hz)\n%s', titles{k}), 'FontSize', 12, 'FontWeight', 'bold');
    hcb = colorbar; ylabel(hcb,'ICPC (% change)');
    print(fullfile(outDir, sprintf('ICPC_theta_topo_%02d_%s.png', k, titles{k})), '-dpng', '-r300');
    close gcf;
end

fprintf('Done: 10 topo figures were exported to:\n%s\n', outDir);





%% Step 5: Channel Selection

clc; clear; close all;

% 1) Load data
load('H:\MT_all(afterqujizhi)\cueoriandbu\correct\TOP\OCD\pfc\tICPC_all.mat', 'tICPC', 'parTF'); % tICPC, parTF
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); % Use freq.elec/freq.label

%% Channel selection — Part A (condition-averaged)
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2  3 121 122];

% ---- Build FieldTrip freq structure (chan × freq × time) ----
ICPC              = struct();
ICPC.label        = freq.label;            % 127 channel labels
ICPC.freq         = parTF.freq4tf;         % 7 frequency points in 4–8 Hz
ICPC.time         = parTF.toi4tf;          % 131 time points
ICPC.powspctrm = squeeze(nanmean(nanmean(tICPC,2))); % chan×freq×time
ICPC.dimord       = 'chan_freq_time';

% ---- Coordinate transform (EGI→FT) ----
elec              = freq.elec;
elec.chanpos(:,1) = -freq.elec.chanpos(:,2);  % x = -y
elec.chanpos(:,2) =  freq.elec.chanpos(:,1);  % y =  x
ICPC.elec         = elec;

% ---- 1) Draw topoplot first (do not use highlight) ----
figure('name','ICPC: condition averaged','Position',[50 100 300 250]);
colormap jet
cfg            = [];
cfg.elec       = ICPC.elec;
cfg.comment    = 'no';
cfg.style      = 'straight';
cfg.marker     = 'on';
cfg.zlim       = [-30 30];
cfg.ylim       = [ICPC.freq(1) ICPC.freq(end)];  % 4–8 Hz
cfg.xlim       = [0.4 0.65];                    % Time window
ft_topoplotTFR(cfg, ICPC);
colorbar; hold on

% ---- 2) Get 2D channel positions on the plot ----
laycfg      = []; 
laycfg.elec = ICPC.elec;               % Use the same electrode coordinates as above
lay         = ft_prepare_layout(laycfg);

% Convert channel indices to logical indices in lay.label
lab     = ICPC.label;
selPFC  = ismember(lay.label, lab([lpfcChan rpfcChan]));
selMF   = ismember(lay.label, lab(MFchan));

% ---- 3) Overlay filled markers (purple PFC + white MF) ----
sz = 90;  % Marker size
% PFC (purple filled, black edge)
scatter(lay.pos(selPFC,1), lay.pos(selPFC,2), sz, [0.55 0 0.9], ...
        'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
% Midfrontal seeds (white filled, black edge)
scatter(lay.pos(selMF,1),  lay.pos(selMF,2),  sz, [1 1 1], ...
        'filled', 'MarkerEdgeColor','k','LineWidth',1.2);

hold off

%% ===== Compute and export the mean ICPC value of each electrode =====

% Check that ICPC exists first
whos ICPC

% 1️⃣ Average across all frequency and time dimensions (keep channels only)
chanMean = squeeze(mean(mean(ICPC.powspctrm, 3, 'omitnan'), 2, 'omitnan'));  % [127×1]
chanMean = chanMean(:);  % Convert to column vector

% 2️⃣ Build table
T = table(ICPC.label(:), chanMean, 'VariableNames', {'Channel', 'Mean_ICPC'});

% 3️⃣ Export to Excel
outXlsx = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\OCD\ICPC_channel_mean.xlsx';
writetable(T, outXlsx);

fprintf('Exported the overall mean of %d electrodes to:\n%s\n', height(T), outXlsx);








%% Channel selection — Part B (condition-averaged)
MFchan = [6 127];
lChan  = [36 42 41 47];
rChan  = [92 97 101 102];

% ---- Build FieldTrip freq structure (chan × freq × time) ----
ICPC              = struct();
ICPC.label        = freq.label;            % 127 channel labels
ICPC.freq         = parTF.freq4tf;         % 7 frequency points in 4–8 Hz
ICPC.time         = parTF.toi4tf;          % 131 time points
ICPC.powspctrm    = squeeze( nanmean( nanmean( ...
                        tICPC(:,[1 3],:,:,:) - tICPC(:,[2 4],:,:,:), 2), 1) );
ICPC.dimord       = 'chan_freq_time';

% ---- Coordinate transform (EGI→FT) ----
elec              = freq.elec;
elec.chanpos(:,1) = -freq.elec.chanpos(:,2);  % x = -y
elec.chanpos(:,2) =  freq.elec.chanpos(:,1);  % y =  x
ICPC.elec         = elec;

% ---- Draw the base topoplot ----
figure('name','ICPC: condition averaged','Position',[50 100 250 220]);
colormap jet

cfg = [];
cfg.elec      = ICPC.elec;
cfg.comment   = 'no';
cfg.style     = 'straight';
cfg.marker    = 'on';
cfg.zlim      = [-30 30];
cfg.ylim      = [ICPC.freq(1) ICPC.freq(end)];  % theta 4–8 Hz
cfg.xlim      = [0.424 0.8];                    % Time window
ft_topoplotTFR(cfg, ICPC);
colorbar; hold on;

% ---- Get layout coordinates ----
laycfg      = [];
laycfg.elec = ICPC.elec;
lay         = ft_prepare_layout(laycfg);

% ---- Get electrode labels ----
lab = ICPC.label;
selChan = ismember(lay.label, lab([lChan rChan]));
selMF   = ismember(lay.label, lab(MFchan));

% ---- Overlay color markers ----
sz = 90;  % Marker size
% Blue: left and right parietal target channels
scatter(lay.pos(selChan,1), lay.pos(selChan,2), sz, [0 0 1], ...
        'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
% White: midfrontal seed channels
scatter(lay.pos(selMF,1), lay.pos(selMF,2), sz, [1 1 1], ...
        'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
hold off;





%% Step 5: Compute Condition Effects
%% ================================================================
%  Condition effects (topographies + motor time courses)
%  Based on Swart et al. (2016) logic + current condition order/channels/time window/custom coordinates
%  Requirement: tICPC: [sub × cond × chan × freq × time], already t-weighted & baseline-corrected
%               parTF.freq4tf ≈ 4–8 Hz (7 points), parTF.toi4tf = -0.25:0.025:3
%               freq.elec/freq.label come from the current data (for real electrode coordinates)
% ================================================================

clc; clear; close all;

%% 1) Load data
load('F:\PIT\correct_all\HC\tICPC_all.mat', 'tICPC','parTF');  % tICPC, parTF
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); % Use freq.elec/label only

fprintf('tICPC size = %s (sub × cond × chan × freq × time)\n', mat2str(size(tICPC)));

% Axes
labels    = freq.label(:);          % 127
freq_axis = parTF.freq4tf(:)';      % 7
time_axis = parTF.toi4tf(:)';       % 131

% Analysis time window (change to [0.45 0.65] if needed to match the paper)
ANALYSIS_WIN = [0.45 0.7];

%% 2) EGI → FieldTrip planar rotation (x=-y, y=x), and prepare layout once
elec = freq.elec;
x0 = elec.chanpos(:,1); y0 = elec.chanpos(:,2);
elec.chanpos(:,1) = -y0; elec.chanpos(:,2) =  x0;
if isfield(elec,'elecpos') && ~isempty(elec.elecpos)
    x0 = elec.elecpos(:,1); y0 = elec.elecpos(:,2);
    elec.elecpos(:,1) = -y0; elec.elecpos(:,2) =  x0;
end
cfgLay      = []; cfgLay.elec = elec;
lay         = ft_prepare_layout(cfgLay);

% Safely align layout order with the data label order
[~, ia, ib] = intersect(lay.label, labels, 'stable');
lay.label   = lay.label(ia);
lay.pos     = lay.pos(ia,:);
lay.width   = lay.width(ia);
lay.height  = lay.height(ia);
labels      = labels(ib);
tICPC       = tICPC(:,:,ib,:,:);  % Reorder channels to match layout

%% 3) ROI channels
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2  3 121 122];
lChan    = [36 42 41 47];
rChan    = [92 97 101 102];

%% 4) One ICPC skeleton (only change powspctrm for each contrast)
ICPC            = struct();
ICPC.label      = labels;
ICPC.freq       = freq_axis;
ICPC.time       = time_axis;
ICPC.dimord     = 'chan_freq_time';

% Common plotting cfg base
basecfg = [];
basecfg.layout   = lay;                  % More stable with layout
basecfg.parameter= 'powspctrm';
basecfg.comment  = 'no';
basecfg.style    = 'straight';
basecfg.marker   = 'on';
basecfg.ylim     = [freq_axis(1) freq_axis(end)]; % theta band
basecfg.xlim     = ANALYSIS_WIN;

%% 5) Condition sets
% 1=GLW, 2=GRW, 3=GLA, 4=GRA, 5=NGW, 6=NGA
CON   = [1 2 6];
INCON = [3 4 5];
LEFT  = [1 3];
RIGHT = [2 4];

%% ================================================================
% Compute condition effects
% ================================================================

%% (1) Congruency effect independent of response: (Avoid − Win)
figure('name','ICPC: congruency effect','Position',[50 100 200 200]); colormap jet
cfg = basecfg; cfg.zlim = [-20 20];
cfg.highlight = 'on'; cfg.highlightchannel = [lpfcChan rpfcChan MFchan];
cfg.highlightsymbol='o'; cfg.highlightsize=12; cfg.highlightcolor=[0 0 0];
ICPC.powspctrm = squeeze( nanmean( nanmean( tICPC(:,INCON,:,:,:) - tICPC(:,CON,:,:,:), 2), 1) );
ft_topoplotTFR(cfg, ICPC); colorbar;

%% (2) Congruency effect for NoGo responses: NGW − NGA = 5 − 6
figure('name','ICPC: win - avoid for nogo responses','Position',[50 100 200 200]); colormap jet
cfg = basecfg; cfg.zlim = [-20 20];
cfg.highlight = 'on'; cfg.highlightchannel = [lChan rChan MFchan];
cfg.highlightsymbol='o'; cfg.highlightsize=12; cfg.highlightcolor=[0 0 0];
ICPC.powspctrm = squeeze( nanmean( nanmean( tICPC(:,5,:,:,:) - tICPC(:,6,:,:,:), 2), 1) );
ft_topoplotTFR(cfg, ICPC); colorbar;

%% (3) Congruency effect for Go left/right separately: Avoid − Win
% Left: 3 − 1, Right: 4 − 2 (under the current condition order)
% === Plot setup ===
cfgBase = basecfg; 
cfgBase.zlim = [-20 20];
cfgBase.highlight = 'on';
cfgBase.highlightchannel = [lChan rChan MFchan];
cfgBase.highlightsymbol = 'o';
cfgBase.highlightsize = 12;
cfgBase.highlightcolor = [0 0 0];

% === Create figure (using tiledlayout) ===
hFig = figure('Position',[600,600,520,250]);
tiledlayout(hFig,1,2,'TileSpacing','compact');
colormap(hFig,'jet');

% === Left: avoid - win ===
nexttile;
cfg = cfgBase;
cfg.figure = 'gca';   % Tell FieldTrip to draw in the current tile
ICPC.powspctrm = squeeze( nanmean( nanmean( ...
                        tICPC(:,3,:,:,:) - tICPC(:,1,:,:,:), 2), 1) );
ft_topoplotTFR(cfg, ICPC);
title('Left: avoid - win');

% === Right: avoid - win ===
nexttile;
cfg = cfgBase;
cfg.figure = 'gca';
ICPC.powspctrm = squeeze( nanmean( nanmean( ...
                        tICPC(:,4,:,:,:) - tICPC(:,2,:,:,:), 2), 1) );
ft_topoplotTFR(cfg, ICPC);
title('Right: avoid - win');

colorbar;

%% (4) Congruency effect for Left - Right Go responses:
% Original paper: [1 4] - [2 3], keep unchanged under the current condition order
cfgBase = basecfg;cfgBase.zlim = [-20 20];cfgBase.highlight = 'on';cfgBase.highlightchannel = [lChan rChan MFchan];
cfg.highlightsymbol='o'; cfg.highlightsize=12; cfg.highlightcolor=[0 0 0];
ICPC.powspctrm = squeeze(nanmean(nanmean(tICPC(:,[1 4],:,:,:)-tICPC(:,[3 2],:,:,:),2))); 
ft_topoplotTFR(cfg, ICPC); colorbar;
title('Go(avoid - win): Left - Right');

%% (4) Congruency effect collapsed over responses, but WITHOUT contralateral motor (Figure 5 in the paper)
figure('name','ICPC: congruency effect (masked motor)', 'Position', [50 100 400 400]);
colormap jet
cfg = basecfg;
cfg.zlim = [-20 20];
cfg.highlight = 'off';  % Turn off default highlight and use custom markers

% ===== Data preparation =====
dat = tICPC;                              % sub × cond × chan × freq × time
dat(:,LEFT,  rChan,:,:) = NaN;            % Left responses → mask right motor area
dat(:,RIGHT, lChan,:,:) = NaN;            % Right responses → mask left motor area
ICPC.powspctrm = squeeze( nanmean( nanmean( dat(:,INCON,:,:,:) - dat(:,CON,:,:,:), 2), 1) );

% ===== Draw topography =====
ft_topoplotTFR(cfg, ICPC);
colorbar; hold on;

% ===== Get electrode layout coordinates =====
if ~exist('lay','var')
    lay = ft_prepare_layout([], ICPC);
end

% ===== Set marker size and color =====
sz = 100;  % Dot size

% ==== Dark red: left/right lateral PFC channels ====
scatter(lay.pos(lpfcChan,1), lay.pos(lpfcChan,2), sz, [0.5 0 0], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth',1.2);
scatter(lay.pos(rpfcChan,1), lay.pos(rpfcChan,2), sz, [0.5 0 0], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth',1.2);

% ==== Dark blue: left/right motor channels ====
scatter(lay.pos(lChan,1), lay.pos(lChan,2), sz, [0 0 0.5], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth',1.2);
scatter(lay.pos(rChan,1), lay.pos(rChan,2), sz, [0 0 0.5], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth',1.2);

% ==== White: midfrontal seed channels ====
scatter(lay.pos(MFchan,1), lay.pos(MFchan,2), sz, [1 1 1], ...
        'filled', 'MarkerEdgeColor','k', 'LineWidth',1.5);

hold off;
title('Incongruent – Congruent');

%% ===== Extract the Incongruent–Congruent mean value of each electrode within timeROI =====

% --- Confirm time-window indices ---
timeROI = [0.4 0.65];
iTime   = dsearchn(parTF.toi4tf', timeROI'); 
iTime   = iTime(1):iTime(2);

% --- Average over the frequency dimension (dim 2) and within timeROI ---
chanDiff = squeeze(mean(mean(ICPC.powspctrm(:, :, iTime), 3, 'omitnan'), 2, 'omitnan'));  % [127×1]

% --- Build table and export ---
T = table(ICPC.label(:), chanDiff(:), ...
    'VariableNames', {'Channel', 'ICPC_Incongruent_minus_Congruent'});

outXlsx = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\TOP\OCD\congruent\ICPC_Congruent.xlsx';
writetable(T, outXlsx);

fprintf('Exported the Incongruent–Congruent mean values of %d electrodes within %.3f–%.3f s to:\n%s\n', ...
        height(T), timeROI(1), timeROI(2), outXlsx);

%% (6) Left / Right / NoGo (collapsed across valence)
% First prepare the shared plotting base (inherits from basecfg)
cfgBase = basecfg;                 
cfgBase.zlim = [-30 30];           
cfgBase.highlight = 'on';
cfgBase.highlightchannel = [lChan rChan lpfcChan rpfcChan MFchan];
cfgBase.highlightsymbol = 'o';
cfgBase.highlightsize   = 12;
cfgBase.highlightcolor  = [0 0 0];

% Canvas and grid
hFig = figure('Name','ICPC: Left / Right / NoGo (valence collapsed)', ...
              'Position',[600,600,1000,280]);   
tiledlayout(hFig,1,3,'TileSpacing','compact','Padding','compact');
colormap(hFig,'jet');

% --- Left responses: mean over cond=[LEFT], then average over subjects → chan×freq×time
nexttile;
cfg = cfgBase; 
cfg.figure = 'gca';                                  
ICPC.powspctrm = squeeze(nanmean(nanmean(tICPC(:,LEFT,:,:,:),  2), 1));
ft_topoplotTFR(cfg, ICPC);
title('Left responses');

% --- Right responses: cond=[RIGHT]
nexttile;
cfg = cfgBase; 
cfg.figure = 'gca';
ICPC.powspctrm = squeeze(nanmean(nanmean(tICPC(:,RIGHT,:,:,:), 2), 1));
ft_topoplotTFR(cfg, ICPC);
title('Right responses');

% --- NoGo responses: cond=[5 6]
nexttile;
cfg = cfgBase; 
cfg.figure = 'gca';
ICPC.powspctrm = squeeze(nanmean(nanmean(tICPC(:,[5 6],:,:,:), 2), 1));
ft_topoplotTFR(cfg, ICPC);
title('NoGo responses');

% Unified colorbar (place on the far right)
cb = colorbar;
cb.Layout.Tile = 'east';           
ylabel(cb,'ICPC (% change)');


%% (7) timeseries motor phase synchrony (same logic as Fig. 6B in the paper)
rightmotor = squeeze( mean( mean(tICPC(:,:,rChan,:,:), 3), 4) ); % sub×cond×time
leftmotor  = squeeze( mean( mean(tICPC(:,:,lChan,:,:), 3), 4) ); % sub×cond×time

contra_win   = ( squeeze(nanmean(rightmotor(:,1,:),1)) + squeeze(nanmean(leftmotor(:,2,:),1)) )/2; % GLW + GRW
contra_avoid = ( squeeze(nanmean(rightmotor(:,3,:),1)) + squeeze(nanmean(leftmotor(:,4,:),1)) )/2; % GLA + GRA
ipsi_win     = ( squeeze(nanmean(rightmotor(:,2,:),1)) + squeeze(nanmean(leftmotor(:,1,:),1)) )/2; % GRW + GLW
ipsi_avoid   = ( squeeze(nanmean(rightmotor(:,4,:),1)) + squeeze(nanmean(leftmotor(:,3,:),1)) )/2; % GRA + GLA
nogo_win     = ( squeeze(nanmean(rightmotor(:,5,:),1)) + squeeze(nanmean(leftmotor(:,5,:),1)) )/2; % NGW
nogo_avoid   = ( squeeze(nanmean(rightmotor(:,6,:),1)) + squeeze(nanmean(leftmotor(:,6,:),1)) )/2; % NGA

figure('name','midfrontal TF power Valence x Required action','Position',[50 100 600 600]); hold on
plot(time_axis, contra_win,   'Color',[0 .8 0], 'LineWidth',2)   % contra-win
plot(time_axis, contra_avoid, 'Color',[.8 0 0], 'LineWidth',2)   % contra-avoid
plot(time_axis, ipsi_win,     'Color',[0 .5 0], 'LineWidth',2)   % ipsi-win
plot(time_axis, ipsi_avoid,   'Color',[.5 0 0], 'LineWidth',2)   % ipsi-avoid
plot(time_axis, nogo_win,     '--','Color',[0 .8 0],'LineWidth',2)
plot(time_axis, nogo_avoid,   '--','Color',[.8 0 0],'LineWidth',2)
set(gca,'ylim',[-10 50],'xlim',[-.25 1.3],'fontsize',11); box off
legend({'contra-win','contra-avoid','ipsi-win','ipsi-avoid','nogo-win','nogo-avoid'},'Location','northwest'); legend boxoff
plot([0 0], get(gca,'ylim'),'k')                        % cue
plot([ANALYSIS_WIN(1) ANALYSIS_WIN(1)],get(gca,'ylim'),':k')  % selected window
plot([ANALYSIS_WIN(2) ANALYSIS_WIN(2)],get(gca,'ylim'),':k')
ylabel('ICPC (% change)'); xlabel('Time (s)');
title('Motor ROIs: theta ICPC time courses (group mean)');

disp('Done: condition-effect topographies + motor time courses');




%% ===============================================================
%  Final pipeline: Export ICPC for SPSS + Replicate bar plots
%  tICPC dims: sub × cond × chan × freq × time
%  Conditions: 1=GLW, 2=GRW, 3=GLA, 4=GRA, 5=NGW, 6=NGA
% ===============================================================

clc; clear; close all;

%% 1) Load data
% load('F:\PIT\correct_all\HC\results\tICPC_all.mat', 'tICPC','parTF');   % tICPC, parTF
% load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); % Use freq.elec/label only
load('F:\PIT\correct_all\OCD\tICPC_all.mat', 'tICPC','parTF');   % tICPC, parTF
load('H:\MT_all(afterqujizhi)\final\correct_cue\cue_epoch\cue_baseline\ft_lap\HC02.mat', 'freq'); % Use freq.elec/label only


%% 2) Params, ROI, time window
dirs.home   = 'F:\PIT\correct_all\HC';
results_dir = fullfile(dirs.home, 'results');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
csvfilename = fullfile(results_dir, 'ICPC_contra_ipsiocd0.4.8plottest.csv');

% Time window → indices
%timeROI = [0.45 0.65];
timeROI = [0.4 0.8];%ocd
iTime   = dsearchn(parTF.toi4tf', timeROI'); 
iTime   = iTime(1):iTime(2);

% ROI channels
lpfcChan = [23 24 26 27]; 
rpfcChan = [2 3 121 122];
lChan    = [36 42 41 47]; 
rChan    = [92 97 101 102];

nSub     = size(tICPC,1);

%% 3) Compute ICPC: Parietal/Motor (midfrontal–motor)
% -----------------------------------------------------
% TF_par(sub, valence[win/avoid], response[contra/ipsi/nogo])
% contra = contralateral motor area; ipsi = ipsilateral motor area; nogo = bilateral average
% -----------------------------------------------------
TF_par = nan(nSub,2,3);

% === Contra ===
TF_par(:,1,1) = mean(mean(mean(nanmean(cat(6,tICPC(:,1,rChan,:,iTime),tICPC(:,2,lChan,:,iTime)),6),5),4),3);  % GLW(left hand→right) + GRW(right hand→left) contra-win
TF_par(:,2,1) = mean(mean(mean(nanmean(cat(6,tICPC(:,3,rChan,:,iTime),tICPC(:,4,lChan,:,iTime)),6),5),4),3);  % GLA(left hand→right) + GRA(right hand→left) contra-avoid

% === Ipsi ===
TF_par(:,1,2) = mean(mean(mean(nanmean(cat(6,tICPC(:,1,lChan,:,iTime),tICPC(:,2,rChan,:,iTime)),6),5),4),3);  % GLW(left hand→left) + GRW(right hand→right) ipsi-win
TF_par(:,2,2) = mean(mean(mean(nanmean(cat(6,tICPC(:,3,lChan,:,iTime),tICPC(:,4,rChan,:,iTime)),6),5),4),3);  % GLA(left hand→left) + GRA(right hand→right) ipsi-avoid

% === NoGo ===
TF_par(:,1,3) = mean(mean(mean(tICPC(:,5,[lChan rChan],:,iTime),5),4),3);  % NGW
TF_par(:,2,3) = mean(mean(mean(tICPC(:,6,[lChan rChan],:,iTime),5),4),3);  % NGA


%% 4) Compute ICPC: Lateral PFC (midfrontal–PFC)
% -----------------------------------------------------
% TF_pfc(sub, valence[win/avoid], response[contra/ipsi/nogo])
% contra = contralateral PFC; ipsi = ipsilateral PFC; nogo = bilateral average
% -----------------------------------------------------
TF_pfc = nan(nSub,2,3);

% === Contra ===
TF_pfc(:,1,1) = mean(mean(mean(nanmean(cat(6,tICPC(:,1,rpfcChan,:,iTime),tICPC(:,2,lpfcChan,:,iTime)),6),5),4),3);  % GLW(left hand→right) + GRW(right hand→left)
TF_pfc(:,2,1) = mean(mean(mean(nanmean(cat(6,tICPC(:,3,rpfcChan,:,iTime),tICPC(:,4,lpfcChan,:,iTime)),6),5),4),3);  % GLA(left hand→right) + GRA(right hand→left)

% === Ipsi ===
TF_pfc(:,1,2) = mean(mean(mean(nanmean(cat(6,tICPC(:,1,lpfcChan,:,iTime),tICPC(:,2,rpfcChan,:,iTime)),6),5),4),3);  % GLW(left hand→left) + GRW(right hand→right)
TF_pfc(:,2,2) = mean(mean(mean(nanmean(cat(6,tICPC(:,3,lpfcChan,:,iTime),tICPC(:,4,rpfcChan,:,iTime)),6),5),4),3);  % GLA(left hand→left) + GRA(right hand→right)

% === NoGo ===
TF_pfc(:,1,3) = mean(mean(mean(tICPC(:,5,[lpfcChan rpfcChan],:,iTime),5),4),3);  % NGW
TF_pfc(:,2,3) = mean(mean(mean(tICPC(:,6,[lpfcChan rpfcChan],:,iTime),5),4),3);  % NGA


%% 5) Export to CSV (SPSS/R friendly)
valString  = {'win';'avoid'};
respString = {'contra';'ipsi';'nogo'};
hdr = {}; data2save = [];

subjects = (1:nSub)'; 
if exist('sub2exclude','var') && ~isempty(sub2exclude)
    subjects(sub2exclude) = [];
    TF_par = TF_par(subjects,:,:);
    TF_pfc = TF_pfc(subjects,:,:);
end
nSub_eff = numel(subjects);

data2save(:,1) = subjects; 
hdr{1} = 'subject,';
varc = 1;

% --- PFC ---
for iVal = 1:2
    for iResp = 1:3
        varc = varc + 1; 
        data2save(:,varc) = TF_pfc(:,iVal,iResp);
        hdr{varc} = sprintf('icpc_pfc_%s_%s,',valString{iVal},respString{iResp});
    end
end

% --- Motor ---
for iVal = 1:2
    for iResp = 1:3
        varc = varc + 1; 
        data2save(:,varc) = TF_par(:,iVal,iResp);
        hdr{varc} = sprintf('icpc_par_%s_%s,',valString{iVal},respString{iResp});
    end
end

hdrAll = cat(2,hdr{:});
fid = fopen(csvfilename,'w'); fprintf(fid,'%s\n',hdrAll(1:end-1)); fclose(fid);
writematrix(data2save,csvfilename,'WriteMode','append');
fprintf('\nCSV saved: %s\n', csvfilename);



%% 6) Barplots.

% -------- A. Midfrontal–Lateral PFC --------
% Go = mean of contra & ipsi; NoGo = nogo
TF_pfc2plot = nan(nSub_eff,2,2);
TF_pfc2plot(:,:,1) = nanmean(TF_pfc(:,:,1:2),3);  % Go (avg contra & ipsi)
TF_pfc2plot(:,:,2) = TF_pfc(:,:,3);               % NoGo

figure('Position',[100 300 250 250]); 
h1 = bar(1:2, squeeze(nanmean(TF_pfc2plot(:,1,:))), 0.4, 'FaceColor',[0 0.6 0]); hold on; % Win
h2 = bar(1.45:2.45, squeeze(nanmean(TF_pfc2plot(:,2,:))), 0.4, 'FaceColor',[0.7 0 0]);     % Avoid
errorbar(1:2, squeeze(nanmean(TF_pfc2plot(:,1,:))), squeeze(nanstd(TF_pfc2plot(:,1,:)))./sqrt(nSub_eff), 'k','linestyle','none');
errorbar(1.45:2.45, squeeze(nanmean(TF_pfc2plot(:,2,:))), squeeze(nanstd(TF_pfc2plot(:,2,:)))./sqrt(nSub_eff), 'k','linestyle','none');
set(gca,'xtick',1.23:2.23,'xticklabel',{'Go','NoGo'},'xlim',[.7 2.7],'ylim',[0 25],'ytick',0:10:20);
ylabel('ICPC (% change)'); legend([h1 h2], {'Win','Avoid'}, 'Location','northeast'); title('Midfrontal–Lateral PFC'); box off;

% -------- B. Midfrontal–Motor (exclude contralateral) --------
figure('Position',[100 300 250 250]);
h1 = bar(1:2,      squeeze(nanmean(TF_par(:,1,2:3))), .4, 'facecolor',[0 0.6 0]); hold on;
h2 = bar(1.45:2.45,squeeze(nanmean(TF_par(:,2,2:3))), .4, 'facecolor',[0.7 0 0]);
errorbar(1:2,       squeeze(nanmean(TF_par(:,1,2:3))), squeeze(nanstd(TF_par(:,1,2:3)))./sqrt(nSub_eff), 'k','linestyle','none');
errorbar(1.45:2.45, squeeze(nanmean(TF_par(:,2,2:3))), squeeze(nanstd(TF_par(:,2,2:3)))./sqrt(nSub_eff), 'k','linestyle','none');
set(gca,'xtick',1.23:2.23,'xticklabel',{'Go','NoGo'},'xlim',[.7 2.7],'ylim',[0 25],'ytick',0:10:20);
ylabel('ICPC (% change)'); legend([h1 h2], {'Win','Avoid'}, 'Location','northeast'); title('Midfrontal–Motor'); box off;

% -------- C. Midfrontal–Motor (include contralateral) --------
figure('Position',[100 300 300 250]);
h1 = bar(1:3,       squeeze(nanmean(TF_par(:,1,:))), .4, 'facecolor',[0 0.6 0]); hold on;
h2 = bar(1.45:3.45, squeeze(nanmean(TF_par(:,2,:))), .4, 'facecolor',[0.7 0 0]);
errorbar(1:3,        squeeze(nanmean(TF_par(:,1,:))), squeeze(nanstd(TF_par(:,1,:)))./sqrt(nSub_eff), 'k','linestyle','none');
errorbar(1.45:3.45,  squeeze(nanmean(TF_par(:,2,:))), squeeze(nanstd(TF_par(:,2,:)))./sqrt(nSub_eff), 'k','linestyle','none');
set(gca,'xtick',1.23:3.23,'xticklabel',{'contra','ipsi','NoGo'},'xlim',[.7 3.7],'ylim',[0 45]);
ylabel('ICPC (% change)'); legend([h1 h2], {'Win','Avoid'}, 'Location','northeast'); title('Midfrontal–Motor (incl. contra)'); box off;

%% 
%% delete6
%% 6) Barplots.

% -------- A. Midfrontal–Lateral PFC --------
% First exclude the specified subjects (rows): 1, 10, 18, 21, 23, 31
sub_excl = [1 10 18 21 23 31];

TF_pfc_excl = TF_pfc;                 % Make a copy
TF_pfc_excl(sub_excl,:,:) = NaN;      % Set excluded subject rows to NaN

% Go = mean of contra & ipsi; NoGo = nogo
TF_pfc2plot = nan(size(TF_pfc_excl,1),2,2);
TF_pfc2plot(:,:,1) = nanmean(TF_pfc_excl(:,:,1:2),3);  % Go (avg contra & ipsi)
TF_pfc2plot(:,:,2) = TF_pfc_excl(:,:,3);               % NoGo

% Compute mean and SEM (NaN-aware, i.e., excluded subjects are ignored)
mean_win   = squeeze(nanmean(TF_pfc2plot(:,1,:)));                     % 1×2
sd_win     = squeeze(nanstd (TF_pfc2plot(:,1,:)));                     % 1×2
n_win      = squeeze(sum(~isnan(TF_pfc2plot(:,1,:))));                 % 1×2
sem_win    = sd_win ./ sqrt(n_win);

mean_avoid = squeeze(nanmean(TF_pfc2plot(:,2,:)));
sd_avoid   = squeeze(nanstd (TF_pfc2plot(:,2,:)));
n_avoid    = squeeze(sum(~isnan(TF_pfc2plot(:,2,:))));
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% Plot
figure('Position',[100 300 250 250]); 

h1 = bar(1:2, mean_win,   0.4, 'FaceColor',[0 0.6 0]);   hold on;  % Win
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'FaceColor',[0.7 0 0]);       % Avoid

errorbar(1:2,          mean_win,   sem_win,   'k','linestyle','none');
errorbar(1.45:2.45,    mean_avoid, sem_avoid, 'k','linestyle','none');

set(gca,'xtick',1.23:2.23,...
        'xticklabel',{'Go','NoGo'},...
        'xlim',[.7 2.7],...
        'ylim',[0 25],...
        'ytick',0:10:20);

ylabel('ICPC (% change)');
legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');
title('Midfrontal–Lateral PFC');
box off;



% -------- B. Midfrontal–Motor (exclude contralateral) --------

% Subjects to exclude
sub_excl = [1 10 18 21 23 31];

% Make a copy and set excluded subject rows to NaN
TF_par_excl = TF_par;
TF_par_excl(sub_excl,:,:) = NaN;

% Win condition (valence = 1), columns 2:3 = Go(ipsi) & NoGo
mean_win = squeeze(nanmean(TF_par_excl(:,1,2:3), 1));          % 1×2
sd_win   = squeeze(nanstd (TF_par_excl(:,1,2:3), 0, 1));       % 1×2
n_win    = squeeze(sum(~isnan(TF_par_excl(:,1,2:3)), 1));      % 1×2
sem_win  = sd_win ./ sqrt(n_win);

% Avoid condition (valence = -1), columns 2:3 = Go(ipsi) & NoGo
mean_avoid = squeeze(nanmean(TF_par_excl(:,2,2:3), 1));
sd_avoid   = squeeze(nanstd (TF_par_excl(:,2,2:3), 0, 1));
n_avoid    = squeeze(sum(~isnan(TF_par_excl(:,2,2:3)), 1));
sem_avoid  = sd_avoid ./ sqrt(n_avoid);

% Plot
figure('Position',[100 300 250 250]);
h1 = bar(1:2,       mean_win,   0.4, 'facecolor',[0 0.6 0]); hold on;
h2 = bar(1.45:2.45, mean_avoid, 0.4, 'facecolor',[0.7 0 0]);

errorbar(1:2,       mean_win,   sem_win,   'k','linestyle','none');
errorbar(1.45:2.45, mean_avoid, sem_avoid, 'k','linestyle','none');

set(gca,'xtick',1.23:2.23,...
        'xticklabel',{'Go','NoGo'},...
        'xlim',[.7 2.7],...
        'ylim',[0 25],...
        'ytick',0:10:20);

ylabel('ICPC (% change)');
legend([h1 h2], {'Win','Avoid'}, 'Location','northeast');
title('Midfrontal–Motor');
box off;


%% Compute t-weights for the lateral PFC cluster (Swart et al., 2017)
% Conditions: 1=GLW, 2=GRW, 3=GLA, 4=GRA, 5=NGW, 6=NGA
% Congruent   = (Go-Win + NoGo-Avoid)
% Incongruent = (Go-Avoid + NoGo-Win)

pfcChans = [23 24 26 27 2 3 121 122];   % Left and right PFC channels
tval_pfc = zeros(1, numel(pfcChans));   % Store t-values for each electrode

% Time window of interest (0.4–0.65 s)
iTime = dsearchn(parTF.toi4tf', [0.4 0.8]');
iTime = iTime(1):iTime(2);

for c = 1:numel(pfcChans)
    iChan = pfcChans(c);

    % ===== 1️⃣ Compute mean ICPC for each condition =====
    GoWin     = nanmean(mean(mean(tICPC(:,[1 2],iChan,:,iTime),5),4),3);  % GLW + GRW
    GoAvoid   = nanmean(mean(mean(tICPC(:,[3 4],iChan,:,iTime),5),4),3);  % GLA + GRA
    NoGoWin   = nanmean(mean(mean(tICPC(:,5,iChan,:,iTime),5),4),3);      % NGW
    NoGoAvoid = nanmean(mean(mean(tICPC(:,6,iChan,:,iTime),5),4),3);      % NGA

    % ===== 2️⃣ Compute Congruent and Incongruent means =====
    Congruent   = mean([GoWin, NoGoAvoid], 2);
    Incongruent = mean([GoAvoid, NoGoWin], 2);

    % ===== 3️⃣ Paired t-test: Incongruent vs Congruent =====
    [~,~,~,stats] = ttest(Incongruent, Congruent);
    tval_pfc(c) = stats.tstat;                     % Save t-value for each channel
end

disp('PFC cluster t-values (Incongruent > Congruent):');
disp(tval_pfc);    
% HC    0.7948    2.4114    1.2002    4.2642    3.9349    0.3577    1.0100   -1.0324 (first run)
% HC    0.7036    3.4736    0.3942    3.1620    2.9601   -0.0516    1.2290    0.1822 (second run)

% OCD   0.4255    2.1589   -0.2652    1.4101    0.8377    2.5389    1.6312    2.8189 (first run)

% OCD   0.0426    1.6083    0.0701    1.1739    0.6759    2.6471    1.3754    2.2868 (second run)

%% === Step 2: Compute t-values for the t-weighted motor cluster ===
% Conditions: 1=GLW, 2=GRW, 3=GLA, 4=GRA, 5=NGW, 6=NGA
% Left hemisphere  -> right-hand responses = contra
% Right hemisphere -> left-hand responses  = contra

lChan = [36 42 41 47];   % Left motor area (right hand = contra)
rChan = [92 97 101 102]; % Right motor area (left hand = contra)

tval_motor_L = [];  % Left hemisphere t-values
tval_motor_R = [];  % Right hemisphere t-values

%% ===== Left hemisphere (Right-hand = Contra) =====
for c = 1:numel(lChan)
    iChan = lChan(c);

    % Contra (right-hand Go): conditions 2,4
    TF4par(:,1,1) = mean(mean(nanmean(tICPC(:,2,iChan,:,iTime),5),4),3);  % contra-win (GRW)
    TF4par(:,2,1) = mean(mean(nanmean(tICPC(:,4,iChan,:,iTime),5),4),3);  % contra-avoid (GRA)

    % Ipsi (left-hand Go): conditions 1,3
    TF4par(:,1,2) = mean(mean(nanmean(tICPC(:,1,iChan,:,iTime),5),4),3);  % ipsi-win (GLW)
    TF4par(:,2,2) = mean(mean(nanmean(tICPC(:,3,iChan,:,iTime),5),4),3);  % ipsi-avoid (GLA)

    % NoGo conditions
    TF4par(:,1,3) = mean(mean(mean(tICPC(:,5,iChan,:,iTime),5),4),3);     % nogo-win (NGW)
    TF4par(:,2,3) = mean(mean(mean(tICPC(:,6,iChan,:,iTime),5),4),3);     % nogo-avoid (NGA)

    % Weight comparison formula matching the original paper
    [~,~,~,stats] = ttest( ...
        mean([2*TF4par(:,1,1) TF4par(:,2,2) TF4par(:,1,3)],2), ...  % preferred pattern
        mean([2*TF4par(:,2,1) TF4par(:,1,2) TF4par(:,2,3)],2));     % non-preferred pattern

    tval_motor_L(c) = stats.tstat;
end

disp('Left motor cluster t-values (left hemisphere, contra = right hand):');
disp(tval_motor_L);
% HC   -1.0139    0.6538    0.8044    0.2330 (first run)
% HC   -0.6664    0.4300    0.7318    0.4082 (second run)

% OCD  -1.1191    0.0210    1.1295    0.1666 (first run)
% OCD  -0.7087    0.0901    0.6992   -0.2053 (second run)

%% ===== Right hemisphere (Left-hand = Contra) =====
for c = 1:numel(rChan)
    iChan = rChan(c);

    % Contra (left-hand Go): conditions 1,3
    TF4par(:,1,1) = mean(mean(nanmean(tICPC(:,1,iChan,:,iTime),5),4),3);  % contra-win (GLW)
    TF4par(:,2,1) = mean(mean(nanmean(tICPC(:,3,iChan,:,iTime),5),4),3);  % contra-avoid (GLA)

    % Ipsi (right-hand Go): conditions 2,4
    TF4par(:,1,2) = mean(mean(nanmean(tICPC(:,2,iChan,:,iTime),5),4),3);  % ipsi-win (GRW)
    TF4par(:,2,2) = mean(mean(nanmean(tICPC(:,4,iChan,:,iTime),5),4),3);  % ipsi-avoid (GRA)

    % NoGo conditions
    TF4par(:,1,3) = mean(mean(mean(tICPC(:,5,iChan,:,iTime),5),4),3);     % nogo-win (NGW)
    TF4par(:,2,3) = mean(mean(mean(tICPC(:,6,iChan,:,iTime),5),4),3);     % nogo-avoid (NGA)

    % Same t-test
    [~,~,~,stats] = ttest( ...
        mean([2*TF4par(:,1,1) TF4par(:,2,2) TF4par(:,1,3)],2), ...
        mean([2*TF4par(:,2,1) TF4par(:,1,2) TF4par(:,2,3)],2));

    tval_motor_R(c) = stats.tstat;
end

disp('Right motor cluster t-values (right hemisphere, contra = left hand):');
disp(tval_motor_R);
% HC   1.2122    1.0770    0.5535    2.1191 (first run)
% HC   1.5457    1.3977    1.3053    2.3354 (second run)

% OCD  1.4087    0.7652    0.1929    0.0753 (first run)
% OCD  1.2397    0.8111    0.1660   -0.4212 (second run)




%%% The code below is not used
% 
% %%
% % Normalized weights
% weights_pfc   = tval_pfc ./ sum(abs(tval_pfc));
% weights_motor = tval_motor ./ sum(abs(tval_motor));
% 
% % t-weighted ICPC (example: PFC)
% tICPC_pfc_weighted = zeros(size(tICPC,1), 6, numel(parTF.freq4tf), numel(parTF.toi4tf));
% for c = 1:numel(pfcChans)
%     tICPC_pfc_weighted = tICPC_pfc_weighted + ...
%         weights_pfc(c) * squeeze(tICPC(:,:,pfcChans(c),:,:));
% end




%%%% all %%%%

%% Export for statistical analysis in SPSS.
% tICPC format: sub × condi × chan × freq × time
% Fixed condition order: GLW, GRW, GLA, GRA, NGW, NGA

clear TF_par TF_pfc

%% === Define ROI time window ===
% Example: 0.35 ~ 0.75 s (commonly used previously)
timeROI = [0.35 0.75];
iTime = dsearchn(parTF.toi4tf', timeROI');
iTime = iTime(1):iTime(2);  % Convert to index range

%% 3) ROI channels
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2  3 121 122];
lChan    = [36 42 41 47];
rChan    = [92 97 101 102];

%% ===== Parietal electrodes =====
TF_par(:,1) = mean(mean(mean(tICPC(:,1,[rChan lChan],:,iTime),5),4),3);  % GLW
TF_par(:,2) = mean(mean(mean(tICPC(:,2,[rChan lChan],:,iTime),5),4),3);  % GRW
TF_par(:,3) = mean(mean(mean(tICPC(:,3,[rChan lChan],:,iTime),5),4),3);  % GLA
TF_par(:,4) = mean(mean(mean(tICPC(:,4,[rChan lChan],:,iTime),5),4),3);  % GRA
TF_par(:,5) = mean(mean(mean(tICPC(:,5,[rChan lChan],:,iTime),5),4),3);  % NGW
TF_par(:,6) = mean(mean(mean(tICPC(:,6,[rChan lChan],:,iTime),5),4),3);  % NGA

%% ===== Parietal electrodes =====
TF_par(:,1) = mean(mean(mean(tICPC(:,1,rChan,:,iTime),5),4),3);  % GLW
TF_par(:,2) = mean(mean(mean(tICPC(:,2,rChan,:,iTime),5),4),3);  % GRW
TF_par(:,3) = mean(mean(mean(tICPC(:,3,rChan,:,iTime),5),4),3);  % GLA
TF_par(:,4) = mean(mean(mean(tICPC(:,4,rChan,:,iTime),5),4),3);  % GRA
TF_par(:,5) = mean(mean(mean(tICPC(:,5,rChan,:,iTime),5),4),3);  % NGW
TF_par(:,6) = mean(mean(mean(tICPC(:,6,rChan,:,iTime),5),4),3);  % NGA

%% ===== Parietal electrodes =====
TF_par(:,1) = mean(mean(mean(tICPC(:,1,lChan,:,iTime),5),4),3);  % GLW
TF_par(:,2) = mean(mean(mean(tICPC(:,2,lChan,:,iTime),5),4),3);  % GRW
TF_par(:,3) = mean(mean(mean(tICPC(:,3,lChan,:,iTime),5),4),3);  % GLA
TF_par(:,4) = mean(mean(mean(tICPC(:,4,lChan,:,iTime),5),4),3);  % GRA
TF_par(:,5) = mean(mean(mean(tICPC(:,5,lChan,:,iTime),5),4),3);  % NGW
TF_par(:,6) = mean(mean(mean(tICPC(:,6,lChan,:,iTime),5),4),3);  % NGA

%% ===== Lateral prefrontal electrodes =====
TF_pfc(:,1) = nanmean(mean(mean(tICPC(:,1,[lpfcChan rpfcChan],:,iTime),5),4),3);  % GLW
TF_pfc(:,2) = nanmean(mean(mean(tICPC(:,2,[lpfcChan rpfcChan],:,iTime),5),4),3);  % GRW
TF_pfc(:,3) = nanmean(mean(mean(tICPC(:,3,[lpfcChan rpfcChan],:,iTime),5),4),3);  % GLA
TF_pfc(:,4) = nanmean(mean(mean(tICPC(:,4,[lpfcChan rpfcChan],:,iTime),5),4),3);  % GRA
TF_pfc(:,5) = nanmean(mean(mean(tICPC(:,5,[lpfcChan rpfcChan],:,iTime),5),4),3);  % NGW
TF_pfc(:,6) = nanmean(mean(mean(tICPC(:,6,[lpfcChan rpfcChan],:,iTime),5),4),3);  % NGA

%% ===== Export to CSV =====
dirs.home = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\HC';
results_dir = fullfile(dirs.home, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

csvfilename = fullfile(results_dir, 'ICPC4SPSS_tweighted_GLW_GRW_GLA_GRA_NGW_NGA3.csv');

nSub = size(TF_pfc,1);
data2save = zeros(nSub, 1 + 12);  % subject + (6 pfc + 6 par)
hdr = cell(1, 1 + 12);

% === First column: subject ===
subjects = (1:nSub)';
data2save(:,1) = subjects;
hdr{1} = 'subject,';

% === Columns 2–7: six PFC conditions ===
condNames = {'GLW','GRW','GLA','GRA','NGW','NGA'};
for i = 1:6
    data2save(:, i+1) = TF_pfc(:, i);
    hdr{i+1} = sprintf('icpc_pfc_%s,', condNames{i});
end

% === Columns 8–13: six parietal conditions ===
for i = 1:6
    data2save(:, i+7) = TF_par(:, i);
    hdr{i+7} = sprintf('icpc_par_%s,', condNames{i});
end

% === Concatenate and write header ===
hdrAll = cat(2, hdr{:});
fid = fopen(csvfilename, 'w');
fprintf(fid, '%s\n', hdrAll(1:end-1));
fclose(fid);

writematrix(data2save, csvfilename, 'WriteMode', 'append');

fprintf('\nData saved to: %s\n', csvfilename);

%%%% all %%%%

%% Export for statistical analysis in SPSS.
% tICPC format: sub × condi × chan × freq × time
% Fixed condition order: GLW, GRW, GLA, GRA, NGW, NGA

clear TF_parL TF_parR TF_pfcL TF_pfcR
%%%% all %%%%

%% Export for statistical analysis in SPSS.
% tICPC format: sub × condi × chan × freq × time
% Fixed condition order: GLW, GRW, GLA, GRA, NGW, NGA

clear TF_parL TF_parR TF_parAll TF_pfc

%% === Define ROI time window ===
timeROI = [0.45 0.7];
iTime = dsearchn(parTF.toi4tf', timeROI');
iTime = iTime(1):iTime(2);

%% === ROI channel definition ===
MFchan   = [6 127];
lpfcChan = [23 24 26 27];
rpfcChan = [2 3 121 122];
lChan    = [36 42 41 47];
rChan    = [92 97 101 102];
% lChan    = [58 65 69 70];
% rChan    = [95 99];

%% ===== (1) Lateral prefrontal PFC electrodes (left-right average) =====
for i = 1:6
    TF_pfc(:, i) = mean(mean(mean(tICPC(:, i, [lpfcChan rpfcChan], :, iTime), 5), 4), 3);
end

%% ===== (2) Motor / parietal electrodes (left / right / grand average) =====
for i = 1:6
    TF_parL(:, i)   = mean(mean(mean(tICPC(:, i, lChan, :, iTime), 5), 4), 3);   % Left motor area
    TF_parR(:, i)   = mean(mean(mean(tICPC(:, i, rChan, :, iTime), 5), 4), 3);   % Right motor area
    TF_parAll(:, i) = mean(mean(mean(tICPC(:, i, [lChan rChan], :, iTime), 5), 4), 3); % Left-right average
end

%% ===== Export to CSV =====
dirs.home = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\change6';
results_dir = fullfile(dirs.home, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

csvfilename = fullfile(results_dir, 'ICPCfinal2.csv');

nSub = size(TF_pfc, 1);
condNames = {'GLW','GRW','GLA','GRA','NGW','NGA'};

% === Build export matrix ===
% subject + (6 PFC + 6 ParL + 6 ParR + 6 ParAll) = 1 + 24 columns
data2save = zeros(nSub, 1 + 24);
hdr = cell(1, 1 + 24);

% === First column: subject ===
subjects = (1:nSub)';
data2save(:,1) = subjects;
hdr{1} = 'subject,';

% === (1) Overall lateral prefrontal average ===
col = 2;
for i = 1:6
    data2save(:, col) = TF_pfc(:, i);
    hdr{col} = sprintf('icpc_pfcAll_%s,', condNames{i});
    col = col + 1;
end

% === (2) Left motor area ===
for i = 1:6
    data2save(:, col) = TF_parL(:, i);
    hdr{col} = sprintf('icpc_parL_%s,', condNames{i});
    col = col + 1;
end

% === (3) Right motor area ===
for i = 1:6
    data2save(:, col) = TF_parR(:, i);
    hdr{col} = sprintf('icpc_parR_%s,', condNames{i});
    col = col + 1;
end

% === (4) Left-right motor average ===
for i = 1:6
    data2save(:, col) = TF_parAll(:, i);
    hdr{col} = sprintf('icpc_parAll_%s,', condNames{i});
    col = col + 1;
end

% === Write header ===
hdrAll = cat(2, hdr{:});
fid = fopen(csvfilename, 'w');
fprintf(fid, '%s\n', hdrAll(1:end-1));
fclose(fid);

% === Write data ===
writematrix(data2save, csvfilename, 'WriteMode', 'append');

fprintf('\nData saved to: %s\n', csvfilename);