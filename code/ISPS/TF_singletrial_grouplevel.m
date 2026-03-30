%%%%%%%%%%%%%%%%%%%%% Vectorize all trials at once %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========================================================================
%  Trial-level midfrontal theta POWER + Inter-channel Phase Coherence (ICPC)
%  Adapted from Swart et al. (2017, DCCN, Nijmegen)
%  For OCD dataset (EEGLAB → FieldTrip)  [vectorized version: ICPC computed for all trials at once]
% =========================================================================
clear; clc; close all;
dbstop if error

%% ====== Paths and dependencies ======
data_dir = 'H:\MT_all(afterqujizhi)\cueoriandbu\error\HC';
fileList = dir(fullfile(data_dir, '*.set'));
nSub     = numel(fileList);

addpath('D:\MATLAB\toolbox\eeglab_current\eeglab2022.0_old');
addpath('D:\MATLAB\toolbox\fieldtrip-20211209\fieldtrip-20211209');
ft_defaults;

%% ====== Analysis parameters ======
freqRange = [4 8];        % theta band
timeWin   = [0.40 0.65];  % time window of interest (keep 0.40–0.65 s as requested)

% ---- Electrode group definitions ----
seedChans = {'E6','E129'};                         % Midfrontal (FCz, Cz)
pfcChans  = {'E23','E24','E26','E27','E2','E3','E123','E124'};  % bilateral PFC
lMotorCh  = {'E36','E42','E41','E47'};             % left motor area
rMotorCh  = {'E93','E98','E103','E104'};           % right motor area

% ---- Group weights (unchanged) ----
seedweights       = [5.3689, 5.0774];               % Midfrontal weights
tweights_pfc      = [0.7948 2.4114 1.2002 4.2642 3.9349 0.3577 1.0100 -1.0324];
tweights_motor_L  = [-1.0139 0.6538 0.8044 0.2330];
tweights_motor_R  = [1.2122 1.0770 0.5535 2.1191];

% Normalize weights (seed & t-weights) for vectorized use
sw = seedweights ./ sum(seedweights);
tw_pfc = tweights_pfc ./ sum(tweights_pfc);
tw_L   = tweights_motor_L ./ sum(tweights_motor_L);
tw_R   = tweights_motor_R ./ sum(tweights_motor_R);

%% ====== Initialize output matrices ======
trial2use     = cell(nSub,1);
MFpower       = cell(nSub,1);
ICPC_pfc      = cell(nSub,1);
ICPC_motor_L  = cell(nSub,1);
ICPC_motor_R  = cell(nSub,1);

%% ====== Main loop ======
for s = 1:nSub
    fprintf('\n=== Processing %s (%d/%d) ===\n', fileList(s).name, s, nSub);

    % --- Step 1: load and convert ---
    EEG = pop_loadset('filename', fileList(s).name, 'filepath', data_dir);
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');

    % --- Step 2: Laplacian (CSD) ---
    cfg = [];
    cfg.method = 'spline';
    cfg.elec   = ft_data.elec;
    scd = ft_scalpcurrentdensity(cfg, ft_data);

    nTrial = numel(scd.trial);
    trial2use{s} = ones(1,nTrial);  % Update here if trial exclusion is needed

    %% ===============================================================
    % ① Midfrontal θ Power (trial-wise) —— keep original implementation
    % ===============================================================
    cfg = [];
    cfg.channel = seedChans;
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    % Hilbert amplitude → square to get power
    cfg = [];
    cfg.hilbert = 'abs';
    hilbertData = ft_preprocessing(cfg, filtered);
    for tr = 1:nTrial
        hilbertData.trial{tr} = hilbertData.trial{tr}.^2;
    end

    % Concatenate + weighted average + zscore
    allP     = [hilbertData.trial{:}];  % chan × (time*trial)
    weights4power  = repmat(seedweights',1,size(allP,2));
    st_power = zscore(mean(weights4power .* allP,1));
    nTime    = size(hilbertData.trial{1},2);
    st_power = reshape(st_power,[nTime,nTrial]);

    % Extract mean within the time window
    iTime = dsearchn(hilbertData.time{1}',timeWin(1)) : ...
            dsearchn(hilbertData.time{1}',timeWin(2));
    MFpower{s} = squeeze(mean(st_power(iTime,:),1));  % 1×nTrial

    %% ===============================================================
    % ② Midfrontal–Lateral PFC ICPC (trial-wise) —— vectorized implementation
    % ===============================================================
    cfg = [];
    cfg.channel = [seedChans pfcChans];
    cfg.bpfilter = 'yes'; cfg.bpfreq = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    cfg = []; cfg.hilbert = 'angle';
    hilbertdata = ft_preprocessing(cfg, filtered);

    % Assemble as [nChan × nTime × nTrial]
    nChan  = numel(hilbertdata.label);
    nTime  = size(hilbertdata.trial{1},2);
    nTrial = numel(hilbertdata.trial);
    ang = reshape([hilbertdata.trial{:}], nChan, nTime, nTrial);

    % Time-window indices (same as θ power)
    iTime = dsearchn(hilbertdata.time{1}',timeWin(1)) : ...
            dsearchn(hilbertdata.time{1}',timeWin(2));

    % Compute all targets (PFC) and all trials at once:
    % diff1/2: [nTarget × |iTime| × nTrial]
    diff1_all = ang(1, iTime, :) - ang(3:end, iTime, :);
    diff2_all = ang(2, iTime, :) - ang(3:end, iTime, :);

    % PLV (mean across time) → [nTarget × 1 × nTrial]
    plv1 = abs(mean(exp(1i*diff1_all), 2));
    plv2 = abs(mean(exp(1i*diff2_all), 2));

    % Weighted combination of the two seeds → [nTarget × nTrial]
    plv_pfc = squeeze(sw(1).*plv1 + sw(2).*plv2);      % after squeeze: [nTarget × nTrial]

    % Linear combination with t-weights (signed) across targets → [1 × nTrial]
    ICPC_pfc_vals = (tw_pfc * plv_pfc);                 % 1×nTrial
    ICPC_pfc{s}   = ICPC_pfc_vals;

    %% ===============================================================
    % ③ Midfrontal–Motor ICPC (Left / Right) —— vectorized implementation
    % ===============================================================
    cfg = [];
    cfg.channel = [seedChans lMotorCh rMotorCh];
    cfg.bpfilter = 'yes'; cfg.bpfreq = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    cfg = []; cfg.hilbert = 'angle';
    hilbertdata = ft_preprocessing(cfg, filtered);

    % Assemble as [nChan × nTime × nTrial]
    nChan  = numel(hilbertdata.label);
    nTime  = size(hilbertdata.trial{1},2);
    nTrial = numel(hilbertdata.trial);
    ang = reshape([hilbertdata.trial{:}], nChan, nTime, nTrial);

    % Time-window indices
    iTime = dsearchn(hilbertdata.time{1}',timeWin(1)) : ...
            dsearchn(hilbertdata.time{1}',timeWin(2));

    % ---- Left motor (target indices 3:6) ----
    diff1_L = ang(1, iTime, :) - ang(3:6, iTime, :);    % [4 × |iTime| × nTrial]
    diff2_L = ang(2, iTime, :) - ang(3:6, iTime, :);
    plv1_L  = abs(mean(exp(1i*diff1_L), 2));            % [4 × 1 × nTrial]
    plv2_L  = abs(mean(exp(1i*diff2_L), 2));
    plv_L   = squeeze(sw(1).*plv1_L + sw(2).*plv2_L);   % [4 × nTrial]
    Lvals   = (tw_L * plv_L);                           % [1 × nTrial]
    ICPC_motor_L{s} = Lvals;

    % ---- Right motor (target indices 7:10) ----
    diff1_R = ang(1, iTime, :) - ang(7:10, iTime, :);   % [4 × |iTime| × nTrial]
    diff2_R = ang(2, iTime, :) - ang(7:10, iTime, :);
    plv1_R  = abs(mean(exp(1i*diff1_R), 2));
    plv2_R  = abs(mean(exp(1i*diff2_R), 2));
    plv_R   = squeeze(sw(1).*plv1_R + sw(2).*plv2_R);   % [4 × nTrial]
    Rvals   = (tw_R * plv_R);                           % [1 × nTrial]
    ICPC_motor_R{s} = Rvals;
end

%% ====== Save results ======
save(fullfile(data_dir,'TrialLevel_MFpower_ICPC.mat'), ...
     'MFpower','ICPC_pfc','ICPC_motor_L','ICPC_motor_R','trial2use');
fprintf('Saved trial-level θ-power and ICPC results (ICPC was computed with vectorization).\n');







%% Export the full CSV data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================== Configuration ==================
data_dir = 'H:\MT_all(afterqujizhi)\cueoriandbu\error\HC';
mat_file = fullfile(data_dir, 'TrialLevel_MFpower_ICPC.mat');
out_csv  = fullfile(data_dir, 'TrialLevel_MFpower_ICPC.csv');

%% ================== Load data ==================
S = load(mat_file, 'MFpower','ICPC_pfc','ICPC_motor_L','ICPC_motor_R','trial2use');
MFpower       = S.MFpower;
ICPC_pfc      = S.ICPC_pfc;
ICPC_motor_L  = S.ICPC_motor_L;
ICPC_motor_R  = S.ICPC_motor_R;
trial2use     = S.trial2use;

% File order consistent with the original computation (.set order in the folder)
fileList = dir(fullfile(data_dir, '*.set'));
nSub     = numel(fileList);

% Basic consistency check (optional)
if nSub ~= numel(MFpower)
    warning('fileList count (%d) does not match MFpower subject count (%d). The smaller one will be used.', nSub, numel(MFpower));
end
nSub_use = min([nSub, numel(MFpower), numel(ICPC_pfc), numel(ICPC_motor_L), numel(ICPC_motor_R), numel(trial2use)]);

%% ================== Count total rows and preallocate ==================
rowsPerSub = zeros(nSub_use,1);
for s = 1:nSub_use
    if ~isempty(MFpower{s})
        rowsPerSub(s) = numel(MFpower{s});
    else
        rowsPerSub(s) = 0;
    end
end
N = sum(rowsPerSub);

sub_col          = strings(N,1);
trial_col        = zeros(N,1);
MFpower_col      = nan(N,1);
ICPCpfc_col      = nan(N,1);
l_ICPCmotor_col  = nan(N,1);
r_ICPCmotor_col  = nan(N,1);
trial2use_col    = nan(N,1);

%% ================== Expand to long table ==================
ptr = 1;
for s = 1:nSub_use
    % Get the four measures and trial2use for this subject (all are 1×nTrial)
    mp  = MFpower{s};
    pf  = ICPC_pfc{s};
    lmo = ICPC_motor_L{s};
    rmo = ICPC_motor_R{s};
    tuse = trial2use{s};

    if isempty(mp)
        continue; % Skip empty subject
    end

    nTrial = numel(mp);

    % sub: first four characters of each filename (without extension)
    fname = fileList(s).name;              % e.g. 'OCD01.set' or 'x130.set'
    base  = erase(string(fname), '.set');  % remove extension
    if strlength(base) >= 4
        subID = extractBetween(base, 1, 4);
    else
        subID = base; % use original if length is < 4
    end

    % Write index range
    idx = ptr:(ptr + nTrial - 1);

    sub_col(idx)         = repmat(subID, nTrial, 1);
    trial_col(idx)       = (1:nTrial).';       % trial index starts from 1
    MFpower_col(idx)     = mp(:);
    ICPCpfc_col(idx)     = pf(:);
    l_ICPCmotor_col(idx) = lmo(:);
    r_ICPCmotor_col(idx) = rmo(:);

    % trial2use may be 1×nTrial or empty; set to 0 if empty
    if ~isempty(tuse)
        trial2use_col(idx) = tuse(:);
    else
        trial2use_col(idx) = 0;
    end

    ptr = ptr + nTrial;
end

%% ================== Build table and write CSV ==================
T = table( ...
    sub_col, trial_col, ...
    MFpower_col, ICPCpfc_col, l_ICPCmotor_col, r_ICPCmotor_col, trial2use_col, ...
    'VariableNames', {'sub','trial','MFpower','ICPCpfc','l_ICPCmotor','r_ICPCmotor','trial2use'} ...
);

writetable(T, out_csv);
fprintf('CSV written: %s\nTotal rows = %d (sub×trial long table)\n', out_csv, height(T));







%% OCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Vectorize all trials at once %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========================================================================
%  Trial-level midfrontal theta POWER + Inter-channel Phase Coherence (ICPC)
%  Adapted from Swart et al. (2017, DCCN, Nijmegen)
%  For OCD dataset (EEGLAB → FieldTrip)  [vectorized version: ICPC computed for all trials at once]
% =========================================================================
clear; clc; close all;
dbstop if error

%% ====== Paths and dependencies ======
data_dir = 'H:\MT_all(afterqujizhi)\cueoriandbu\error\OCD';
fileList = dir(fullfile(data_dir, '*.set'));
nSub     = numel(fileList);

addpath('D:\MATLAB\toolbox\eeglab_current\eeglab2022.0_old');
addpath('D:\MATLAB\toolbox\fieldtrip-20211209\fieldtrip-20211209');
ft_defaults;

%% ====== Analysis parameters ======
freqRange = [4 8];        % theta band
timeWin   = [0.424 0.8];  % time window of interest (keep 0.40–0.65 s as requested)

% ---- Electrode group definitions ----
seedChans = {'E6','E129'};                         % Midfrontal (FCz, Cz)
pfcChans  = {'E23','E24','E26','E27','E2','E3','E123','E124'};  % bilateral PFC
lMotorCh  = {'E36','E42','E41','E47'};             % left motor area
rMotorCh  = {'E93','E98','E103','E104'};           % right motor area

% ---- Group weights (unchanged) ----
seedweights       = [3.1478, 4.4980];               % Midfrontal weights
tweights_pfc      = [0.4255 2.1589 -0.2652 1.4101 0.8377 2.5389 1.6312 2.8189];
tweights_motor_L  = [-1.1191 0.0210 1.1295 0.1666];
tweights_motor_R  = [1.4087 0.7652 0.1929 0.0753];

% Normalize weights (seed & t-weights) for vectorized use
sw = seedweights ./ sum(seedweights);
tw_pfc = tweights_pfc ./ sum(tweights_pfc);
tw_L   = tweights_motor_L ./ sum(tweights_motor_L);
tw_R   = tweights_motor_R ./ sum(tweights_motor_R);

%% ====== Initialize output matrices ======
trial2use     = cell(nSub,1);
MFpower       = cell(nSub,1);
ICPC_pfc      = cell(nSub,1);
ICPC_motor_L  = cell(nSub,1);
ICPC_motor_R  = cell(nSub,1);

%% ====== Main loop ======
for s = 1:nSub
    fprintf('\n=== Processing %s (%d/%d) ===\n', fileList(s).name, s, nSub);

    % --- Step 1: load and convert ---
    EEG = pop_loadset('filename', fileList(s).name, 'filepath', data_dir);
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');

    % --- Step 2: Laplacian (CSD) ---
    cfg = [];
    cfg.method = 'spline';
    cfg.elec   = ft_data.elec;
    scd = ft_scalpcurrentdensity(cfg, ft_data);

    nTrial = numel(scd.trial);
    trial2use{s} = ones(1,nTrial);  % Update here if trial exclusion is needed

    %% ===============================================================
    % ① Midfrontal θ Power (trial-wise) —— keep original implementation
    % ===============================================================
    cfg = [];
    cfg.channel = seedChans;
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    % Hilbert amplitude → square to get power
    cfg = [];
    cfg.hilbert = 'abs';
    hilbertData = ft_preprocessing(cfg, filtered);
    for tr = 1:nTrial
        hilbertData.trial{tr} = hilbertData.trial{tr}.^2;
    end

    % Concatenate + weighted average + zscore
    allP     = [hilbertData.trial{:}];  % chan × (time*trial)
    weights4power  = repmat(seedweights',1,size(allP,2));
    st_power = zscore(mean(weights4power .* allP,1));
    nTime    = size(hilbertData.trial{1},2);
    st_power = reshape(st_power,[nTime,nTrial]);

    % Extract mean within the time window
    iTime = dsearchn(hilbertData.time{1}',timeWin(1)) : ...
            dsearchn(hilbertData.time{1}',timeWin(2));
    MFpower{s} = squeeze(mean(st_power(iTime,:),1));  % 1×nTrial

    %% ===============================================================
    % ② Midfrontal–Lateral PFC ICPC (trial-wise) —— vectorized implementation
    % ===============================================================
    cfg = [];
    cfg.channel = [seedChans pfcChans];
    cfg.bpfilter = 'yes'; cfg.bpfreq = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    cfg = []; cfg.hilbert = 'angle';
    hilbertdata = ft_preprocessing(cfg, filtered);

    % Assemble as [nChan × nTime × nTrial]
    nChan  = numel(hilbertdata.label);
    nTime  = size(hilbertdata.trial{1},2);
    nTrial = numel(hilbertdata.trial);
    ang = reshape([hilbertdata.trial{:}], nChan, nTime, nTrial);

    % Time-window indices (same as θ power)
    iTime = dsearchn(hilbertdata.time{1}',timeWin(1)) : ...
            dsearchn(hilbertdata.time{1}',timeWin(2));

    % Compute all targets (PFC) and all trials at once:
    % diff1/2: [nTarget × |iTime| × nTrial]
    diff1_all = ang(1, iTime, :) - ang(3:end, iTime, :);
    diff2_all = ang(2, iTime, :) - ang(3:end, iTime, :);

    % PLV (mean across time) → [nTarget × 1 × nTrial]
    plv1 = abs(mean(exp(1i*diff1_all), 2));
    plv2 = abs(mean(exp(1i*diff2_all), 2));

    % Weighted combination of the two seeds → [nTarget × nTrial]
    plv_pfc = squeeze(sw(1).*plv1 + sw(2).*plv2);      % after squeeze: [nTarget × nTrial]

    % Linear combination with t-weights (signed) across targets → [1 × nTrial]
    ICPC_pfc_vals = (tw_pfc * plv_pfc);                 % 1×nTrial
    ICPC_pfc{s}   = ICPC_pfc_vals;

    %% ===============================================================
    % ③ Midfrontal–Motor ICPC (Left / Right) —— vectorized implementation
    % ===============================================================
    cfg = [];
    cfg.channel = [seedChans lMotorCh rMotorCh];
    cfg.bpfilter = 'yes'; cfg.bpfreq = freqRange;
    filtered = ft_preprocessing(cfg, scd);

    cfg = []; cfg.hilbert = 'angle';
    hilbertdata = ft_preprocessing(cfg, filtered);

    % Assemble as [nChan × nTime × nTrial]
    nChan  = numel(hilbertdata.label);
    nTime  = size(hilbertdata.trial{1},2);
    nTrial = numel(hilbertdata.trial);
    ang = reshape([hilbertdata.trial{:}], nChan, nTime, nTrial);

    % Time-window indices
    iTime = dsearchn(hilbertdata.time{1}',timeWin(1)) : ...
            dsearchn(hilbertdata.time{1}',timeWin(2));

    % ---- Left motor (target indices 3:6) ----
    diff1_L = ang(1, iTime, :) - ang(3:6, iTime, :);    % [4 × |iTime| × nTrial]
    diff2_L = ang(2, iTime, :) - ang(3:6, iTime, :);
    plv1_L  = abs(mean(exp(1i*diff1_L), 2));            % [4 × 1 × nTrial]
    plv2_L  = abs(mean(exp(1i*diff2_L), 2));
    plv_L   = squeeze(sw(1).*plv1_L + sw(2).*plv2_L);   % [4 × nTrial]
    Lvals   = (tw_L * plv_L);                           % [1 × nTrial]
    ICPC_motor_L{s} = Lvals;

    % ---- Right motor (target indices 7:10) ----
    diff1_R = ang(1, iTime, :) - ang(7:10, iTime, :);   % [4 × |iTime| × nTrial]
    diff2_R = ang(2, iTime, :) - ang(7:10, iTime, :);
    plv1_R  = abs(mean(exp(1i*diff1_R), 2));
    plv2_R  = abs(mean(exp(1i*diff2_R), 2));
    plv_R   = squeeze(sw(1).*plv1_R + sw(2).*plv2_R);   % [4 × nTrial]
    Rvals   = (tw_R * plv_R);                           % [1 × nTrial]
    ICPC_motor_R{s} = Rvals;
end

%% ====== Save results ======
save(fullfile(data_dir,'TrialLevel_MFpower_ICPC.mat'), ...
     'MFpower','ICPC_pfc','ICPC_motor_L','ICPC_motor_R','trial2use');
fprintf('Saved trial-level θ-power and ICPC results (ICPC was computed with vectorization).\n');







clc; clear;
%% Export the full CSV data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================== Configuration ==================
data_dir = 'H:\MT_all(afterqujizhi)\cueoriandbu\error\OCD';
mat_file = fullfile(data_dir, 'TrialLevel_MFpower_ICPC.mat');
out_csv  = fullfile(data_dir, 'TrialLevel_MFpower_ICPC.csv');

%% ================== Load data ==================
S = load(mat_file, 'MFpower','ICPC_pfc','ICPC_motor_L','ICPC_motor_R','trial2use');
MFpower       = S.MFpower;
ICPC_pfc      = S.ICPC_pfc;
ICPC_motor_L  = S.ICPC_motor_L;
ICPC_motor_R  = S.ICPC_motor_R;
trial2use     = S.trial2use;

% File order consistent with the original computation (.set order in the folder)
fileList = dir(fullfile(data_dir, '*.set'));
nSub     = numel(fileList);

% Basic consistency check (optional)
if nSub ~= numel(MFpower)
    warning('fileList count (%d) does not match MFpower subject count (%d). The smaller one will be used.', nSub, numel(MFpower));
end
nSub_use = min([nSub, numel(MFpower), numel(ICPC_pfc), numel(ICPC_motor_L), numel(ICPC_motor_R), numel(trial2use)]);

%% ================== Count total rows and preallocate ==================
rowsPerSub = zeros(nSub_use,1);
for s = 1:nSub_use
    if ~isempty(MFpower{s})
        rowsPerSub(s) = numel(MFpower{s});
    else
        rowsPerSub(s) = 0;
    end
end
N = sum(rowsPerSub);

sub_col          = strings(N,1);
trial_col        = zeros(N,1);
MFpower_col      = nan(N,1);
ICPCpfc_col      = nan(N,1);
l_ICPCmotor_col  = nan(N,1);
r_ICPCmotor_col  = nan(N,1);
trial2use_col    = nan(N,1);

%% ================== Expand to long table ==================
ptr = 1;
for s = 1:nSub_use
    % Get the four measures and trial2use for this subject (all are 1×nTrial)
    mp  = MFpower{s};
    pf  = ICPC_pfc{s};
    lmo = ICPC_motor_L{s};
    rmo = ICPC_motor_R{s};
    tuse = trial2use{s};

    if isempty(mp)
        continue; % Skip empty subject
    end

    nTrial = numel(mp);

    % sub: first four characters of each filename (without extension)
    fname = fileList(s).name;              % e.g. 'OCD01.set' or 'x130.set'
    base  = erase(string(fname), '.set');  % remove extension
    if strlength(base) >= 5
        subID = extractBetween(base, 1, 5);
    else
        subID = base; % use original if length is < 4
    end

    % Write index range
    idx = ptr:(ptr + nTrial - 1);

    sub_col(idx)         = repmat(subID, nTrial, 1);
    trial_col(idx)       = (1:nTrial).';       % trial index starts from 1
    MFpower_col(idx)     = mp(:);
    ICPCpfc_col(idx)     = pf(:);
    l_ICPCmotor_col(idx) = lmo(:);
    r_ICPCmotor_col(idx) = rmo(:);

    % trial2use may be 1×nTrial or empty; set to 0 if empty
    if ~isempty(tuse)
        trial2use_col(idx) = tuse(:);
    else
        trial2use_col(idx) = 0;
    end

    ptr = ptr + nTrial;
end

%% ================== Build table and write CSV ==================
T = table( ...
    sub_col, trial_col, ...
    MFpower_col, ICPCpfc_col, l_ICPCmotor_col, r_ICPCmotor_col, trial2use_col, ...
    'VariableNames', {'sub','trial','MFpower','ICPCpfc','l_ICPCmotor','r_ICPCmotor','trial2use'} ...
);

writetable(T, out_csv);
fprintf('CSV written: %s\nTotal rows = %d (sub×trial long table)\n', out_csv, height(T));







clc; clear;
%% Generate CSV and MAT data for model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================== Configuration ==================
xlsx_file = 'H:\MT_all(afterqujizhi)\PIT_model.xlsx';
in_sheet  = 'model';
out_csv   = 'H:\MT_all(afterqujizhi)\TrialLevel_MFpower_ICPC.csv';
out_mat   = 'H:\MT_all(afterqujizhi)\TrialLevel_MFpower_ICPC.mat';
TARGET_N  = 73;   % Target output is 73×1 cell

%% ================== Load and basic check ==================
T = readtable(xlsx_file, 'Sheet', in_sheet);

% Required columns (case-insensitive)
needCols = {'sub','group','trial','stimulus','MFpower','ICPCpfc', ...
            'l_ICPCmotor','r_ICPCmotor','IVleft','IVright','trial2use'};

% Helper function for case-insensitive column access
getVar = @(tbl,name) tbl.(tbl.Properties.VariableNames{ ...
                      find(strcmpi(tbl.Properties.VariableNames,name),1)});

% Verify that all columns exist
for k = 1:numel(needCols)
    if ~any(strcmpi(T.Properties.VariableNames, needCols{k}))
        error('Column "%s" was not found (sheet %s).', needCols{k}, in_sheet);
    end
end

%% ===== Helper: safely convert any column to double; unparseable values (empty/non-numeric) -> NaN =====
numify = @(v) ( ...
    (isnumeric(v) .* double(v)) + ...
    (~isnumeric(v) .* str2double(string(v))) );
% Note: for strings/categorical/cells, use string(...) then str2double(...);
% blanks/unparseable values become NaN

%% ===== Clean the four measure columns: empty/missing -> NaN (used for both CSV and MAT) =====
MFpower_all   = numify(getVar(T,'MFpower'));
ICPCpfc_all   = numify(getVar(T,'ICPCpfc'));
lICPC_all     = numify(getVar(T,'l_ICPCmotor'));
rICPC_all     = numify(getVar(T,'r_ICPCmotor'));

%% ================== 1) Export CSV (in the specified column order) ==================
Tout = table();
for k = 1:numel(needCols)
    Tout.(needCols{k}) = getVar(T, needCols{k});
end
% Overwrite the original four columns with cleaned versions so empty values are written as NaN in CSV
Tout.MFpower      = MFpower_all;
Tout.ICPCpfc      = ICPCpfc_all;
Tout.l_ICPCmotor  = lICPC_all;
Tout.r_ICPCmotor  = rICPC_all;

writetable(Tout, out_csv);
fprintf('CSV exported: %s (%d rows × %d columns)\n', out_csv, height(Tout), width(Tout));

%% ================== 2) Build 73×1 cell structure and save MAT ==================
% —— Get and convert the other columns —— 
sub_raw        = getVar(T, 'sub');
group_raw      = getVar(T, 'group');
trial_all      = numify(getVar(T, 'trial'));
stim_all       = numify(getVar(T, 'stimulus'));
IVleft_all     = numify(getVar(T, 'IVleft'));
IVright_all    = numify(getVar(T, 'IVright'));
trial2use_all  = numify(getVar(T, 'trial2use'));

% —— Convert sub to string for grouping; also generate a numeric sub label for each row —— 
sub_str = string(sub_raw);
sub_num_all = nan(height(T),1);
tryDigits = regexp(sub_str, '(\d+)', 'match', 'once');     % extract digits, e.g. 'HC30'->'30'
hasDigits = ~cellfun(@isempty, tryDigits);
sub_num_all(hasDigits) = str2double(string(tryDigits(hasDigits)));
% For rows where digits cannot be extracted, try direct numeric conversion
isnanMask = isnan(sub_num_all);
sub_num_all(isnanMask) = str2double(sub_str(isnanMask));

% —— Map group to numeric: HC->0, OCD->1; otherwise try numeric conversion —— 
group_str = upper(strtrim(string(group_raw)));
group_num_all = nan(height(T),1);
group_num_all(group_str=="HC")  = 0;
group_num_all(group_str=="OCD") = 1;
mask_nan = isnan(group_num_all);
group_num_all(mask_nan) = str2double(group_str(mask_nan));

% —— Group by sub_str (stable order) —— 
[subs_unique, ~, G] = unique(sub_str, 'stable');   % G: group index for each row
nSub_actual = numel(subs_unique);

% Preallocate cell: first by actual subject count, then align to 73
fields = {'sub','group','trial','stimulus','MFpower','ICPCpfc', ...
          'l_ICPCmotor','r_ICPCmotor','IVleft','IVright','trial2use'};
S = struct();
for f = 1:numel(fields)
    S.(fields{f}) = cell(nSub_actual,1);
end

% —— Fill each cell by subject (1×trial row) —— 
for i = 1:nSub_actual
    idx  = (G == i);
    % Sort by trial to keep consistent order
    [trial_i, order] = sort(trial_all(idx));
    pick = find(idx);
    pick = pick(order);

    % sub: use the numeric sub label; if all NaN, fall back to group index i
    sub_num_i = sub_num_all(pick);
    if all(isnan(sub_num_i))
        sub_num_i = repmat(double(i), numel(pick), 1);
    end
    % group: numeric value (HC=0/OCD=1/otherwise try numeric conversion/still NaN -> NaN)
    group_num_i = group_num_all(pick);

    % Write into structure (convert to 1×N row); all fields are double row vectors
    S.sub{i}         = sub_num_i(:).';
    S.group{i}       = group_num_i(:).';
    S.trial{i}       = trial_i(:).';
    S.stimulus{i}    = stim_all(pick).';
    S.MFpower{i}     = MFpower_all(pick).';
    S.ICPCpfc{i}     = ICPCpfc_all(pick).';
    S.l_ICPCmotor{i} = lICPC_all(pick).';
    S.r_ICPCmotor{i} = rICPC_all(pick).';
    S.IVleft{i}      = IVleft_all(pick).';
    S.IVright{i}     = IVright_all(pick).';
    S.trial2use{i}   = trial2use_all(pick).';
end

% —— Align to 73×1 cell: pad with empty cells if fewer, keep all if more —— 
if nSub_actual < TARGET_N
    padN = TARGET_N - nSub_actual;
    for f = 1:numel(fields)
        S.(fields{f}) = [S.(fields{f}); repmat({[]}, padN, 1)];
    end
    fprintf('Actual subject count %d < %d, padded with empty cells to %d×1.\n', nSub_actual, TARGET_N, TARGET_N);
elseif nSub_actual > TARGET_N
    fprintf('Actual subject count %d > %d, saving the actual subject count without truncation.\n', nSub_actual, TARGET_N);
end

% —— Save MAT —— 
save(out_mat, '-struct', 'S');
fprintf('MAT saved: %s\nFields: %s\nSize: %d×1 cell (if >73, actual subject count is kept)\n', ...
    out_mat, strjoin(fields, ', '), max(TARGET_N, nSub_actual));