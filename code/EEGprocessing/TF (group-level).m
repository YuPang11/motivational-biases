
%% Data organization
clc;
clear;

% Set data folder paths
golefttowin = 'F:\PIT\correct_all\TF\GLW\OCD';
golefttoavoid = 'F:\PIT\correct_all\TF\GLA\OCD';
nogotowin = 'F:\PIT\correct_all\TF\NGW\OCD';
gorighttowin = 'F:\PIT\correct_all\TF\GRW\OCD';
gorighttoavoid = 'F:\PIT\correct_all\TF\GRA\OCD';
nogotoavoid = 'F:\PIT\correct_all\TF\NGA\OCD'; 

% Set output folder path
output_folder = 'F:\PIT\correct_all\TF\OCD';

% Get all .mat files in each folder
golefttowin_files = dir(fullfile(golefttowin, '*.mat'));
golefttoavoid_files = dir(fullfile(golefttoavoid, '*.mat'));
nogotowin_files = dir(fullfile(nogotowin, '*.mat'));
gorighttowin_files = dir(fullfile(gorighttowin, '*.mat'));
gorighttoavoid_files = dir(fullfile(gorighttoavoid, '*.mat'));
nogotoavoid_files = dir(fullfile(nogotoavoid, '*.mat'));

% Preallocate arrays for TF data of each condition
numFiles = length(golefttowin_files);
gtwl_sub = NaN(numFiles, 127, 39, 131);
gtal_sub = NaN(numFiles, 127, 39, 131);
ngtw_sub = NaN(numFiles, 127, 39, 131);
gtwr_sub = NaN(numFiles, 127, 39, 131);
gtar_sub = NaN(numFiles, 127, 39, 131);
ngta_sub = NaN(numFiles, 127, 39, 131);

% Loop through all subjects
for i = 1:numFiles  % Iterate over files
    % Get file names
    golefttowin_filename = fullfile(golefttowin, golefttowin_files(i).name);
    golefttoavoid_filename = fullfile(golefttoavoid, golefttoavoid_files(i).name);
    nogotowin_filename = fullfile(nogotowin, nogotowin_files(i).name);
    gorighttowin_filename = fullfile(gorighttowin, gorighttowin_files(i).name);
    gorighttoavoid_filename = fullfile(gorighttoavoid, gorighttoavoid_files(i).name);
    nogotoavoid_filename = fullfile(nogotoavoid, nogotoavoid_files(i).name);
    
    disp(['Current filename: ', golefttowin_files(i).name]);  % Show current file name

    % Load data for each condition
    try
        data_gtwl = load(golefttowin_filename);
        data_gtal = load(golefttoavoid_filename);
        data_ngtw = load(nogotowin_filename);
        data_gtwr = load(gorighttowin_filename);
        data_gtar = load(gorighttoavoid_filename);
        data_ngta = load(nogotoavoid_filename);
    catch ME
        warning('Failed to load data for subject %s: %s', golefttowin_files(i).name, ME.message);
        continue;  % Skip current file if loading fails
    end
    
    % Extract power spectrum data
    power_gtwl = data_gtwl.basefreq.powspctrm;  
    power_gtal = data_gtal.basefreq.powspctrm;  
    power_ngtw = data_ngtw.basefreq.powspctrm;  
    power_gtwr = data_gtwr.basefreq.powspctrm;  
    power_gtar = data_gtar.basefreq.powspctrm;  
    power_ngta = data_ngta.basefreq.powspctrm;  

    % Store data for each condition
    gtwl_sub(i, :, :, :) = power_gtwl;  
    gtal_sub(i, :, :, :) = power_gtal;
    ngtw_sub(i, :, :, :) = power_ngtw;
    gtwr_sub(i, :, :, :) = power_gtwr;
    gtar_sub(i, :, :, :) = power_gtar;
    ngta_sub(i, :, :, :) = power_ngta;
end

% Compute mean across relevant conditions
gtw_fin = (gtwl_sub + gtwr_sub)./2;
ngtw_fin = ngtw_sub;
gta_fin = (gtal_sub + gtar_sub)./2;
ngta_fin = ngta_sub;

% Compute congruent and incongruent conditions
congruent = (gtw_fin + ngta_sub)./2;
incongruent = (gta_fin + ngtw_sub)./2;
differ = incongruent - congruent;

% Compute grand averages across subjects
fre_gtw = squeeze(mean(gtw_fin, 1)); 
fre_ngtw = squeeze(mean(ngtw_fin, 1));
fre_gta = squeeze(mean(gta_fin, 1));
fre_ngta = squeeze(mean(ngta_fin, 1));
fre_congruent = squeeze(mean(congruent, 1));
fre_incongruent = squeeze(mean(incongruent, 1));
fre_differ = squeeze(mean(differ, 1));
fre_all = (fre_gta + fre_gtw + fre_ngta + fre_ngtw)./4;

% Save processed data
save(fullfile(output_folder, 'fre_all.mat'), 'fre_all');
save(fullfile(output_folder, 'gtwl_sub.mat'), 'gtwl_sub');
save(fullfile(output_folder, 'gtal_sub.mat'), 'gtal_sub');
save(fullfile(output_folder, 'ngtw_sub.mat'), 'ngtw_sub');
save(fullfile(output_folder, 'gtwr_sub.mat'), 'gtwr_sub');
save(fullfile(output_folder, 'gtar_sub.mat'), 'gtar_sub');
save(fullfile(output_folder, 'ngta_sub.mat'), 'ngta_sub');
save(fullfile(output_folder, 'gtw_fin.mat'), 'gtw_fin');
save(fullfile(output_folder, 'ngtw_fin.mat'), 'ngtw_fin');
save(fullfile(output_folder, 'gta_fin.mat'), 'gta_fin');
save(fullfile(output_folder, 'ngta_fin.mat'), 'ngta_fin');
save(fullfile(output_folder, 'congruent.mat'), 'congruent');
save(fullfile(output_folder, 'incongruent.mat'), 'incongruent');
save(fullfile(output_folder, 'differ.mat'), 'differ');
save(fullfile(output_folder, 'fre_gtw.mat'), 'fre_gtw');
save(fullfile(output_folder, 'fre_ngtw.mat'), 'fre_ngtw');
save(fullfile(output_folder, 'fre_gta.mat'), 'fre_gta');
save(fullfile(output_folder, 'fre_ngta.mat'), 'fre_ngta');
save(fullfile(output_folder, 'fre_congruent.mat'), 'fre_congruent');
save(fullfile(output_folder, 'fre_incongruent.mat'), 'fre_incongruent');
save(fullfile(output_folder, 'fre_differ.mat'), 'fre_differ');

disp('All data processing completed');




%% PERMUTATION TESTING (congruent & incongruent)
MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range

% Set random seed for reproducibility
rng(70);

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];  % Midfrontal channels

% Initialize TF1 and TF2
TF1 = cell(1, 36);
TF2 = cell(1, 36);

for iSub = 1:36
    % Create TF1{iSub} using condition 1 data
    TF1{iSub} = struct();
    TF1{iSub}.powspctrm = squeeze(congruent(iSub, :, :, :));
    TF1{iSub}.label = freq.label;
    TF1{iSub}.time = freq.time;
    TF1{iSub}.freq = freq.freq;
    TF1{iSub}.elec = freq.elec;
    TF1{iSub}.dimord = 'chan_freq_time';
    
    % Create TF2{iSub} using condition 2 data
    TF2{iSub} = struct();
    TF2{iSub}.powspctrm = squeeze(incongruent(iSub, :, :, :));
    TF2{iSub}.label = freq.label;
    TF2{iSub}.time = freq.time;
    TF2{iSub}.freq = freq.freq;
    TF2{iSub}.elec = freq.elec;
    TF2{iSub}.dimord = 'chan_freq_time';
end

%% Set permutation test configuration
cfg = [];
cfg.channel = MFchans;
cfg.avgoverchan = 'yes';
cfg.frequency = MFfreq;
cfg.avgoverfreq = 'yes';
cfg.latency = [0 1.3];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 1;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 500;

% Prepare neighbourhood structure
cfg_neighb.method = 'distance';
cfg_neighb.elecfile = 'easycap-M10.txt';
cfg.neighbours = ft_prepare_neighbours(cfg_neighb);

% Build design matrix
subj = size(congruent, 1);
design = zeros(2, 2 * subj);
design(1,:) = repmat(1:subj, [1, 2]);
design(2,:) = [ones(1, subj), 2 * ones(1, subj)];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

% Run permutation test
[stat] = ft_freqstatistics(cfg, TF1{:}, TF2{:});

% Show significant cluster information
if ~isempty(stat.negclusters)
    disp('Negative clusters:');
    disp([stat.negclusters.prob]);
end
if ~isempty(stat.posclusters)
    disp('Positive clusters:');
    disp([stat.posclusters.prob]);
end

%% Further analysis of significant clusters
% Get time indices of the significant negative cluster
negClusterTimeIdx = find(stat.negclusterslabelmat == 1);

% Extract time points of the significant negative cluster
tCluster1 = stat.time(negClusterTimeIdx);
disp(['Significant negative cluster time range: ', num2str(tCluster1)]);

%% Set significant cluster time window
MFtime = [0.424 0.8];
iTime = dsearchn(freq.time', MFtime');
iTime = iTime(1):iTime(2);

%% Midfrontal channels
MFchans = [6 11 127];
chanLabels = freq.label(MFchans);

%% Trial-level congruent vs incongruent difference
congruence = [];
for iSub = 1:36
    congruence(iSub,:) = mean(mean(TF2{iSub}.powspctrm(MFchans, iFreq, iTime) - ...
        TF1{iSub}.powspctrm(MFchans, iFreq, iTime), 2), 3);
end

%% Channel-wise t-test
[~, p, ~, stats] = ttest(congruence);

%% Print results
fprintf('===== Midfrontal Theta (%.3f–%.3f s, %d–%d Hz) =====\n', ...
        MFtime(1), MFtime(2), iFreq(1), iFreq(end));

for i = 1:length(MFchans)
    fprintf('%s (%d): t = %.4f, p = %.4f\n', ...
        chanLabels{i}, MFchans(i), stats.tstat(i), p(i));
end




%% resp-locked
clc;
clear;

% Set data folder paths
golefttowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\final_correct\resp\TF\GLW\OCD';
golefttoavoid = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\final_correct\resp\TF\GLA\OCD';
gorighttowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\final_correct\resp\TF\GRW\OCD';
gorighttoavoid = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\final_correct\resp\TF\GRA\OCD';

% Set output folder path
output_folder = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\final_correct\resp\TF\OCD';

% Get all .mat files
golefttowin_files = dir(fullfile(golefttowin, '*.mat'));
golefttoavoid_files = dir(fullfile(golefttoavoid, '*.mat'));
gorighttowin_files = dir(fullfile(gorighttowin, '*.mat'));
gorighttoavoid_files = dir(fullfile(gorighttoavoid, '*.mat'));

% Preallocate arrays for TF data
numFiles = length(golefttowin_files);
gtwl_sub = NaN(numFiles, 127, 37, 81);
gtal_sub = NaN(numFiles, 127, 37, 81);
gtwr_sub = NaN(numFiles, 127, 37, 81);
gtar_sub = NaN(numFiles, 127, 37, 81);

% Loop through all subjects
for i = 1:numFiles
    % Get file names
    golefttowin_filename = fullfile(golefttowin, golefttowin_files(i).name);
    golefttoavoid_filename = fullfile(golefttoavoid, golefttoavoid_files(i).name);
    gorighttowin_filename = fullfile(gorighttowin, gorighttowin_files(i).name);
    gorighttoavoid_filename = fullfile(gorighttoavoid, gorighttoavoid_files(i).name);
    
    disp(['Current filename: ', golefttowin_files(i).name]);

    % Load each condition
    try
        data_gtwl = load(golefttowin_filename);
        data_gtal = load(golefttoavoid_filename);
        data_gtwr = load(gorighttowin_filename);
        data_gtar = load(gorighttoavoid_filename);
    catch ME
        warning('Failed to load data for subject %s: %s', golefttowin_files(i).name, ME.message);
        continue;
    end
    
    % Extract power spectrum data
    power_gtwl = data_gtwl.basefreq.powspctrm;  
    power_gtal = data_gtal.basefreq.powspctrm;  
    power_gtwr = data_gtwr.basefreq.powspctrm;  
    power_gtar = data_gtar.basefreq.powspctrm;  

    % Store data for each condition
    gtwl_sub(i, :, :, :) = power_gtwl;
    gtal_sub(i, :, :, :) = power_gtal;
    gtwr_sub(i, :, :, :) = power_gtwr;
    gtar_sub(i, :, :, :) = power_gtar;
end

% Compute condition averages
gtw_fin = (gtwl_sub + gtwr_sub)./2;
gta_fin = (gtal_sub + gtar_sub)./2;

% Compute condition difference
differ = gta_fin - gtw_fin;

% Compute subject averages
fre_gtw = squeeze(mean(gtw_fin, 1)); 
fre_gta = squeeze(mean(gta_fin, 1));
fre_differ = squeeze(mean(differ, 1));

% Save processed data
save(fullfile(output_folder, 'gtwl_sub.mat'), 'gtwl_sub');
save(fullfile(output_folder, 'gtal_sub.mat'), 'gtal_sub');
save(fullfile(output_folder, 'gtwr_sub.mat'), 'gtwr_sub');
save(fullfile(output_folder, 'gtar_sub.mat'), 'gtar_sub');
save(fullfile(output_folder, 'gtw_fin.mat'), 'gtw_fin');
save(fullfile(output_folder, 'gta_fin.mat'), 'gta_fin');
save(fullfile(output_folder, 'differ.mat'), 'differ');
save(fullfile(output_folder, 'fre_gtw.mat'), 'fre_gtw');
save(fullfile(output_folder, 'fre_gta.mat'), 'fre_gta');
save(fullfile(output_folder, 'fre_differ.mat'), 'fre_differ');

disp('All data processing completed');


%% PERMUTATION TESTING (congruent & incongruent)
MFfreq = [4 8];
iFreq = dsearchn(freq.freq', MFfreq');
iFreq = iFreq(1):iFreq(2);

% Set random seed for reproducibility
rng(70);

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];

% Initialize TF1 and TF2
TF1 = cell(1, 37);
TF2 = cell(1, 37);

for iSub = 1:37
    % Create TF1{iSub} using condition 1 data
    TF1{iSub} = struct();
    TF1{iSub}.powspctrm = squeeze(gtw_fin(iSub, :, :, :));
    TF1{iSub}.label = freq.label;
    TF1{iSub}.time = freq.time;
    TF1{iSub}.freq = freq.freq;
    TF1{iSub}.elec = freq.elec;
    TF1{iSub}.dimord = 'chan_freq_time';
    
    % Create TF2{iSub} using condition 2 data
    TF2{iSub} = struct();
    TF2{iSub}.powspctrm = squeeze(gta_fin(iSub, :, :, :));
    TF2{iSub}.label = freq.label;
    TF2{iSub}.time = freq.time;
    TF2{iSub}.freq = freq.freq;
    TF2{iSub}.elec = freq.elec;
    TF2{iSub}.dimord = 'chan_freq_time';
end

%% Set permutation test configuration
cfg = [];
cfg.channel = MFchans;
cfg.avgoverchan = 'yes';
cfg.frequency = MFfreq;
cfg.avgoverfreq = 'yes';
cfg.latency = [-1 1];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 1;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 500;

% Prepare neighbourhood structure
cfg_neighb.method = 'distance';
cfg_neighb.elecfile = elec;
cfg.neighbours = ft_prepare_neighbours(cfg_neighb);

% Build design matrix
subj = 37;
design = zeros(2, 2 * subj);
design(1,:) = repmat(1:subj, [1, 2]);
design(2,:) = [ones(1, subj), 2 * ones(1, subj)];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

% Run permutation test
[stat] = ft_freqstatistics(cfg, TF1{:}, TF2{:});

% Number of clusters with p < .05
if ~isempty(stat.negclusters); disp('neg:');disp([stat.negclusters.prob]); end
if ~isempty(stat.posclusters); disp('pos'); disp([stat.posclusters.prob]); end
% Only the first negative cluster is significant (p = .002).

% Get negative cluster labels
negClusters = stat.negclusterslabelmat;
time = stat.time;

% Show time windows for all negative clusters
for iCluster = 1:max(negClusters(:))
    clusterTime = time(negClusters == iCluster);
    fprintf('Time window for negative cluster %d: %f to %f\n', iCluster, min(clusterTime), max(clusterTime));
end

tCluster1 = stat.time(stat.negclusterslabelmat==1);
disp('Time of negative clusters:');
disp(tCluster1);

%% Set significant cluster time window
MFtime = [-0.676  -0.250];
%MFtime = [-0.550  -0.126];
iTime = dsearchn(freq.time', MFtime');
iTime = iTime(1):iTime(2);

%% Midfrontal channels
MFchans = [6 11 127];
chanLabels = freq.label(MFchans);

%% Trial-level congruent vs incongruent difference
congruence = [];
for iSub = 1:36
    congruence(iSub,:) = mean(mean(TF2{iSub}.powspctrm(MFchans, iFreq, iTime) - ...
        TF1{iSub}.powspctrm(MFchans, iFreq, iTime), 2), 3);
end

%% Channel-wise t-test
[~, p, ~, stats] = ttest(congruence);

%% Print results
fprintf('===== Midfrontal Theta (%.3f–%.3f s, %d–%d Hz) =====\n', ...
        MFtime(1), MFtime(2), iFreq(1), iFreq(end));

for i = 1:length(MFchans)
    fprintf('%s (%d): t = %.4f, p = %.4f\n', ...
        chanLabels{i}, MFchans(i), stats.tstat(i), p(i));
end



%% feedback-locked
clc;
clear;

% Set data folder paths
gonopun = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\gonopun\OCD';
gonowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\gonowin\OCD';
gopun = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\gopun\OCD';
gowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\gowin\OCD';
nogonopun = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\nogonopun\OCD';
nogonowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\nogonowin\OCD'; 
nogopun = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\nogopun\OCD';
nogowin = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\nogowin\OCD'; 

% Set output folder path
output_folder = 'H:\MT_all(afterqujizhi)\cueoriandbu\feedback\TF\OCD';

% Get all .mat files in each folder
gonopun_files = dir(fullfile(gonopun, '*.mat'));
gonowin_files = dir(fullfile(gonowin, '*.mat'));
gopun_files = dir(fullfile(gopun, '*.mat'));
gowin_files = dir(fullfile(gowin, '*.mat'));
nogonopun_files = dir(fullfile(nogonopun, '*.mat'));
nogonowin_files = dir(fullfile(nogonowin, '*.mat'));
nogopun_files = dir(fullfile(nogopun, '*.mat'));
nogowin_files = dir(fullfile(nogowin, '*.mat'));

% Preallocate arrays for TF data of each condition
numFiles = length(gonopun_files);
gonopun_sub = NaN(numFiles, 127, 37, 81);
gonowin_sub = NaN(numFiles, 127, 37, 81);
gopun_sub = NaN(numFiles, 127, 37, 81);
gowin_sub = NaN(numFiles, 127, 37, 81);
nogonopun_sub = NaN(numFiles, 127, 37, 81);
nogonowin_sub = NaN(numFiles, 127, 37, 81);
nogopun_sub = NaN(numFiles, 127, 37, 81);
nogowin_sub = NaN(numFiles, 127, 37, 81);

% Loop through all subjects
for i = 1:numFiles  % Iterate over files
    % Get file names
    gonopun_filename = fullfile(gonopun, gonopun_files(i).name);
    gonowin_filename = fullfile(gonowin, gonowin_files(i).name);
    gopun_filename = fullfile(gopun, gopun_files(i).name);
    gowin_filename = fullfile(gowin, gowin_files(i).name);
    nogonopun_filename = fullfile(nogonopun, nogonopun_files(i).name);
    nogonowin_filename = fullfile(nogonowin, nogonowin_files(i).name);
    nogopun_filename = fullfile(nogopun, nogopun_files(i).name);
    nogowin_filename = fullfile(nogowin, nogowin_files(i).name);
    
    disp(['Current filename: ', gonopun_files(i).name]);  % Show current file name

    % Load data for each condition
    try
        data_gonopun = load(gonopun_filename);
        data_gonowin = load(gonowin_filename);
        data_gopun = load(gopun_filename);
        data_gowin = load(gowin_filename);
        data_nogonopun = load(nogonopun_filename);
        data_nogonowin = load(nogonowin_filename);
        data_nogopun = load(nogopun_filename);
        data_nogowin = load(nogowin_filename);
    catch ME
        warning('Failed to load data for subject %s: %s', golefttowin_files(i).name, ME.message);
        continue;  % Skip current file if loading fails
    end
    
    % Extract power spectrum data
    power_gonopun = data_gonopun.basefreq.powspctrm;  
    power_gonowin = data_gonowin.basefreq.powspctrm;  
    power_gopun = data_gopun.basefreq.powspctrm;  
    power_gowin = data_gowin.basefreq.powspctrm;  
    power_nogonopun = data_nogonopun.basefreq.powspctrm;  
    power_nogonowin = data_nogonowin.basefreq.powspctrm;  
    power_nogopun = data_nogopun.basefreq.powspctrm;  
    power_nogowin = data_nogowin.basefreq.powspctrm;  

    % Store data for each condition
    gonopun_sub(i, :, :, :) = power_gonopun;  % i is the index of the valid file
    gonowin_sub(i, :, :, :) = power_gonowin;
    gopun_sub(i, :, :, :) = power_gopun;
    gowin_sub(i, :, :, :) = power_gowin;
    nogonopun_sub(i, :, :, :) = power_nogonopun;
    nogonowin_sub(i, :, :, :) = power_nogonowin;
    nogopun_sub(i, :, :, :) = power_nogopun;
    nogowin_sub(i, :, :, :) = power_nogowin;
end

% Compute preferred and non-preferred conditions
punishment_sub = (gopun_sub + nogopun_sub)./2;
neuwin_sub = (gonowin_sub + nogonowin_sub)./2;
neuavoid_sub = (gonopun_sub + nogonopun_sub)./2;
reward_sub = (gowin_sub + nogowin_sub)./2;

nonpre_sub = (gonowin_sub + nogonowin_sub + gopun_sub + nogopun_sub)./4;
pre_sub = (gowin_sub + nogowin_sub + gonopun_sub + nogonopun_sub)./4;

differ_sub = nonpre_sub - pre_sub;

% Compute subject-level averages for all conditions
fre_gonopun = squeeze(mean(gonopun_sub, 1));
fre_gonowin = squeeze(mean(gonowin_sub, 1)); 
fre_gopun = squeeze(mean(gopun_sub, 1)); 
fre_gowin = squeeze(mean(gowin_sub, 1)); 
fre_nogonopun = squeeze(mean(nogonopun_sub, 1)); 
fre_nogonowin = squeeze(mean(nogonowin_sub, 1)); 
fre_nogopun = squeeze(mean(nogopun_sub, 1)); 
fre_nogowin = squeeze(mean(nogowin_sub, 1)); 

fre_punishment = squeeze(mean(punishment_sub, 1)); 
fre_neuwin = squeeze(mean(neuwin_sub, 1)); 
fre_neuavoid = squeeze(mean(neuavoid_sub, 1)); 
fre_reward = squeeze(mean(reward_sub, 1)); 

fre_nonpre = squeeze(mean(nonpre_sub, 1)); 
fre_pre = squeeze(mean(pre_sub, 1)); 
fre_differ = squeeze(mean(differ_sub, 1));  

% Save processed data
save(fullfile(output_folder, 'gonopun_sub.mat'), 'gonopun_sub');
save(fullfile(output_folder, 'gonowin_sub.mat'), 'gonowin_sub');
save(fullfile(output_folder, 'gopun_sub.mat'), 'gopun_sub');
save(fullfile(output_folder, 'gowin_sub.mat'), 'gowin_sub');
save(fullfile(output_folder, 'nogonopun_sub.mat'), 'nogonopun_sub');
save(fullfile(output_folder, 'nogonowin_sub.mat'), 'nogonowin_sub');
save(fullfile(output_folder, 'nogopun_sub.mat'), 'nogopun_sub');
save(fullfile(output_folder, 'nogowin_sub.mat'), 'nogowin_sub');

save(fullfile(output_folder, 'punishment_sub.mat'), 'punishment_sub');
save(fullfile(output_folder, 'neuwin_sub.mat'), 'neuwin_sub');
save(fullfile(output_folder, 'neuavoid_sub.mat'), 'neuavoid_sub');
save(fullfile(output_folder, 'reward_sub.mat'), 'reward_sub');

save(fullfile(output_folder, 'nonpre_sub.mat'), 'nonpre_sub');
save(fullfile(output_folder, 'pre_sub.mat'), 'pre_sub');
save(fullfile(output_folder, 'differ_sub.mat'), 'differ_sub');

save(fullfile(output_folder, 'fre_gonopun.mat'), 'fre_gonopun');
save(fullfile(output_folder, 'fre_gonowin.mat'), 'fre_gonowin');
save(fullfile(output_folder, 'fre_gopun.mat'), 'fre_gopun');
save(fullfile(output_folder, 'fre_gowin.mat'), 'fre_gowin');
save(fullfile(output_folder, 'fre_nogonopun.mat'), 'fre_nogonopun');
save(fullfile(output_folder, 'fre_nogonowin.mat'), 'fre_nogonowin');
save(fullfile(output_folder, 'fre_nogopun.mat'), 'fre_nogopun');
save(fullfile(output_folder, 'fre_nogowin.mat'), 'fre_nogowin');

save(fullfile(output_folder, 'fre_punishment.mat'), 'fre_punishment');
save(fullfile(output_folder, 'fre_neuwin.mat'), 'fre_neuwin');
save(fullfile(output_folder, 'fre_neuavoid.mat'), 'fre_neuavoid');
save(fullfile(output_folder, 'fre_reward.mat'), 'fre_reward');

save(fullfile(output_folder, 'fre_nonpre.mat'), 'fre_nonpre');
save(fullfile(output_folder, 'fre_pre.mat'), 'fre_pre');
save(fullfile(output_folder, 'fre_differ.mat'), 'fre_differ');

disp('All data processing completed');


%%
%% Plot TF map
presti_all_ROI = [6 11 127]; % Define ROI electrodes
figure;
x = freq.time;  
y = freq.freq;  
mu_ngtw = squeeze(mean(fre_pre(presti_all_ROI, :, :), 1));

imagesc(x, y, mu_ngtw(:,:));  % mean_freqs should match the dimensions of freq.time and freq.freq
axis xy;  
colorbar;  

%set(gca, 'clim', [-2 -1], 'xlim', [-0.25 3], 'yscale', 'log', 'ytick', [2 4 8 16 32]);
set(gca,  'xlim', [-1 1], 'yscale', 'log', 'ytick', [2 4 8 16 32]);

% Set x-axis range and ticks
xlim([-0.25 1]);  % Set x-axis range
xticks(-0.25:0.25:1);  % Set x-axis ticks every 0.25

title(['nonpre-pre'], 'fontsize', 16);
xlabel('Latence (ms)', 'fontsize', 16);
ylabel('Frequency (Hz)', 'fontsize', 16);

hold on;  % Keep current figure
contourf(x, y, mu_ngtw, 36, 'linestyle', 'none');  % Draw contour plot




% CONTRAST PLOTS.
figure('name','midfrontal TF power Valence x Required action','Position',[100 500 1000 250])
colormap('jet')
MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range
MFtime = [0.15 0.65];
iTime = dsearchn(freq.time',MFtime');
iTime = iTime(1):iTime(2);

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];

%% Difference topoplot
chan2plot = [6 11 127];
figure('name','midfrontal TF power Valence x Required action','Position',[100 500 500 200]) % Increase width
colormap('jet')

% Build freq_differ structure and ensure it is a valid freq datatype
freq_differ = struct();
freq_differ.powspctrm = fre_differ;  % Spectral data
freq_differ.time = freq.time;        % Time information
freq_differ.freq = freq.freq;        % Frequency information
freq_differ.label = freq.label;      % Electrode labels
freq_differ.elec = freq.elec;        % Electrode information

% Draw topoplot with ft_topoplotTFR
cfg = []; 
cfg.ylim = MFfreq; 
cfg.zlim = [-1 1]; 
cfg.marker = 'on'; 
cfg.style = 'straight'; 

% Manually adjust electrode coordinates
elec = freq.elec;  
elec.chanpos(:, 1) = -freq.elec.chanpos(:, 2);  % Reverse x-coordinate
elec.chanpos(:, 2) = freq.elec.chanpos(:, 1);   % Set y-coordinate to original x-coordinate
cfg.elec = elec;  % Update electrode layout
cfg.comment = 'no'; 
cfg.xlim = MFtime;
cfg.highlight = 'on'; 
cfg.highlightchannel = MFchans; 
cfg.highlightsymbol = 'o';

% Draw topoplot
ft_topoplotTFR(cfg, freq_differ);


%% Difference TF plot
plot_differ = squeeze(mean(fre_differ(chan2plot,:,:), 1));
contourf(freq.time, freq.freq, plot_differ, 36, 'linestyle', 'none'); hold on
set(gca, 'ytick', [2 4 8 16 32], 'yscale', 'log', 'xtick', [0 .5 1],...
    'ylim', [3 50], 'fontsize', 11, 'clim', [-1.2 1.2], 'xlim', [-0.25 1]);
box off;
plot([MFtime(2) MFtime(2)], MFfreq, 'k');
plot([MFtime(1) MFtime(1)], MFfreq, 'k');
plot(MFtime, [MFfreq(1) MFfreq(1)], 'k');
plot(MFtime, [MFfreq(2) MFfreq(2)], 'k');
colorbar;
ylabel('Frequency (Hz)', 'fontsize', 11, 'fontweight', 'bold');
xlabel('Time (s)', 'fontsize', 11, 'fontweight', 'bold');


clc;
clear;


%% PERMUTATION TESTING (congruent & incongruent)
MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range

% Set random seed for reproducibility
rng(70);

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];

% Initialize TF1 and TF2
TF1 = cell(1, 36);  % 27 subjects, each stored as a structure
TF2 = cell(1, 36);

for iSub = 1:36
    % Create TF1{iSub} using condition 1 data
    TF1{iSub} = struct();
    TF1{iSub}.powspctrm = squeeze(nonpre_sub(iSub, :, :, :));
    TF1{iSub}.label = freq.label;
    TF1{iSub}.time = freq.time;
    TF1{iSub}.freq = freq.freq;
    TF1{iSub}.elec = freq.elec;
    TF1{iSub}.dimord = 'chan_freq_time';
    
    % Create TF2{iSub} using condition 2 data
    TF2{iSub} = struct();
    TF2{iSub}.powspctrm = squeeze(pre_sub(iSub, :, :, :));
    TF2{iSub}.label = freq.label;
    TF2{iSub}.time = freq.time;
    TF2{iSub}.freq = freq.freq;
    TF2{iSub}.elec = freq.elec;
    TF2{iSub}.dimord = 'chan_freq_time';
end

%% Set permutation test configuration
cfg = [];
cfg.channel = MFchans;  % Channels of interest
cfg.avgoverchan = 'yes';  % Average over channels
cfg.frequency = MFfreq;  % Frequency band (Theta: 4-8 Hz)
cfg.avgoverfreq = 'yes';  % Average over frequency
cfg.latency = [0 1];  % Time window
cfg.method = 'montecarlo';  % Use Monte Carlo method
cfg.statistic = 'ft_statfun_depsamplesT';  % Paired t-test
cfg.correctm = 'cluster';  % Cluster correction
cfg.clusteralpha = 0.05;  % Cluster significance level
cfg.clusterstatistic = 'maxsum';  % Cluster statistic method
cfg.minnbchan = 1;  % Minimum neighboring channels
cfg.tail = 0;  % Two-tailed test
cfg.clustertail = 0;  % Two-tailed cluster test
cfg.alpha = 0.05;  % Significance level
cfg.numrandomization = 500;  % Number of randomizations

% Prepare neighbourhood structure for cluster correction
cfg_neighb.method = 'distance';  % Define neighbours by distance
cfg_neighb.elecfile = elec;  % Electrode layout
cfg.neighbours = ft_prepare_neighbours(cfg_neighb);  % Prepare neighbourhood info

% Build design matrix (2 conditions)
subj = 36;  % Number of subjects
design = zeros(2, 2 * subj);  % 2 rows, 2 * subj columns
design(1,:) = repmat(1:subj, [1, 2]);  % Subject IDs
design(2,:) = [ones(1, subj), 2 * ones(1, subj)];  % Conditions

cfg.design = design;
cfg.uvar = 1;  % Subject variable
cfg.ivar = 2;  % Condition variable

% Run permutation test
[stat] = ft_freqstatistics(cfg, TF1{:}, TF2{:});

% Number of clusters with p < .05
if ~isempty(stat.negclusters); disp('neg:');disp([stat.negclusters.prob]); end
if ~isempty(stat.posclusters); disp('pos'); disp([stat.posclusters.prob]); end
% Only the first negative cluster is significant (p = .002).

posClusters = stat.posclusterslabelmat;  % Positive cluster label matrix
for iCluster = 1:max(posClusters(:))  % Get time window of each positive cluster
    clusterTime = stat.time(posClusters == iCluster);
    fprintf('Time window for positive cluster %d: %f to %f\n', iCluster, min(clusterTime), max(clusterTime));
end

% Time of significant cluster
tCluster1 = stat.time(stat.negclusterslabelmat==1); % ranges from -.826 to -.15.
MFtime = [0.2 0.45];
iTime = dsearchn(freq.time',MFtime');
iTime = iTime(1):iTime(2);

% Which electrodes contribute to this effect? Post-hoc alpha = .05/3 = .017
congruence = [];
for iSub = 1:36
    congruence(iSub,:) = mean(mean(TF2{iSub}.powspctrm(MFchans,iFreq,iTime)- ...
        TF1{iSub}.powspctrm(MFchans,iFreq,iTime),3),2);
end
[~,p] = ttest(congruence) % all electrodes [6 11 127] show the effect significantly.


%% Plotting of permutation results

% Two-condition line plot
for iSub = 1:36
    non_preferred(iSub,:) = mean(mean(TF1{iSub}.powspctrm(MFchans,iFreq,:),2));
    preferred(iSub,:) = mean(mean(TF2{iSub}.powspctrm(MFchans,iFreq,:),2));
end
figure
plot(freq.time,mean(preferred),'b','linewidth',2); hold on
plot(freq.time,mean(non_preferred),'r','linewidth',2); 
plot(freq.time(iTime),mean(preferred(:,iTime)),'linewidth',5); 
plot(freq.time(iTime),mean(non_preferred(:,iTime)),'r','linewidth',5); 
ylabel('Power (dB)','fontsize',11,'fontweight','bold')
title('Midfrontal [6 11 127] theta (4-8Hz) power','fontsize',11)
box off; hold on
plot([2 2], get(gca,'ylim'),'k')
plot([3 3], get(gca,'ylim'),'k')
set(gca,'xtick',[0 1],'xtickLabel',{'FB','ITI'},...
    'ylim',[-1 2.5],'fontsize',11,'xlim',[0 1])
plot(MFtime,[0.8 0.8],'k','LineWidth',2) % Increase line width
% Add ** above the horizontal line
text(mean(MFtime), 0.5, 'p=0.01(0.2s-0.45s)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
legend('preferred','non-preferred'); legend boxoff


% Four-condition line plot: reward / punishment / neutral(2x) outcomes
%% ineplot of separate conditions.(gta_fin\gtw_fin\ngta_fin\ngtw_fin)

chan2plot = [6 11 127]; 

MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range

MFtime = [0.15 0.576];
iTime = dsearchn(freq.time',MFtime');
iTime = iTime(1):iTime(2);

plot_punishment = squeeze(mean(mean(mean(punishment_sub(:,chan2plot,iFreq,:), 1), 2), 3));
plot_neuwin = squeeze(mean(mean(mean(neuwin_sub(:,chan2plot,iFreq,:),1),2),3)) ;
plot_neuavoid = squeeze(mean(mean(mean(neuavoid_sub(:,chan2plot,iFreq,:),1),2),3)) ;
plot_reward = squeeze(mean(mean(mean(reward_sub(:,chan2plot,iFreq,:),1),2),3)) ;

figure('Position',[800 500 450 250]); hold on
plot(freq.time,plot_punishment,'-','linewidth',2,'Color',[.8 0 0])
plot(freq.time, plot_neuwin, '-', 'linewidth', 2, 'Color', [0.5 1 0.5])   
plot(freq.time, plot_neuavoid, '-', 'linewidth', 2, 'Color', [1 0.5 0.5])  
plot(freq.time,plot_reward,'-','linewidth',2,'Color',[0 .5 0])
set(gca,'xtick',[0 1],'ylim',[-1 3],'xtickLabel',{'FB','1'},'fontsize',11,'xlim',[0 1])
legend({'punishment','neutral (Win cues)','neutral (Avoid cue)','reward'},'Location','NorthEastOutside' ); 
legend boxoff 
plot(MFtime, [0.8 0.8], 'k', 'linewidth', 2, 'HandleVisibility', 'off');
% Add ** above the horizontal line
text(mean(MFtime), 0.9, '**', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ylabel('Power (dB)','fontweight','bold'); xlabel('Time (s)','fontweight','bold')




%%
% Plot only the difference TF map
figure('name','midfrontal TF power Bad-Good outcome','Position',[200 500 400 300])

% Define midfrontal channels
MFchans = [6 11 127];  % Midfrontal channels
MFtime  = [0.2 0.45];
MFfreq  = [4 8];  % Theta band

colormap jet

% ===== TF map =====
dat2plot = squeeze(mean(fre_differ(MFchans,:,:), 1));  % Average across midfrontal channels
title('Non-preferred - preferred outcome','fontsize',11); hold on
contourf(freq.time, freq.freq, dat2plot, 36, 'linestyle','none');

set(gca,'ytick',[4 8 16 32],'yscale','log','xlim',[-0.1 1],...
    'xtick',[0 1],'xtickLabel',{'FB','1'},'ylim',[3 50],...
    'fontsize',11,'clim',[-1 1])

% Mark analysis time window and frequency range (rectangle)
rectangle('Position',[MFtime(1), MFfreq(1), diff(MFtime), diff(MFfreq)], ...
          'EdgeColor','k','LineWidth',1.5);

colorbar
ylabel('Frequency (Hz)','fontsize',11,'fontweight','bold');
xlabel('Time (s)','fontsize',11,'fontweight','bold');




%% Difference topoplot
% Define midfrontal channels
MFchans = [6 11 127];  % Midfrontal channels
figure('name','midfrontal TF power Valence x Required action','Position',[100 500 500 200])

cmap = jet(64);        % Keep fully consistent with above
colormap(cmap);

% Build freq_differ structure
freq_differ = struct();
freq_differ.powspctrm = fre_differ;
freq_differ.time      = freq.time;
freq_differ.freq      = freq.freq;
freq_differ.label     = freq.label;
freq_differ.elec      = freq.elec;

cfg = [];
cfg.ylim      = MFfreq;
cfg.xlim      = MFtime;
cfg.zlim      = [-1.5 1.5];   % Match the TF plot clim
cfg.marker    = 'on';
cfg.style     = 'straight';
cfg.comment   = 'no';
cfg.highlight = 'on';
cfg.highlightchannel = MFchans;
cfg.highlightsymbol  = 'o';

cfg.colormap = cmap;  % Key step: make ft_topoplotTFR use the same colormap

% Manually adjust electrode coordinates
elec = freq.elec;
elec.chanpos(:,1) = -freq.elec.chanpos(:,2);
elec.chanpos(:,2) =  freq.elec.chanpos(:,1);
cfg.elec = elec;

% Draw topo
ft_topoplotTFR(cfg, freq_differ);

% Extra safeguard: force matching color range and colormap
caxis([-1.5 1.5]);    % Equivalent to set(gca,'clim',[-0.6 0.6]);
colorbar;
set(gca,'fontsize',11);





%% OCD
% Plot bar chart for eight conditions
chan2plot = [6 11 127]; 

MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range

MFtime = [0.2 0.45];
iTime = dsearchn(freq.time',MFtime');
iTime = iTime(1):iTime(2);

gonowin = zeros(36, 1);  % Store power for gonowin
nogonowin = zeros(36, 1);  % Store power for nogonowin
gopun = zeros(36, 1);  % Store power for gopun
nogopun = zeros(36, 1);  % Store power for nogopun
gowin = zeros(36, 1);  % Store power for gowin
nogowin = zeros(36, 1);  % Store power for nogowin
gonopun = zeros(36, 1);  % Store power for gonopun
nogonopun = zeros(36, 1);  % Store power for nogonopun

for iSub = 1:36
    gonowin(iSub) = squeeze(mean(mean(mean(gonowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogonowin(iSub) = squeeze(mean(mean(mean(nogonowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gopun(iSub) = squeeze(mean(mean(mean(gopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gowin(iSub) = squeeze(mean(mean(mean(gowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gonopun(iSub) = squeeze(mean(mean(mean(gonopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogonopun(iSub) = squeeze(mean(mean(mean(nogonopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
end
valxaccAction = [gonowin, nogonowin, gopun, nogopun, gowin, nogowin, gonopun, nogonopun];  % Combined data

%% Assume the following 8 vectors have already been computed:
% gonowin, nogonowin, gopun, nogopun, gowin, nogowin, gonopun, nogonopun
% Each is [nSub × 1], here nSub = 36

% Subject IDs (replace with actual IDs if available)
nSub    = size(gonowin, 1);
subject = (1:nSub)';

% Build table with variable names matching the SPSS columns:
% subject, goNoWin, nogoNoWin, goPun, nogoPun, goWin, nogoWin, goNoPun, nogoNoPun
T = table(subject, ...
          gonowin, ...
          nogonowin, ...
          gopun, ...
          nogopun, ...
          gowin, ...
          nogowin, ...
          gonopun, ...
          nogonopun, ...
          'VariableNames', {'subject', ...
                            'goNoWin', 'nogoNoWin', ...
                            'goPun',   'nogoPun', ...
                            'goWin',   'nogoWin', ...
                            'goNoPun', 'nogoNoPun'});

% Output path
outDir  = 'H:\MT_all(afterqujizhi)\feedback\feedback2\OCD';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outFile = fullfile(outDir, 'fb4spss.csv');

% Write CSV
writetable(T, outFile);   % Default: include header, comma-separated

fprintf('Saved to: %s\n', outFile);


%% HC
% Plot bar chart for eight conditions
chan2plot = [6 11 127]; 

MFfreq = [4 8];  % Theta band
iFreq = dsearchn(freq.freq', MFfreq');  % Find frequency range indices
iFreq = iFreq(1):iFreq(2);  % Frequency index range

MFtime = [0.15 0.576];
iTime = dsearchn(freq.time',MFtime');
iTime = iTime(1):iTime(2);

gonowin = zeros(37, 1);  % Store power for gonowin
nogonowin = zeros(37, 1);  % Store power for nogonowin
gopun = zeros(37, 1);  % Store power for gopun
nogopun = zeros(37, 1);  % Store power for nogopun
gowin = zeros(37, 1);  % Store power for gowin
nogowin = zeros(37, 1);  % Store power for nogowin
gonopun = zeros(37, 1);  % Store power for gonopun
nogonopun = zeros(37, 1);  % Store power for nogonopun

for iSub = 1:37
    gonowin(iSub) = squeeze(mean(mean(mean(gonowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogonowin(iSub) = squeeze(mean(mean(mean(nogonowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gopun(iSub) = squeeze(mean(mean(mean(gopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gowin(iSub) = squeeze(mean(mean(mean(gowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    gonopun(iSub) = squeeze(mean(mean(mean(gonopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
    nogonopun(iSub) = squeeze(mean(mean(mean(nogonopun_sub(iSub,chan2plot,iFreq,iTime), 2), 3), 4));
end
valxaccAction = [gonowin, nogonowin, gopun, nogopun, gowin, nogowin, gonopun, nogonopun];  % Combined data

%% Assume the following 8 vectors have already been computed:
% gonowin, nogonowin, gopun, nogopun, gowin, nogowin, gonopun, nogonopun
% Each is [nSub × 1], here nSub = 36

% Subject IDs (replace with actual IDs if available)
nSub    = size(gonowin, 1);
subject = (1:nSub)';

% Build table with variable names matching the SPSS columns:
% subject, goNoWin, nogoNoWin, goPun, nogoPun, goWin, nogoWin, goNoPun, nogoNoPun
T = table(subject, ...
          gonowin, ...
          nogonowin, ...
          gopun, ...
          nogopun, ...
          gowin, ...
          nogowin, ...
          gonopun, ...
          nogonopun, ...
          'VariableNames', {'subject', ...
                            'goNoWin', 'nogoNoWin', ...
                            'goPun',   'nogoPun', ...
                            'goWin',   'nogoWin', ...
                            'goNoPun', 'nogoNoPun'});

% Output path
outDir  = 'H:\MT_all(afterqujizhi)\feedback\feedback2\HC';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outFile = fullfile(outDir, 'fb4spss.csv');

% Write CSV
writetable(T, outFile);   % Default: include header, comma-separated

fprintf('Saved to: %s\n', outFile);



figure;
hold on;

% Draw the first four bars (red)
bar(1:4, mean(valxaccAction(:, 1:4)), 'FaceColor', [1 0 0]);  % Red

% Draw the last four bars (blue)
bar(5:8, mean(valxaccAction(:, 5:8)), 'FaceColor', [0 0 1]);  % Blue

% Draw error bars
errorbar(1:8, mean(valxaccAction), std(valxaccAction)./sqrt(37), 'k', 'LineStyle', 'none');

% Draw a black line between bars 1 and 4
plot([1 4], [max(mean(valxaccAction)) + 0.3, max(mean(valxaccAction)) + 0.3], 'k', 'LineWidth', 2);
text(2.5, max(mean(valxaccAction)) + 0.4, 'Non-preferred', 'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold', 'FontSize', 14);  % Larger font

% Draw a black line between bars 5 and 8
plot([5 8], [max(mean(valxaccAction)) + 0.05, max(mean(valxaccAction)) + 0.05], 'k', 'LineWidth', 2);
text(6.5, max(mean(valxaccAction)) + 0.15, 'Preferred', 'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold', 'FontSize', 14);  % Larger font

% Set figure labels
ylabel('Power (dB)', 'FontWeight', 'bold', 'FontSize', 22);
set(gca, 'XTick', 1:8, 'XTickLabel', {'goNoWin', 'nogoNoWin', 'goPun', 'nogoPun', ...
    'goWin', 'nogoWin', 'goNoPun', 'nogoNoPun'}, 'FontSize', 12, 'FontWeight', 'bold');

hold off;



clc;
clear;
%% Relative enhanced vs. reduced learning contrast TF map
MFchans = [6 11 127];   % Midfrontal channels
MFfreq  = [4 8];        % Theta band
iFreq   = dsearchn(freq.freq', MFfreq');  
iFreq   = iFreq(1):iFreq(2);  

MFtime  = [0.2  0.45];   % Feedback time window
iTime   = dsearchn(freq.time',MFtime');
iTime   = iTime(1):iTime(2);

figure('name','midfrontal TF power relative enhanced vs. reduced learning','Position',[200 500 700 500])
colormap jet

% ===== Build FieldTrip freq data structure =====
data = struct();
data.time  = freq.time;   
data.freq  = freq.freq;   
data.label = freq.label;  
data.elec  = freq.elec;   

% ===== Build contrast formula: (GoWin - NoGoWin) + (GoPun - NoGoPun) =====
contrast_data = (fre_gowin - fre_nogowin) + (fre_gopun - fre_nogopun);
data.powspctrm = contrast_data ./ 2;   % Divide by 2 to keep the same magnitude

% ===== Midfrontal average =====
dat2plot = squeeze(mean(data.powspctrm(MFchans,:,:), 1));

% ===== Draw TF map =====
contourf(freq.time, freq.freq, dat2plot, 40, 'linestyle','none')
set(gca,'ytick',[2 4 8 16 32],'yscale','log','xlim',[-0.1 1],...
    'xtick',[0 1],'xtickLabel',{'FB','1'},'ylim',[3 50],...
    'fontsize',11,'clim',[-1 1])

% ===== Mark ROI =====
rectangle('Position',[MFtime(1), MFfreq(1), diff(MFtime), diff(MFfreq)],...
          'EdgeColor','k','LineWidth',1.5)

colorbar
title('Relative enhanced vs. reduced learning','fontsize',12)
ylabel('Frequency (Hz)','fontsize',11,'fontweight','bold')
xlabel('Time (s)','fontsize',11,'fontweight','bold')



%% Relative enhanced vs. reduced learning contrast topoplot

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];

figure('name','midfrontal TF power Valence x Required action','Position',[100 500 500 200]) % Increase width

% Use the same colormap and color range as above
cmap = jet(64);
colormap(cmap);

% Draw topoplot with ft_topoplotTFR
cfg = []; 
cfg.ylim   = MFfreq; 
cfg.zlim   = [-1.5 1.5];      % Keep consistent with the above section
cfg.marker = 'on'; 
cfg.style  = 'straight'; 

% Manually adjust electrode coordinates
elec = freq.elec;
elec.chanpos(:, 1) = -freq.elec.chanpos(:, 2);  % Reverse x-coordinate
elec.chanpos(:, 2) =  freq.elec.chanpos(:, 1);  % Set y-coordinate to original x-coordinate
cfg.elec  = elec;                                % Update electrode layout

cfg.comment          = 'no'; 
cfg.xlim             = MFtime;
cfg.highlight        = 'on'; 
cfg.highlightchannel = MFchans; 
cfg.highlightsymbol  = 'o';

cfg.colormap = cmap;   % Make ft_topoplotTFR use the same colormap

% Draw topo
ft_topoplotTFR(cfg, data);

% Force the same color range and show colorbar
caxis([-1.5 1.5]);      % Equivalent to set(gca,'clim',[-1.5 1.5]);
colorbar;
set(gca,'fontsize',11);


%% Compute the difference between relative enhanced vs. reduced learning and run one-sample t-test, then plot

gopun = zeros(36, 1);  % Store power for gopun
nogopun = zeros(36, 1);  % Store power for nogopun
gowin = zeros(36, 1);  % Store power for gowin
nogowin = zeros(36, 1);  % Store power for nogowin

for iSub = 1:36   
    gopun(iSub) = squeeze(mean(mean(mean(gopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    gowin(iSub) = squeeze(mean(mean(mean(gowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));    
end

% ===== Build contrast: (GoWin + GoPun) - (NoGoWin + NoGoPun) =====
contrast_eachSub = (gowin + gopun) - (nogowin + nogopun);

% ===== One-sample t-test (test whether the contrast differs from 0) =====
[~, p] = ttest(contrast_eachSub);

% ===== Compute mean and standard error (for plotting) =====
means = mean(contrast_eachSub);
sems  = std(contrast_eachSub) / sqrt(37);

% ===== Draw bar plot =====
figure;
b = bar(means, 'FaceColor', 'flat'); 

% Use a soft color
b.CData = [0.1 0.5 0.9];  % Lake blue

ylim([-1 3]);  
alpha(0.85);
hold on;

% Add error bar
errorbar(1, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 8);

% Add p-value label and horizontal line
y_line = 2;  
y_pval = 2.1; 
plot([0.8, 1.2], [y_line, y_line], 'k', 'LineWidth', 1.5);  
text(1, y_pval, sprintf('p = %.3f', p), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');  

% Beautify axes
xticks(1);
xticklabels({'Contrast (Enhanced - Reduced)'});
ylabel('Theta Power (dB)', 'FontWeight', 'bold', 'FontSize', 16);
title('Relative Enhanced vs. Reduced Learning (Single-Sample T-Test)', ...
      'FontWeight', 'bold', 'FontSize', 16);
set(gca,'FontSize',14,'FontWeight','bold','XColor','k','YColor','k');
grid off;
hold off;



%% Compute t-test between relative enhanced vs. reduced learning and draw bar plot
gopun   = zeros(36, 1);
nogopun = zeros(36, 1);
gowin   = zeros(36, 1);
nogowin = zeros(36, 1);

% Compute mean power for each subject
for iSub = 1:36   
    gopun(iSub)   = squeeze(mean(mean(mean(gopun_sub(iSub,MFchans,iFreq,iTime),   2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    gowin(iSub)   = squeeze(mean(mean(mean(gowin_sub(iSub,MFchans,iFreq,iTime),   2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));    
end

% ===== Build difference scores =====
valxaccAction1 = gowin - nogowin;   % Enhanced Learning
valxaccAction2 = gopun - nogopun;   % Reduced Learning

% ===== Compute means and standard errors (SE) =====
means = [mean(valxaccAction1), mean(valxaccAction2)];
sems  = [std(valxaccAction1)/sqrt(37), std(valxaccAction2)/sqrt(37)];

% ===== Paired-sample t-test (Enhanced vs. Reduced) =====
[~, p] = ttest(valxaccAction1, valxaccAction2);

%% ===== Draw bar plot =====
figure;
b = bar(means, 'FaceColor', 'flat'); 
b.CData = [0.1 0.5 0.9;    % Lake blue (Enhanced Learning)
           0.9 0.4 0.3];   % Coral red (Reduced Learning)
ylim([-0.5 0.5]);  
alpha(0.85);
hold on;

% Error bars
errorbar(1:2, means, sems, 'k', 'LineStyle','none', 'LineWidth', 2, 'CapSize', 8);

% Draw connecting line and p-value
y_line = max(means + sems) + 0.2;
y_pval = y_line + 0.1;
plot([1, 2], [y_line, y_line], 'k', 'LineWidth', 1.5);  
text(1.5, y_pval, sprintf('p = %.3f', p), ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment','center');  

% Beautify axes
xticks([1 2]);
xticklabels({'Enhanced Learning', 'Reduced Learning'});
ylabel('Theta Power (dB)', 'FontWeight','bold','FontSize',14);
title('Comparison of Relative Enhanced vs. Reduced Learning', 'FontWeight','bold','FontSize',16);
set(gca,'FontSize',12,'FontWeight','bold');
grid off;
hold off;



clc;
clear;
%% Biased vs. unbiased learning contrast TF map
MFchans = [6 11 127];   % Midfrontal channels
MFfreq  = [4 8];        % Theta band
iFreq   = dsearchn(freq.freq', MFfreq');  
iFreq   = iFreq(1):iFreq(2);  

MFtime  = [0.2 0.45];   % Feedback time window
iTime   = dsearchn(freq.time',MFtime');
iTime   = iTime(1):iTime(2);

figure('name','midfrontal TF power biased-unbiased learning','Position',[200 500 700 500])
colormap jet

% ===== Build FieldTrip freq data structure =====
data = struct();
data.time  = freq.time;   
data.freq  = freq.freq;   
data.label = freq.label;  
data.elec  = freq.elec;   

% ===== Build contrast formula: (GoWin + NoGoPun) - (NoGoWin + GoPun) =====
contrast_data = (fre_gowin + fre_nogopun) - (fre_nogowin + fre_gopun);
data.powspctrm = contrast_data ./ 2;   % Divide by 2 to keep the same magnitude

% ===== Midfrontal average =====
dat2plot = squeeze(mean(data.powspctrm(MFchans,:,:), 1));

% ===== Draw TF map =====
contourf(freq.time, freq.freq, dat2plot, 40, 'linestyle','none')
set(gca,'ytick',[2 4 8 16 32],'yscale','log','xlim',[-0.1 1],...
    'xtick',[0 1],'xtickLabel',{'FB','1'},'ylim',[3 50],...
    'fontsize',11,'clim',[-1 1])

% ===== Mark ROI =====
rectangle('Position',[MFtime(1), MFfreq(1), diff(MFtime), diff(MFfreq)],...
          'EdgeColor','k','LineWidth',1.5)

colorbar
title('Biased vs. Unbiased learning','fontsize',12)
ylabel('Frequency (Hz)','fontsize',11,'fontweight','bold')
xlabel('Time (s)','fontsize',11,'fontweight','bold')


%% Biased vs. unbiased learning contrast topoplot

% Get electrode information
elec = freq.elec;

% Define midfrontal channels
MFchans = [6 11 127];

figure('name','midfrontal TF power Valence x Required action','Position',[100 500 500 200]) % Increase width

% Use the same colormap and color range as above
cmap = jet(64);
colormap(cmap);

% Draw topoplot with ft_topoplotTFR
cfg = []; 
cfg.ylim   = MFfreq; 
cfg.zlim   = [-1.5 1.5];      % Keep consistent with the above section
cfg.marker = 'on'; 
cfg.style  = 'straight'; 

% Manually adjust electrode coordinates
elec = freq.elec;                     % Get electrode information
elec.chanpos(:, 1) = -freq.elec.chanpos(:, 2);  % Reverse x-coordinate
elec.chanpos(:, 2) =  freq.elec.chanpos(:, 1);  % Set y-coordinate to original x-coordinate
cfg.elec  = elec;                      % Update electrode layout

cfg.comment          = 'no'; 
cfg.xlim             = MFtime;
cfg.highlight        = 'on'; 
cfg.highlightchannel = MFchans; 
cfg.highlightsymbol  = 'o';

cfg.colormap = cmap;   % Make ft_topoplotTFR use the same colormap

% Draw topo
ft_topoplotTFR(cfg, data);

% Force the same color range and show colorbar
caxis([-1.5 1.5]);      % Equivalent to set(gca,'clim',[-1.5 1.5]);
colorbar;
set(gca,'fontsize',11);



%% One-sample t-test for biased vs. unbiased learning contrast

gopun = zeros(36, 1);  % Store power for gopun
nogopun = zeros(36, 1);  % Store power for nogopun
gowin = zeros(36, 1);  % Store power for gowin
nogowin = zeros(36, 1);  % Store power for nogowin

for iSub = 1:36   
    gopun(iSub) = squeeze(mean(mean(mean(gopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    gowin(iSub) = squeeze(mean(mean(mean(gowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));    
end

% ===== Build contrast: (GoWin + GoPun) - (NoGoWin + NoGoPun) =====
contrast_eachSub = (gowin + nogopun) - (nogowin + gopun);

% ===== One-sample t-test (test whether the contrast differs from 0) =====
[~, p] = ttest(contrast_eachSub);

% ===== Compute mean and standard error (for plotting) =====
means = mean(contrast_eachSub);
sems  = std(contrast_eachSub) / sqrt(36);

% ===== Draw bar plot =====
figure;
b = bar(means, 'FaceColor', 'flat'); 

% Use a soft color
b.CData = [0.1 0.5 0.9];  % Lake blue

ylim([-2 2.5]);  
alpha(0.85);
hold on;

% Add error bar
errorbar(1, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 8);

% Add p-value label and horizontal line
y_line = 2;  
y_pval = 2.1; 
plot([0.8, 1.2], [y_line, y_line], 'k', 'LineWidth', 1.5);  
text(1, y_pval, sprintf('p = %.3f', p), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');  

% Beautify axes
xticks(1);
xticklabels({'Contrast (Biased - Unbiased)'});
ylabel('Theta Power (dB)', 'FontWeight', 'bold', 'FontSize', 16);
title('Relative Biased vs. Unbiased Learning (Single-Sample T-Test)', ...
      'FontWeight', 'bold', 'FontSize', 16);
set(gca,'FontSize',14,'FontWeight','bold','XColor','k','YColor','k');
grid off;
hold off;


%% Compute t-test between biased vs. unbiased learning contrasts and draw bar plot
gopun = zeros(36, 1);  
nogopun = zeros(36, 1);  
gowin = zeros(36, 1);  
nogowin = zeros(36, 1);  

% Compute mean power for each subject
for iSub = 1:36   
    gopun(iSub) = squeeze(mean(mean(mean(gopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogopun(iSub) = squeeze(mean(mean(mean(nogopun_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    gowin(iSub) = squeeze(mean(mean(mean(gowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));
    nogowin(iSub) = squeeze(mean(mean(mean(nogowin_sub(iSub,MFchans,iFreq,iTime), 2), 3), 4));    
end

% Compute and combine data
valxaccAction1 = [nogopun, gowin];  % First contrast group
valxaccAction2 = [gopun, nogowin];  % Second contrast group

% Compute means and standard errors (SE)
means = [mean(mean(valxaccAction1)), mean(mean(valxaccAction2))]; % Group means
sems = [std(mean(valxaccAction1))/sqrt(37), std(mean(valxaccAction2))/sqrt(36)]; % Standard errors

% Run t-test
[~, p] = ttest(mean(valxaccAction1 - valxaccAction2));

% Draw bar plot
figure;
b = bar(means, 'FaceColor', 'flat'); 

% Use softer gradient colors
b.CData = [0.1 0.5 0.9;  % Lake blue (Biased Learning)
           0.9 0.4 0.3]; % Coral red (Unbiased Learning)

% Fix Y-axis range
ylim([0 2.5]);  

% Add transparency for softer colors
alpha(0.85);

% Add error bars
hold on;
errorbar(1:2, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 8);

% Add p-value label and horizontal line
y_line = 2;
y_pval = 2.1;

plot([1, 2], [y_line, y_line], 'k', 'LineWidth', 1.5);
text(1.5, y_pval, sprintf('p = %.3f', p), 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Beautify x-axis labels
xticks([0.8, 2.2]); % Increase spacing between bars
xticklabels({'Biased Learning', 'Unbiased Learning'});
ylabel('Theta Power (dB)', 'FontWeight', 'bold', 'FontSize', 16);
title('Comparison of Learning Bias Effects', 'FontWeight', 'bold', 'FontSize', 16);

% Remove gridlines for a cleaner look
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'XColor', 'k', 'YColor', 'k');
grid off;

hold off;