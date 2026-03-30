%% cue-locked
%% step 1: source estimation: ERP-based source activity (to ft)
clear all; close all; clc
dbstop if error

% Add FieldTrip path
addpath('D:\MATLAB    oolbox\fieldtrip-20211209\fieldtrip-20211209'); % Make sure the path is correct
ft_defaults;

% Set data path
path = 'G:\MT\reference\rename\cue_epoch\out_put_ICA\cue\'; % main path
tic;

% Set subject range
range = 1:37;  % Subject IDs from 02 to 16

%% TF parameters
parTF.minFreq         = 1;
parTF.maxFreq         = 50;
parTF.nFreq           = 40;
parTF.freq4tf         = logspace(log10(parTF.minFreq),log10(parTF.maxFreq),parTF.nFreq);
parTF.width4tf        = 4;
parTF.toi4tf          = -.25:.025:3;
parTF.baselinetime    = [-.25 -.05];  % Define baseline time window

% Loop over subjects
for ns = range
    filename = sprintf('HC%02dgtw.set', ns);  % Build filename
    inpath = path;  % Input path

    cd(inpath);  % Change to data directory

    % Load data
    EEG = pop_loadset('filename', filename, 'filepath', inpath);  % Load data with EEGLAB

    % Convert EEG data to FieldTrip structure
    cfg = [];
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');  % Convert to FieldTrip format

    % Compute scalp current density
    cfg                     = [];
    cfg.method              = 'spline';  % Or use another method
    cfg.elec                = ft_data.elec;  % Assume electrode positions are defined
    scd                     = ft_scalpcurrentdensity(cfg, ft_data);  % Compute CSD

%     %% perform wavelet convolution
    cfg = [];
    cfg.method          = 'wavelet';
    cfg.output          = 'pow'; % 'fourier'; % to also keep the phase.
    cfg.keeptrials      = 'yes';
    cfg.foi             = parTF.freq4tf;      % frequencies of interest
    cfg.width           = parTF.width4tf;     % number of cycles
    cfg.toi             = parTF.toi4tf;       % times (s) to center the analysis window
    
    % Run frequency analysis
     freq                = ft_freqanalysis(cfg, scd);
    
%     % Baseline correction
%     cfg_baseline              = [];
%     cfg_baseline.baseline     = parTF.baselinetime;  % Set baseline time window
%     cfg_baseline.baselinetype = 'db';       
%     freq                      = ft_freqbaseline(cfg_baseline, freq);

    clear scd  % Clear intermediate variable

    % Create output path
    outpath1 = fullfile(path, 'ft_lap');  % Create output folder path
    if ~exist(outpath1, 'dir')  % Check whether the folder exists
        mkdir(outpath1);  % Create directory
    end
    
    % Save processed data
    save_filename = fullfile(outpath1, [filename(1:end-4) '.mat']);
    save(save_filename, 'freq', 'ft_data', '-v7.3');  % Save as .mat file
end



%%
%% response-locked
%% step 1: source estimation: ERP-based source activity (to ft)
clear all; close all; clc
dbstop if error

% Add FieldTrip path
addpath('D:\MATLAB    oolbox\fieldtrip-20211209\fieldtrip-20211209'); % Make sure the path is correct
ft_defaults;

% Set data path
path = 'G:\MT\reference\rename\cue_epoch\out_put_ICA\resp\'; % main path
tic;

 % Set subject range
 range = 2:2;  % Subject IDs from 02 to 16

%% TF parameters
parTF.minFreq         = 1;
parTF.maxFreq         = 50;
parTF.nFreq           = 40;
parTF.freq4tf         = logspace(log10(parTF.minFreq),log10(parTF.maxFreq),parTF.nFreq);
parTF.width4tf        = 4;
parTF.toi4tf          = -1:.025:1;
parTF.baselinetime    = [-.25 -.05];  % Define baseline time window

% Loop over subjects
for ns = range
    filename = sprintf('HC%02d.set', ns);  % Build filename
    inpath = path;  % Input path

    cd(inpath);  % Change to data directory

    % Load data
    EEG = pop_loadset('filename', filename, 'filepath', inpath);  % Load data with EEGLAB

    % Convert EEG data to FieldTrip structure
    cfg = [];
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');  % Convert to FieldTrip format

    % Compute scalp current density
    cfg                     = [];
    cfg.method              = 'spline';  % Or use another method
    cfg.elec                = ft_data.elec;  % Assume electrode positions are defined
    scd                     = ft_scalpcurrentdensity(cfg, ft_data);  % Compute CSD

%     %% perform wavelet convolution
    cfg = [];
    cfg.method          = 'wavelet';
    cfg.output          = 'pow'; % 'fourier'; % to also keep the phase.
    cfg.keeptrials      = 'yes';
    cfg.foi             = parTF.freq4tf;      % frequencies of interest
    cfg.width           = parTF.width4tf;     % number of cycles
    cfg.toi             = parTF.toi4tf;       % times (s) to center the analysis window
    
    % Run frequency analysis
     freq                = ft_freqanalysis(cfg, scd);
    
%     % Baseline correction
%     cfg_baseline              = [];
%     cfg_baseline.baseline     = parTF.baselinetime;  % Set baseline time window
%     cfg_baseline.baselinetype = 'db';       
%     freq                      = ft_freqbaseline(cfg_baseline, freq);

    clear scd  % Clear intermediate variable

    % Create output path
    outpath1 = fullfile(path, 'ft_lap');  % Create output folder path
    if ~exist(outpath1, 'dir')  % Check whether the folder exists
        mkdir(outpath1);  % Create directory
    end
    
    % Save processed data
    save_filename = fullfile(outpath1, [filename(1:end-4) '.mat']);
    save(save_filename, 'freq', 'ft_data', '-v7.3');  % Save as .mat file
end

% Set subject range
range = 1:36;  % Subject IDs from 02 to 16
% Loop over subjects
for ns = range
    filename = sprintf('OCD%02dgtw.set', ns);  % Build filename
    inpath = path;  % Input path

    cd(inpath);  % Change to data directory

    % Load data
    EEG = pop_loadset('filename', filename, 'filepath', inpath);  % Load data with EEGLAB

    % Convert EEG data to FieldTrip structure
    cfg = [];
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');  % Convert to FieldTrip format

    % Compute scalp current density
    cfg                     = [];
    cfg.method              = 'spline';  % Or use another method
    cfg.elec                = ft_data.elec;  % Assume electrode positions are defined
    scd                     = ft_scalpcurrentdensity(cfg, ft_data);  % Compute CSD

%     %% perform wavelet convolution
    cfg = [];
    cfg.method          = 'wavelet';
    cfg.output          = 'pow'; % 'fourier'; % to also keep the phase.
    cfg.keeptrials      = 'yes';
    cfg.foi             = parTF.freq4tf;      % frequencies of interest
    cfg.width           = parTF.width4tf;     % number of cycles
    cfg.toi             = parTF.toi4tf;       % times (s) to center the analysis window
    
    % Run frequency analysis
     freq                = ft_freqanalysis(cfg, scd);
    
%     % Baseline correction
%     cfg_baseline              = [];
%     cfg_baseline.baseline     = parTF.baselinetime;  % Set baseline time window
%     cfg_baseline.baselinetype = 'db';       
%     freq                      = ft_freqbaseline(cfg_baseline, freq);

    clear scd  % Clear intermediate variable

    % Create output path
    outpath1 = fullfile(path, 'ft_lap');  % Create output folder path
    if ~exist(outpath1, 'dir')  % Check whether the folder exists
        mkdir(outpath1);  % Create directory
    end
    
    % Save processed data
    save_filename = fullfile(outpath1, [filename(1:end-4) '.mat']);
    save(save_filename, 'freq', 'ft_data', '-v7.3');  % Save as .mat file
end




%% baseline
clc;
clear;

% Compute baseline, apply baseline correction, and save corrected data for each subject
baseline_dir = 'G:\MT\reference\rename\cue_epoch\out_put_ICA\cue';  % Baseline data folder
golefttowin_dir = 'G:\MT\reference\rename\cue_epoch\out_put_ICA\cue\ft_lap';  % Original data folder
output_dir = 'G:\MT\reference\rename\cue_epoch\out_put_ICA\cue\ft_lap\baseline';  % Output folder

% Create output folder if it does not exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Get all .mat files in baseline and original data folders
baseline_files = dir(fullfile(baseline_dir, '*.mat'));  % Get all baseline data files
golefttowin_files = dir(fullfile(golefttowin_dir, '*.mat'));  % Get all original data files

% Loop through each file in the baseline folder
for i = 1:length(baseline_files)
    % Get file names
    baseline_filename = fullfile(baseline_dir, baseline_files(i).name);  % Baseline data file path
    golefttowin_filename = fullfile(golefttowin_dir, golefttowin_files(i).name);  % Original data file path
    
    % Load baseline data
    load(baseline_filename);  % Assume the file contains freq.powspctrm data
    baseline = mean_baseline;
    
    % Load original data (golefttowin)
    load(golefttowin_filename);  % Assume the file contains freq.powspctrm data
    oridata = permute(freq.powspctrm, [2, 3, 4, 1]);  % Adjust dimension order
    oridata_trial = squeeze(mean(oridata(:, :, :, :), 4));
    
    % Check the dimensions of oridata and baseline
    disp(size(oridata_trial));  % Show size of oridata
    disp(size(baseline));  % Show size of baseline
    
    % Expand baseline according to the third dimension of oridata
    baselineExpanded = repmat(baseline, [1, 1, size(oridata_trial, 3)]);
    
    % Check the size of baselineExpanded
    disp(size(baselineExpanded));  % Show size of baselineExpanded
    
    % Compute baseline-corrected data
    basefreq.powspctrm = 10 * log10(oridata_trial ./ baselineExpanded);  % Apply baseline correction
    
    % Save results
    output_filename = fullfile(output_dir, sprintf('%s_base_correct.mat', baseline_files(i).name(1:end-4)));
    save(output_filename, 'basefreq');  % Save basefreq.powspctrm
end

disp('All data processing completed');


