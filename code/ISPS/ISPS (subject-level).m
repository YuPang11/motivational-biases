%% ISPC calculation - batch for HC + OCD (all subjects found in folders)
clear; clc; close all;
dbstop if error;

% Add FieldTrip path
addpath('D:\MATLAB\toolbox\fieldtrip-20211209\fieldtrip-20211209');
ft_defaults;

% ===================== Parameter settings =====================
parTF.minFreq = 1;
parTF.maxFreq = 50;
parTF.nFreq   = 40;
parTF.freq4tf = logspace(log10(parTF.minFreq), log10(parTF.maxFreq), parTF.nFreq);
parTF.freq4tf(parTF.freq4tf < 4) = [];   % Keep only 4–8 Hz
parTF.freq4tf(parTF.freq4tf > 8) = [];
parTF.nFreq   = numel(parTF.freq4tf);
parTF.width4tf = 4;
parTF.toi4tf = -0.25:0.025:3;            % Time window
parTF.baselinetime = [-0.25 -0.05];
parTF.MFchan = [6 127];                  % Midfrontal seeds

nPerm   = 100;    % Number of permutations
nCutoff = 20;     % Minimum trial count

% Folder paths for the 6 conditions
condNames = {'GLW','GRW','GLA','GRA','NGW','NGA'};% Condition order
baseDir   = 'H:\MT_all(afterqujizhi)\cueoriandbu\correct\TOP\HC\CONGRUENT';

% ===================== Get all subject IDs =====================
exampleDir = fullfile(baseDir, condNames{1});
fileList   = dir(fullfile(exampleDir,'*.set'));

subIDs = {};
for i = 1:length(fileList)
    fname = fileList(i).name;  % e.g. "HC02.set" or "OCD15.set"
    subIDs{i} = regexprep(fname, '\.set$', '');
end
subIDs = unique(subIDs); % Remove duplicates and sort

fprintf('Total subjects found: %d (HC + OCD)\n', numel(subIDs));

% ===================== Main loop - subjects =====================
for s = 1:numel(subIDs)
    subID = subIDs{s};
    fprintf('=== Processing %s ===\n', subID);

    ICPCall = cell(6,1);  % Store results for 6 conditions
    nTrial  = NaN(6,1);

    % ---- Step 1: first count trials in each condition ----
    for iCondi = 1:6
        condDir = fullfile(baseDir, condNames{iCondi});
        filename = fullfile(condDir, sprintf('%s.set', subID));

        if ~exist(filename,'file')
            continue;
        end

        EEG = pop_loadset('filename', sprintf('%s.set', subID), 'filepath', condDir);
        nTrial(iCondi) = EEG.trials;
    end

    % Mark conditions below cutoff as NaN
    nTrial(nTrial < nCutoff) = NaN;

    % If all 6 conditions are below cutoff -> skip this subject
    if all(isnan(nTrial))
        fprintf('Skipped %s (all conditions < %d trials)\n', subID, nCutoff);
        continue;
    end

    % Find the minimum valid trial count for this subject
    minTrialN = min(nTrial, [], 'omitnan');
    fprintf('%s minTrialN = %d (after cutoff)\n', subID, minTrialN);

    % ---- Step 2: loop through 6 conditions and sample to minTrialN ----
    for iCondi = 1:6
        condDir = fullfile(baseDir, condNames{iCondi});
        filename = fullfile(condDir, sprintf('%s.set', subID));

        if ~exist(filename,'file') || isnan(nTrial(iCondi))
            continue;  % Skip missing conditions or those below cutoff
        end

        % Load data
        EEG = pop_loadset('filename', sprintf('%s.set', subID), 'filepath', condDir);

        % Convert to FieldTrip structure
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing');

        % Laplacian
        cfg = [];
        cfg.method = 'spline';
        cfg.elec   = ft_data.elec;
        scd        = ft_scalpcurrentdensity(cfg, ft_data);

        % Time-frequency analysis
        cfg = [];
        cfg.method     = 'wavelet';
        cfg.output     = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.foi        = parTF.freq4tf;
        cfg.width      = parTF.width4tf;
        cfg.toi        = parTF.toi4tf;
        freq = ft_freqanalysis(cfg, scd);

        % Phase extraction
        seedPhase   = angle(freq.fourierspctrm(:, parTF.MFchan, :, :));
        targetPhase = angle(freq.fourierspctrm);
        parTF.nChan = size(EEG.chanlocs, 2);

        % ----------------- Phase synchrony calculation -----------------
        ICPC = NaN(numel(parTF.MFchan), parTF.nChan, parTF.nFreq, length(parTF.toi4tf), nPerm);

        for iPerm = 1:nPerm
            permidx = randperm(nTrial(iCondi));
            permidx = permidx(1:minTrialN);   % Key: resample to minTrialN for all conditions

            for iSeed = 1:numel(parTF.MFchan)
                for iTarget = 1:parTF.nChan
                    if parTF.MFchan(iSeed) == iTarget; continue; end
                    phasediff = seedPhase(permidx, iSeed, :, :) - targetPhase(permidx, iTarget, :, :);
                    ICPC(iSeed, iTarget, :, :, iPerm) = abs(mean(exp(1i * phasediff),1));
                end
            end
        end

        % Average across permutations
        ICPCall{iCondi} = mean(ICPC, 5);
        clear EEG ft_data scd freq ICPC
    end

    % ---- Save subject result ----
    saveFile = fullfile(baseDir, sprintf('%s_ISPC.mat', subID));
    save(saveFile, 'ICPCall','nTrial','minTrialN','parTF');
    fprintf('Saved results: %s\n', saveFile);
end
