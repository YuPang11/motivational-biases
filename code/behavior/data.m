clc; clear;

% 1) Read data
filePath = "C:\Users\16055\behavior\oridata.xlsx";
T = readtable(filePath);

% 2) Get all unique subject IDs
subject_list = unique(T.subject);
nSub = numel(subject_list);

% 3) Preallocate structure
data = struct();
data.go      = [];
data.acc     = [];
data.stim    = [];
data.RT      = [];
data.outcom  = [];
data.resp    = [];

% 4) Loop over subjects and build matrices
for i = 1:nSub
    subj_id = subject_list(i);
    Tsub = T(T.subject == subj_id, :);

    % Trial-wise data for each subject
    data.go(i,:)     = Tsub.GO';
    data.acc(i,:)    = Tsub.ACC';
    data.stim(i,:)   = Tsub.stimulus';
    data.RT(i,:)     = Tsub.RT';
    data.outcom(i,:) = Tsub.outcom';
    data.resp(i,:)   = Tsub.RESP';
end

% 5) Save parameter info
par = struct();
par.nSub = nSub;
par.subjList = subject_list;

% 6) Save as one combined file
save_dir = "H:\mt_model\OCD";
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save(fullfile(save_dir, "EEGsingletrial.mat"), 'data', 'par');

disp("Done: EEGsingletrial.mat has been created.");





%%
%% Includes group labels and exports data for behavioral linear models
clc; clear;

mat_path = "H:\mt_behavalluse\EEGsingletrial.mat";
load(mat_path, 'data', 'par');  % -> data.go/acc/stim/RT/outcom/resp, par.nSub, par.subjList

%% ===== Define group lists =====
HC_list  = [102 105 106 107 108 109 110 113 115 116 117 118 124 125 130 132 133 134 ...
             136 137 138 139 140 141 142 143 145 146 147 148 150 151 152 154 155 156 157];
OCD_list = [201 207 209 210 211 213 214 216 217 218 219 220 221 224 225 227 231 232 ...
             233 234 238 241 243 245 247 249 251 253 254 255 256 258 259 260 261 262];

%% ===== Map stim to cue valence and required action =====
valence_map = containers.Map('KeyType','double','ValueType','double');
action_map  = containers.Map('KeyType','double','ValueType','double');
for s = 1:8
    if ismember(s,[1 2 5 6])  % win cues
        valence_map(s) = 1;
    else                      % avoid cues (3,4,7,8)
        valence_map(s) = 0;
    end
    if ismember(s,[1 2 3 4])  % go cues
        action_map(s) = 1;
    else                      % nogo cues (5,6,7,8)
        action_map(s) = 0;
    end
end

%% ===== Loop over all subjects =====
rows = {};  % collect all trial rows
for iSub = 1:par.nSub
    sID = par.subjList(iSub);
    nTrial = size(data.go, 2);

    stim_i   = data.stim(iSub, :);
    DVgo_i   = data.go(iSub, :);
    DVacc_i  = data.acc(iSub, :);
    RTms_i   = data.RT(iSub, :);
    outcom_i = data.outcom(iSub, :);
    resp_i   = data.resp(iSub, :);

    % Determine group
    if ismember(sID, OCD_list)
        group_i = 1;  % OCD
    elseif ismember(sID, HC_list)
        group_i = 0;  % HC
    else
        warning('Subject %d is not classified as HC or OCD. Set to NaN.', sID);
        group_i = NaN;
    end

    % Convert RT from ms to seconds
    DVrt_i = RTms_i ./ 1000;

    % Derive valence and action from stim
    val_i = arrayfun(@(x) valence_map(x), stim_i);
    act_i = arrayfun(@(x) action_map(x),  stim_i);

    % Organize as column vectors plus group
    rows{end+1} = table( ...
        repmat(sID, [nTrial,1]), ...
        repmat(group_i, [nTrial,1]), ...
        DVgo_i(:), ...
        DVacc_i(:), ...
        DVrt_i(:), ...
        val_i(:), ...
        act_i(:), ...
        stim_i(:), ...
        outcom_i(:), ...
        resp_i(:), ...
        'VariableNames', {'sID','group','DVgo','DVacc','DVrt','valence','action','stim','outcom','resp'});
end

%% ===== Concatenate all trials and export =====
EEG4mixedmodelR = vertcat(rows{:});

out_dir = "H:\mt_behavalluse";
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

out_file = fullfile(out_dir, "EEG4mixedmodelR_withGroup.csv");
writetable(EEG4mixedmodelR, out_file);

disp("Done: CSV has been created: " + out_file + " (group included: HC=0, OCD=1)");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavior + EEG model data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; dbstop if error;

% 1) Read grouped data, OCD example
filePath = "F:\PIT\model\OCD3.xlsx";
T = readtable(filePath);

% 2) Get subject list and trial count
subject_list = unique(T.subject, 'stable');   % keep original order
nSub = numel(subject_list);

if ismember('trial', T.Properties.VariableNames)
    nTrial = max(T.trial);                    % assume trial starts at 1 and is continuous
else
    error('Excel is missing the "trial" column.');
end

% 3) Preallocate structure (rows = subjects, cols = trials)
data = struct();

% Behavioral variables
data.go     = nan(nSub, nTrial);
data.acc    = nan(nSub, nTrial);
data.stim   = nan(nSub, nTrial);
data.RT     = nan(nSub, nTrial);
data.outcom = nan(nSub, nTrial);
data.resp   = nan(nSub, nTrial);

% EEG and indicator variables
data.MFpower      = nan(nSub, nTrial);
data.ICPCpfc      = nan(nSub, nTrial);
data.l_ICPCmotor  = nan(nSub, nTrial);
data.r_ICPCmotor  = nan(nSub, nTrial);
data.IVleft       = nan(nSub, nTrial);
data.IVright      = nan(nSub, nTrial);
data.trial2use    = zeros(nSub, nTrial);  % 0 = missing/rejected, 1 = usable

% 4) Loop over subjects and fill matrices by trial index
for i = 1:nSub
    subj_id = subject_list(i);
    Tsub = T(T.subject == subj_id, :);

    % Sort by trial to keep column order consistent
    Tsub = sortrows(Tsub, 'trial');

    % Trial index should fall within 1..nTrial
    idx = Tsub.trial(:)';
    if any(idx < 1 | idx > nTrial)
        error('Subject %g has out-of-range trial indices.', subj_id);
    end

    % Behavioral variables
    data.go(i,idx)     = Tsub.GO';
    data.acc(i,idx)    = Tsub.ACC';
    data.stim(i,idx)   = Tsub.stimulus';
    data.RT(i,idx)     = Tsub.RT';
    data.outcom(i,idx) = Tsub.outcom';
    data.resp(i,idx)   = Tsub.RESP';

    % EEG and indicator variables
    data.MFpower(i,idx)     = Tsub.MFpower';
    data.ICPCpfc(i,idx)     = Tsub.ICPCpfc';
    data.l_ICPCmotor(i,idx) = Tsub.l_ICPCmotor';
    data.r_ICPCmotor(i,idx) = Tsub.r_ICPCmotor';
    data.IVleft(i,idx)      = Tsub.IVleft';
    data.IVright(i,idx)     = Tsub.IVright';
    data.trial2use(i,idx)   = Tsub.trial2use';
end

% 5) Save parameter info
par = struct();
par.nSub     = nSub;
par.nTrial   = nTrial;
par.subjList = subject_list;
if ismember('stimulus', T.Properties.VariableNames)
    par.nStim = numel(unique(T.stimulus));
else
    par.nStim = [];
end
par.note = 'trial2use: 1 = usable EEG, 0 = artifact/missing; missing EEG values remain NaN.';

% 6) Save as a .mat file compatible with the original structure
save_dir = "F:\PIT\model\EEGmodelHC";
if ~exist(save_dir, 'dir'); mkdir(save_dir); end
save(fullfile(save_dir, "EEGsingletrialocd3.mat"), 'data', 'par');

disp("Done: EEGsingletrialocd3.mat has been created (behavior + EEG variables included).");





%% Stan input data with EEG metrics: strict original-name version
clc; clear; dbstop if error;

% Paths
dirs.singletrial = 'F:\PIT\model\EEGmodelHC\EEGsingletrialocd3.mat';
dirs.stan        = 'F:\PIT\model\EEGmodelHC\MFpower4Rocd3.mat';

% Generate Stan input data
EEGpav_prepStan4Reeg(dirs);

%% Quick self-check
fprintf('\n=== File contents ===\n');
whos -file(dirs.stan)

load(dirs.stan)  % load ya/a/s/r/rew/N/Nsub and EEG metrics if present

% Small helper functions
printMatInfo = @(name, M) fprintf('%-12s  size=%s  NaN count=%d  NaN ratio=%.6f\n', ...
    name, mat2str(size(M)), sum(~isfinite(M(:))), mean(~isfinite(M(:))));
checkSize = @(name,M,refsz) assert(isequal(size(M),refsz), ...
    '%s size does not match reference [%d %d] (current %s)', name, refsz(1), refsz(2), mat2str(size(M)));

fprintf('\n=== Basic behavioral checks ===\n');

if exist('N','var') && exist('Nsub','var')
    fprintf('Nsub=%d, N=%d\n', Nsub, N);
else
    error('Missing N or Nsub.');
end

if exist('s','var')
    u_s = unique(s(:))';      fprintf('unique(s): ');  disp(u_s);
    if ~all(ismember(u_s,1:8)), warning('s contains values outside 1..8.'); end
    checkSize('s', s, [Nsub, N]);
end

if exist('rew','var')
    u_rew = unique(rew(:))';  fprintf('unique(rew): ');  disp(u_rew);  % expected [0 1]
    if ~all(ismember(u_rew,[0 1])), warning('rew is not binary (0/1).'); end
    checkSize('rew', rew, [Nsub, N]);
end

if exist('r','var')
    u_r = unique(r(:))';      fprintf('unique(r): ');    disp(u_r);    % expected [-1 0 1]
    if ~all(ismember(u_r,[-1 0 1])), warning('r is not in {-1,0,1}.'); end
    checkSize('r', r, [Nsub, N]);
end

if exist('ya','var')
    fprintf('size(ya): ');  disp(size(ya));               % expected [Nsub N]
    u_ya = unique(ya(~isnan(ya)))';
    fprintf('unique(ya(~isnan)): '); disp(u_ya);          % expected [0 1 2]
    if ~all(ismember(u_ya,[0 1 2])), warning('ya contains values outside {0,1,2}.'); end
    checkSize('ya', ya, [Nsub, N]);
end

if exist('a','var')
    u_a = unique(a(~isnan(a)))';
    fprintf('unique(a(~isnan)): '); disp(u_a);            % expected [1 2 3]
    if ~all(ismember(u_a,[1 2 3])), warning('a contains values outside {1,2,3}.'); end
    checkSize('a', a, [Nsub, N]);
end

% sub / subIDs
if exist('sub','var')
    fprintf('sub (first 10): ');  disp(sub(1:min(10,numel(sub)))');
    if numel(sub) ~= Nsub, warning('sub length (%d) ~= Nsub (%d).', numel(sub), Nsub); end
end
if exist('subIDs','var')
    fprintf('subIDs (first 10): ');  disp(subIDs(1:min(10,numel(subIDs)))');
    if numel(subIDs) ~= Nsub, warning('subIDs length (%d) ~= Nsub (%d).', numel(subIDs), Nsub); end
end

%% EEG metric checks
fprintf('\n=== EEG metric checks (if present) ===\n');
hasTrial2Use = exist('trial2use','var') == 1;

if exist('MFpower','var'),        printMatInfo('MFpower', MFpower);             checkSize('MFpower', MFpower, [Nsub, N]); end
if exist('ICPCpfc','var'),        printMatInfo('ICPCpfc', ICPCpfc);             checkSize('ICPCpfc', ICPCpfc, [Nsub, N]); end
if exist('l_ICPCmotor','var'),    printMatInfo('l_ICPCmotor', l_ICPCmotor);     checkSize('l_ICPCmotor', l_ICPCmotor, [Nsub, N]); end
if exist('r_ICPCmotor','var'),    printMatInfo('r_ICPCmotor', r_ICPCmotor);     checkSize('r_ICPCmotor', r_ICPCmotor, [Nsub, N]); end

% trial2use
if hasTrial2Use
    trial2use = double(trial2use);
    trial2use(~isfinite(trial2use)) = 0;
    trial2use(trial2use ~= 0) = 1;
    tot = numel(trial2use); ok = nnz(trial2use == 1);
    fprintf('trial2use: total=%d, usable=%d (%.3f)\n', tot, ok, ok/max(1,tot));

    if exist('MFpower','var')
        fprintf('  MFpower NaN count on usable trials = %d\n', sum(isnan(MFpower(trial2use==1))));
    end
    if exist('ICPCpfc','var')
        fprintf('  ICPCpfc NaN count on usable trials = %d\n', sum(isnan(ICPCpfc(trial2use==1))));
    end
    if exist('l_ICPCmotor','var')
        fprintf('  l_ICPCmotor NaN count on usable trials = %d\n', sum(isnan(l_ICPCmotor(trial2use==1))));
    end
    if exist('r_ICPCmotor','var')
        fprintf('  r_ICPCmotor NaN count on usable trials = %d\n', sum(isnan(r_ICPCmotor(trial2use==1))));
    end
end

% IVleft / IVright
if exist('IVleft','var')
    u = unique(IVleft(:))'; fprintf('unique(IVleft): '); disp(u);
    fprintf('sum(IVleft==1)=%d\n', nnz(IVleft==1));
    checkSize('IVleft', IVleft, [Nsub, N]);
end
if exist('IVright','var')
    u = unique(IVright(:))'; fprintf('unique(IVright): '); disp(u);
    fprintf('sum(IVright==1)=%d\n', nnz(IVright==1));
    checkSize('IVright', IVright, [Nsub, N]);
end

%% Size consistency check against behavioral matrices
if exist('ya','var')
    if exist('MFpower','var'),     assert(all(size(MFpower)     == size(ya)), 'MFpower size does not match ya'); end
    if exist('ICPCpfc','var'),     assert(all(size(ICPCpfc)     == size(ya)), 'ICPCpfc size does not match ya'); end
    if exist('l_ICPCmotor','var'), assert(all(size(l_ICPCmotor) == size(ya)), 'l_ICPCmotor size does not match ya'); end
    if exist('r_ICPCmotor','var'), assert(all(size(r_ICPCmotor) == size(ya)), 'r_ICPCmotor size does not match ya'); end
    if exist('IVleft','var'),      assert(all(size(IVleft)      == size(ya)), 'IVleft size does not match ya'); end
    if exist('IVright','var'),     assert(all(size(IVright)     == size(ya)), 'IVright size does not match ya'); end
    if exist('trial2use','var'),   assert(all(size(trial2use)   == size(ya)), 'trial2use size does not match ya'); end
end

fprintf('\nDone: self-check completed.\n');