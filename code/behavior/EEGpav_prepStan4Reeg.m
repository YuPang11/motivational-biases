function EEGpav_prepStan4Reeg(dirs)
% Generate EEG4stan.mat
% Input:
%   dirs.singletrial = '.../EEGsingletrial.mat'
%   dirs.stan        = '.../EEG4stan.mat'

    % ---- Load single-trial data ----
    S = load(dirs.singletrial, 'data', 'par');
    data = S.data;
    par  = S.par;

    % ---- Feedback r: support outcome/outcom ----
    if isfield(data, 'outcome')
        r = double(data.outcome);
    elseif isfield(data, 'outcom')
        r = double(data.outcom);
    else
        error('data.outcome or data.outcom was not found.');
    end

    % ---- Observed response ya: expected {0,1,2} = {NoGo, Left, Right} ----
    % ---- Unify key codes to 1/2 ----
    ya = double(data.resp);
    ya(ya==65 | ya==97)  = 1;  % A/a -> left
    ya(ya==69 | ya==101) = 2;  % E/e -> right
    valid = isnan(ya) | ya==0 | ya==1 | ya==2;
    ya(~valid) = NaN;

    % ---- a: recode NoGo as 3 for Stan indexing ----
    a = ya;
    a(ya==0) = 3;

    % ---- s: stimulus ID (expected 1..8) ----
    s = double(data.stim);
    if any(s(:) < 1 | s(:) > 8)
        warning('Stimulus IDs outside 1..8 were detected; values will be remapped to 1..8 using mod8.');
        s = mod(s-1, 8) + 1;
    end

    % ---- rew: cue valence class (Win=1 / Avoid=0) ----
    % ---- based on {1,2,5,6} / {3,4,7,8} ----
    rew = double(ismember(s, [1 2 5 6]));  % Win=1, otherwise 0

    % ---- N / Nsub ----
    [Nsub, N] = size(ya);

    % ---- sub / subIDs: use par.subjList if available ----
    % ---- otherwise use 1..Nsub ----
    if isfield(par,'subjList') && ~isempty(par.subjList)
        sub    = double(par.subjList(:));
        subIDs = sub;  % duplicate for compatibility with later R scripts
        if numel(sub) ~= Nsub
            warning('Length of sub (%d) does not match Nsub (%d); using 1..Nsub instead.', numel(sub), Nsub);
            sub    = (1:Nsub).';
            subIDs = sub;
        end
    else
        warning('par.subjList is missing or empty; using 1..Nsub as sub / subIDs.');
        sub    = (1:Nsub).';
        subIDs = sub;
    end

    % ---- EEG metrics: load if available ----
    % ---- keep the same size as behavioral matrices ----
    if isfield(data,'MFpower'),      MFpower      = double(data.MFpower);      end
    if isfield(data,'ICPCpfc'),      ICPCpfc      = double(data.ICPCpfc);      end
    if isfield(data,'l_ICPCmotor'),  l_ICPCmotor  = double(data.l_ICPCmotor);  end
    if isfield(data,'r_ICPCmotor'),  r_ICPCmotor  = double(data.r_ICPCmotor);  end
    if isfield(data,'IVleft'),       IVleft       = double(data.IVleft);       end
    if isfield(data,'IVright'),      IVright      = double(data.IVright);      end
    if isfield(data,'trial2use')
        trial2use = double(data.trial2use);
        trial2use(~isfinite(trial2use)) = 0;
        trial2use(trial2use ~= 0) = 1;
    end

    % ---- Size check: transpose back to Nsub x N if needed ----
    if exist('MFpower','var')     && isequal(size(MFpower),     [N, Nsub]), MFpower     = MFpower.';     end
    if exist('ICPCpfc','var')     && isequal(size(ICPCpfc),     [N, Nsub]), ICPCpfc     = ICPCpfc.';     end
    if exist('l_ICPCmotor','var') && isequal(size(l_ICPCmotor), [N, Nsub]), l_ICPCmotor = l_ICPCmotor.'; end
    if exist('r_ICPCmotor','var') && isequal(size(r_ICPCmotor), [N, Nsub]), r_ICPCmotor = r_ICPCmotor.'; end
    if exist('IVleft','var')      && isequal(size(IVleft),      [N, Nsub]), IVleft      = IVleft.';      end
    if exist('IVright','var')     && isequal(size(IVright),     [N, Nsub]), IVright     = IVright.';     end
    if exist('trial2use','var')   && isequal(size(trial2use),   [N, Nsub]), trial2use   = trial2use.';   end

    % ---- Save output ----
    standir  = dirs.stan;
    baseVars = {'ya','a','s','r','rew','N','Nsub','sub','subIDs'};
    eegVars  = {};
    if exist('MFpower','var'),      eegVars{end+1} = 'MFpower';      end
    if exist('ICPCpfc','var'),      eegVars{end+1} = 'ICPCpfc';      end
    if exist('l_ICPCmotor','var'),  eegVars{end+1} = 'l_ICPCmotor';  end
    if exist('r_ICPCmotor','var'),  eegVars{end+1} = 'r_ICPCmotor';  end
    if exist('IVleft','var'),       eegVars{end+1} = 'IVleft';       end
    if exist('IVright','var'),      eegVars{end+1} = 'IVright';      end
    if exist('trial2use','var'),    eegVars{end+1} = 'trial2use';    end

    save(standir, baseVars{:}, eegVars{:}, '-v7');

    fprintf('EEG4stan.mat has been written to: %s\n', standir);
    whos('-file', standir)
end