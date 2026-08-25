function trials = importDB_fraser_v001(dataRoot, opts)
% importDB_fraser_v001  Import Fraser (Analysis_2026) data for Stage 5
% empirical validation. Modelled directly on importDB_pilot_v001 -- same
% structure, same conventions, differences noted inline.
%
% Brain2Bee MSc confirmatory test-retest study (Analysis_2026), a separate,
% later-recruited cohort from Pilot (same StimuliApp software, same study
% design). 48 raw participant folders copied (45 complete 6-session/2-visit,
% 3 excluded: ID101 5/6, ID102 5/6, ID214 3/6 -- see data/fraser/README.md
% point 2). Device: iPad Pro 13" (M4), confirmed 2026-08-06 to be the same
% physical unit used for Pilot (previously both this project's device labels
% were wrong -- see data/pilot/README.md and data/fraser/README.md for the
% correction trail).
%
% Source data/fraser/automatic/ is Ellipse-only (imageNumber==3) already --
% extractFraserEllipseData_v001.m filtered the other 4 shapes out at copy
% time (they exist in Analysis_2026/data/raw/automatic/ if ever needed).
% Requesting Shapes other than 3 here will simply return zero trials for
% that shape, not error -- the Shapes option is kept for interface parity
% with importDB_pilot_v001 and in case more shapes get added later.
%
% PX_PER_MM
%   10.41793 -- derived directly from the confirmed M4 hardware geometry
%   (2752 x 2064 px, 4:3, 13.0" diagonal exactly, ppi=264.615 unrounded;
%   see data/fraser/README.md point 3). This is the ONLY defensible value:
%   there was never a genuinely different, lower-resolution device in this
%   project -- Pilot's "12.9" 4th gen" label was always wrong, the same M4
%   hardware was used throughout. Do NOT use Pilot's 9.73 here.
%
%   A 2026-08-06 check comparing raw-px ellipse trace size between Pilot and
%   Fraser samples (2-3% apart) was initially misread as evidence for 9.73.
%   It isn't: since both datasets come from the identical device, that
%   comparison could never discriminate between candidate px/mm values --
%   there was no second, differently-scaled device to contrast against. What
%   it actually shows is that x,y are recorded 1:1 with real hardware pixels
%   (no hidden points/hi-DPI transform), which supports using the hardware-
%   derived figure rather than undermining it. The original ~120mm physical
%   check (2026-07-24, Pilot) was a coarse sanity check only (trial traces
%   span 101.7-134.1mm) -- not precise enough to arbitrate between 9.73 and
%   10.41793 either.
%
% SAMPLING RATE ASSUMPTION
%   Same as Pilot: touch_index is the reliable integer sample counter;
%   builds a uniform time vector t = (0:N-1)'/240. Analysis_2026's claude.md
%   documents a real pen-lift clock bug affecting this assumption for
%   Cursive trials specifically (single naive touch_index/fs clock
%   compresses pen-up gaps into one sample) -- confirmed NOT applicable to
%   Ellipse (single continuous stroke, 0% of templated trials affected), so
%   the Pilot-style uniform clock is safe to reuse here unmodified.
%
% SUBSAMPLING (caller's responsibility -- NOT done here)
%   Same as Pilot -- see importDB_pilot_v001 header for the stride pattern.
%
% COORDINATE CONVENTION
%   Raw iOS coordinates have y increasing downward (origin top-left). This
%   function negates y so that y increases upward, matching Pilot / Cook /
%   Hickman / Zarandi / the synthetic generator.
%
% SUBJECT ID SCHEME -- differs from Pilot, read before cross-referencing
%   Pilot's ID* folders are 6-digit RPS numbers (ID112182 -> subjectID
%   'P112182'). This cohort's are 3-digit, cohort-coded (ID101-ID109,
%   ID201-ID219, ID300-ID321 -- leading digit marks which MSc student's
%   participant block, not a magnitude). subjectID here uses an 'F' prefix
%   (e.g. 'F101') rather than Pilot's 'P' prefix, specifically so the two
%   cohorts are never visually confusable in downstream pooled tables even
%   though the numeric ranges don't actually overlap (100-999 vs 100000+).
%
% SHAPE LOOKUP (Huh & Sejnowski 2015, same taxonomy as Pilot / Cook / Hickman)
%   imageNumber  nu      Description
%   3            2       Ellipse (canonical power-law shape) -- ONLY shape
%                        present in data/fraser/automatic/ currently.
%
%   Theoretical beta: (2/3)*(1 + nu^2/2) / (1 + nu^2 + nu^4/15)
%
% USAGE
%   trials = importDB_fraser_v001()
%   trials = importDB_fraser_v001(dataRoot)
%   trials = importDB_fraser_v001(dataRoot, Visits=[1 2], Verbose=true)
%
% INPUTS
%   dataRoot       Path to data/fraser/.  Defaults to auto-detect via mfilename.
%
% NAME-VALUE OPTIONS
%   Shapes         imageNumber(s) to import (default: 3, ellipse only --
%                  the only shape this folder currently contains).
%   Visits         Which visits to include: [1], [2], or [1 2] (default).
%   Sessions       sessionInVisit values to include: subset of 1:3 (default: 1:3).
%   MinSamples     Minimum samples after edge clipping (default: 200).
%   EdgeClip       Samples removed from each end per trial (default: 50).
%   MaxSubjects    Cap number of participants (default: Inf; useful for smoke tests).
%   Verbose        Print progress (default: true).
%
% OUTPUT
%   trials   [M x 1] struct array.  Fields: same as importDB_pilot_v001,
%     with .subjectID  'F<subjNum>'  (e.g. 'F101')
%          .database    'fraser'
%          .SigmaToMM   1/10.41793  (M4 hardware-derived -- see PX_PER_MM note above)
%
% See also: importDB_pilot_v001, differentiateKinematicsEBR, regressDataEBR,
%           extractFraserEllipseData_v001
%
% Fraser, D.S. (2026)  v001

    arguments
        dataRoot         (1,:) char   = defaultFraserDataRoot()
        opts.Shapes      (1,:) double = 3
        opts.Visits      (1,:) double = [1 2]
        opts.Sessions    (1,:) double = 1:3
        opts.MinSamples  (1,1) double = 200
        opts.EdgeClip    (1,1) double = 50
        opts.MaxSubjects (1,1) double = Inf
        opts.Verbose     (1,1) logical = true
    end

    %% Constants
    FS          = 240;
    PX_PER_MM   = 10.41793;   % M4 hardware-derived -- see header note. NOT Pilot's 9.73.
    SIGMA_TO_MM = 1 / PX_PER_MM;

    % Angular frequency lookup per Huh & Sejnowski (2015) taxonomy
    % Keys: imageNumber 1-5 (only 3 has data in this folder currently)
    nuByImage = [2/5, 4/5, 2, 4, NaN];

    %% Locate automatic/ subfolder
    dataDir = fullfile(dataRoot, 'automatic');
    if ~isfolder(dataDir)
        error('importDB_fraser:noDir', '%s', ...
            ['automatic/ subfolder not found in: ' dataRoot]);
    end

    %% Enumerate participant ID folders
    entries   = dir(dataDir);
    isDirMask = [entries.isdir];
    isIDMask  = startsWith({entries.name}, 'ID');
    idFolders = sort({entries(isDirMask & isIDMask).name});

    nSubjTotal = numel(idFolders);
    if opts.MaxSubjects < nSubjTotal
        idFolders = idFolders(1 : opts.MaxSubjects);
    end

    if opts.Verbose
        fprintf('[importDB_fraser] Found %d participant folders', nSubjTotal);
        if opts.MaxSubjects < nSubjTotal
            fprintf(' (capped at %d)', opts.MaxSubjects);
        end
        fprintf('\n  Shapes: [%s]  Visits: [%s]  Sessions: [%s]\n', ...
            num2str(opts.Shapes), num2str(opts.Visits), num2str(opts.Sessions));
        if any(opts.Shapes ~= 3)
            fprintf('  NOTE: only shapeCode 3 (Ellipse) has data in this folder;\n');
            fprintf('        other requested shapes will silently return zero trials.\n');
        end
    end

    %% Collect trials
    trialCell  = {};
    nSubjects  = 0;
    nRejected  = 0;
    rejectLog  = {};

    for iSubj = 1:numel(idFolders)
        folderName = idFolders{iSubj};
        subjNum    = str2double(folderName(3:end));   % strip 'ID' prefix
        subjID     = sprintf('F%d', subjNum);          % 'F' prefix -- see header note
        subjDir    = fullfile(dataDir, folderName);

        %% Find session datetime stems from metadata CSVs
        metaPattern = fullfile(subjDir, 'betaTestv002 * tracingSection.csv');
        metaFiles   = dir(metaPattern);

        if isempty(metaFiles)
            continue
        end

        stems = cell(numel(metaFiles), 1);
        for f = 1:numel(metaFiles)
            tok = regexp(metaFiles(f).name, ...
                '^betaTestv002 (\d{4}-\d{2}-\d{2} \d{6}) tracingSection\.csv$', ...
                'tokens', 'once');
            if ~isempty(tok)
                stems{f} = tok{1};
            end
        end
        stems = unique(stems(~cellfun(@isempty, stems)));
        stems = sort(stems);    % chronological (ISO datetime sorts correctly)

        nSessions = numel(stems);
        if nSessions == 0
            continue
        end

        %% Assign visit and sessionInVisit
        % Sessions 1-3 -> Visit 1, sessions 4-6 -> Visit 2. Works unchanged
        % for the 3 incomplete subjects here (nSessions = 3 or 5): whatever
        % sessions exist are simply assigned in chronological order, same
        % as importDB_pilot_v001 -- no special-casing needed.
        visitNum      = zeros(nSessions, 1);
        sessionInVisit = zeros(nSessions, 1);
        for s = 1:nSessions
            if s <= 3
                visitNum(s)       = 1;
                sessionInVisit(s) = s;
            else
                visitNum(s)       = 2;
                sessionInVisit(s) = s - 3;
            end
        end

        nSubjects = nSubjects + 1;

        %% Loop over sessions
        for iSess = 1:nSessions
            vNum   = visitNum(iSess);
            sivNum = sessionInVisit(iSess);

            if ~ismember(vNum,   opts.Visits)   || ...
               ~ismember(sivNum, opts.Sessions)
                continue
            end

            stem     = stems{iSess};
            metaFile = fullfile(subjDir, ...
                ['betaTestv002 ' stem ' tracingSection.csv']);
            trajFile = fullfile(subjDir, ...
                ['betaTestv002 ' stem ' tracingSection-enhanced.csv']);

            if ~isfile(metaFile)
                rejectLog{end+1} = sprintf('%s sess%d: metadata CSV missing', ...
                    subjID, iSess);  %#ok<AGROW>
                continue
            end
            if ~isfile(trajFile)
                rejectLog{end+1} = sprintf('%s sess%d: trajectory CSV missing', ...
                    subjID, iSess);  %#ok<AGROW>
                continue
            end

            %% Load session files
            opts_read = detectImportOptions(metaFile);
            metaTbl   = readtable(metaFile, opts_read);

            opts_traj = detectImportOptions(trajFile);
            trajTbl   = readtable(trajFile, opts_traj);

            %% Process each requested shape
            for iShape = 1:numel(opts.Shapes)
                shapeCode = opts.Shapes(iShape);

                if shapeCode < 1 || shapeCode > 5
                    error('importDB_fraser:badShape', '%s', ...
                        sprintf('imageNumber must be 1-5, got %d', shapeCode));
                end

                nuVal     = nuByImage(shapeCode);
                shapeStr  = sprintf('nu%.4f_img%d', nuVal, shapeCode);
                if isnan(nuVal)
                    shapeStr = sprintf('cursive_img%d', shapeCode);
                end

                % Trial numbers in metadata matching this shape
                shapeMask = metaTbl.tracingScene_object1_imageNumber == shapeCode;
                trialNums = metaTbl.trial(shapeMask);

                for iTr = 1:numel(trialNums)
                    trNum = trialNums(iTr);

                    %% Extract trajectory rows for this trial
                    mask  = trajTbl.response_index == trNum;
                    tData = trajTbl(mask, :);

                    if isempty(tData)
                        rejectLog{end+1} = sprintf( ...
                            '%s V%dS%d T%02d sh%d: no trajectory rows', ...
                            subjID, vNum, sivNum, trNum, shapeCode); %#ok<AGROW>
                        nRejected = nRejected + 1;
                        continue
                    end

                    %% Sort by touch_index (no deduplication needed)
                    % Same assumption as importDB_pilot_v001: touch_index is
                    % a per-trial monotonic counter; uniform 240 Hz assumed.
                    tData = sortrows(tData, 'touch_index');
                    N     = height(tData);

                    %% Edge clip
                    clip = opts.EdgeClip;
                    if N <= 2 * clip
                        rejectLog{end+1} = sprintf( ...
                            '%s V%dS%d T%02d sh%d: too short before clip (%d samp)', ...
                            subjID, vNum, sivNum, trNum, shapeCode, N); %#ok<AGROW>
                        nRejected = nRejected + 1;
                        continue
                    end
                    tData = tData(clip+1 : end-clip, :);
                    N     = height(tData);

                    if N < opts.MinSamples
                        rejectLog{end+1} = sprintf( ...
                            '%s V%dS%d T%02d sh%d: too short after clip (%d samp)', ...
                            subjID, vNum, sivNum, trNum, shapeCode, N); %#ok<AGROW>
                        nRejected = nRejected + 1;
                        continue
                    end

                    %% Build coordinates
                    % Uniform 240 Hz time axis
                    tVec = (0:N-1)' / FS;

                    % Negate y: iOS y increases downward; standard convention y-up
                    xPx  = double(tData.x);
                    yPx  = -double(tData.y);

                    %% Assemble trial struct
                    trial = struct();
                    trial.x               = xPx;
                    trial.y               = yPx;
                    trial.z               = [];
                    trial.t               = tVec;
                    trial.fs              = FS;
                    trial.subjectID       = subjID;
                    trial.trialID         = sprintf('%s_V%d_S%d_T%03d_sh%d', ...
                                               subjID, vNum, sivNum, trNum, shapeCode);
                    trial.condition       = sprintf('VISIT%d', vNum);
                    trial.shape           = shapeStr;
                    trial.database        = 'fraser';
                    trial.notes           = sprintf( ...
                        'stem=%s visit=%d sessionInVisit=%d sessionOverall=%d', ...
                        stem, vNum, sivNum, iSess);
                    trial.units_verified  = false;
                    trial.SigmaToMM       = SIGMA_TO_MM;
                    trial.visit           = vNum;
                    trial.sessionInVisit  = sivNum;
                    trial.sessionOverall  = iSess;
                    trial.shapeCode       = shapeCode;
                    trial.nu              = nuVal;
                    trial.edgeClip        = clip;

                    trialCell{end+1} = trial;  %#ok<AGROW>

                end  % iTr
            end  % iShape
        end  % iSess
    end  % iSubj

    %% Assemble output
    if isempty(trialCell)
        warning('importDB_fraser:noTrials', ...
            'No valid trials extracted for Shapes=[%s] Visits=[%s].', ...
            num2str(opts.Shapes), num2str(opts.Visits));
        trials = struct([]);
        return
    end

    trials = [trialCell{:}]';

    %% Summary
    if opts.Verbose
        fprintf('\n[importDB_fraser] %d trials from %d participants', ...
            numel(trials), nSubjects);
        fprintf('  (%d rejected)\n', nRejected);

        uSubj  = unique({trials.subjectID});
        counts = arrayfun(@(s) sum(strcmp({trials.subjectID}, s)), uSubj);
        fprintf('  Trials per subject: %.1f mean, %d-%d range\n', ...
            mean(counts), min(counts), max(counts));

        uVisit = unique([trials.visit]);
        for v = uVisit
            nv = sum([trials.visit] == v);
            fprintf('  Visit %d: %d trials\n', v, nv);
        end

        if ~isempty(rejectLog)
            fprintf('  Rejection log (%d entries):\n', numel(rejectLog));
            for rr = 1:min(numel(rejectLog), 20)
                fprintf('    %s\n', rejectLog{rr});
            end
            if numel(rejectLog) > 20
                fprintf('    ... and %d more\n', numel(rejectLog) - 20);
            end
        end
    end

end

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function dr = defaultFraserDataRoot()
% DEFAULTFRASERDATAROOT  Resolve default path relative to this file.
% Layout: src/functions/importDB_fraser_v001.m -> ../../data/fraser
    thisDir  = fileparts(mfilename('fullpath'));   % .../src/functions
    srcDir   = fileparts(thisDir);                 % .../src
    repoRoot = fileparts(srcDir);                  % .../PowerLawSimulationPreReg
    dr       = fullfile(repoRoot, 'data', 'fraser');
end
