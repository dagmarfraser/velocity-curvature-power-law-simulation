function trials = importDB_pilot_v001(dataRoot, opts)
% importDB_pilot_v001  Import iPad Pro pilot data for Stage 5 empirical validation.
%
% Brain2Bee MSc test-retest study (Fraser, University of Birmingham).
% N = 61 healthy adults, iPad Pro 13" (M4), Apple Pencil, 240 Hz nominal.
% (Device label corrected 2026-08-06; was mislabelled 12.9" 4th gen. The
% PX_PER_MM=9.73 constant below is unaffected -- independently re-confirmed
% by direct physical measurement, not derived from a spec sheet.)
% 2 visits ~1 week apart, 3 sessions per visit, 50 trials per session.
% 5 shapes per H&S angular-frequency taxonomy (imageNumbers 1-5, nu = 2/5,
% 4/5, 2, 4, NaN).  No intervention: Visit 1 / Visit 2 are test and retest.
%
% SAMPLING RATE ASSUMPTION
%   The Apple Pencil hardware runs at 240 Hz.  iOS timestamps (time_seconds)
%   are truncated to ~3 significant figures, causing occasional duplicate
%   timestamp rows for distinct samples.  touch_index is the reliable integer
%   sample counter.  This function deduplicates on touch_index, sorts, then
%   builds a uniform time vector  t = (0:N-1)'/240.
%
% SUBSAMPLING (caller's responsibility -- NOT done here)
%   To examine the dataset at lower rates without resampling artefacts:
%     120 Hz subset A:  x(1:2:end), y(1:2:end), t(1:2:end)
%     120 Hz subset B:  x(2:2:end), y(2:2:end), t(2:2:end)
%      60 Hz subset A:  x(1:4:end), ...
%      60 Hz subset B:  x(2:4:end), ...
%      60 Hz subset C:  x(3:4:end), ...
%      60 Hz subset D:  x(4:4:end), ...
%   The simulation grid contains fs = [60 120 240] Hz, so all three rates
%   are directly comparable to the synthetic parameter space.
%
% COORDINATE CONVENTION
%   Raw iOS coordinates have y increasing downward (origin top-left).
%   This function negates y so that y increases upward, matching the
%   convention used by Cook, Hickman, Zarandi, and the synthetic generator.
%
% UNITS
%   Coordinates are stored in pixels.  Beta is scale-invariant; no
%   px-to-mm conversion is applied.  trial.SigmaToMM = 1/9.73 is stored
%   so that noise characterisation (characteriseBiologicalNoise_v003) can
%   convert sigma from px to mm before constellation lookup.
%
% SHAPE LOOKUP (Huh & Sejnowski 2015, same taxonomy as Cook / Hickman)
%   imageNumber  nu      Description
%   1            2/5     Lissajous (low angular frequency)
%   2            4/5     Lissajous (mid angular frequency)
%   3            2       Ellipse (canonical power-law shape)
%   4            4       Superellipse
%   5            NaN     Cursive 'l'
%
%   Theoretical beta: (2/3)*(1 + nu^2/2) / (1 + nu^2 + nu^4/15)
%   For cursive 'l' the reference is beta = 1/3.
%
% USAGE
%   trials = importDB_pilot_v001()
%   trials = importDB_pilot_v001(dataRoot)
%   trials = importDB_pilot_v001(dataRoot, Shapes=3)
%   trials = importDB_pilot_v001(dataRoot, Shapes=3, Visits=[1 2], Verbose=true)
%
% INPUTS
%   dataRoot       Path to data/pilot/.  Defaults to auto-detect via mfilename.
%
% NAME-VALUE OPTIONS
%   Shapes         imageNumber(s) to import (default: 3, ellipse only).
%   Visits         Which visits to include: [1], [2], or [1 2] (default).
%   Sessions       sessionInVisit values to include: subset of 1:3 (default: 1:3).
%   MinSamples     Minimum samples after edge clipping (default: 200).
%   EdgeClip       Samples removed from each end per trial (default: 50).
%   MaxSubjects    Cap number of participants (default: Inf; useful for smoke tests).
%   Verbose        Print progress (default: true).
%
% OUTPUT
%   trials   [M x 1] struct array.  Fields:
%     .x, .y           Position in pixels; y negated to standard convention.
%     .z               [] (2-D tablet data)
%     .t               Uniform time vector, seconds (0:N-1)'/240
%     .fs              240
%     .subjectID       'P<RPS_ID>'  (e.g. 'P112321')
%     .trialID         'P112321_V1_S2_T045_sh3'  (unique per trial)
%     .condition       'VISIT1' or 'VISIT2'
%     .shape           'nu2.0000_img3'  (nu-coded, matching Cook naming)
%     .database        'pilot'
%     .notes           stem / visit / session info string
%     .units_verified  false  (native px)
%     .SigmaToMM       1/9.73  (~0.1028 mm/px)
%     .visit           1 or 2
%     .sessionInVisit  1, 2, or 3
%     .sessionOverall  1-6  (1-3 = Visit 1, 4-6 = Visit 2)
%     .shapeCode       imageNumber (1-5)
%     .nu              angular frequency per H&S (NaN for cursive)
%     .edgeClip        samples clipped per end (for reference)
%
% See also: characteriseBiologicalNoise_v003, constellationHickman_v001,
%           differentiateKinematicsEBR, regressDataEBR
%
% Fraser, D.S. (2026)  v001

    arguments
        dataRoot         (1,:) char   = defaultPilotDataRoot()
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
    PX_PER_MM   = 9.73;   % iPad Pro 13" M4 (Fraser et al. calibration; label corrected 2026-08-06)
    SIGMA_TO_MM = 1 / PX_PER_MM;

    % Angular frequency lookup per Huh & Sejnowski (2015) taxonomy
    % Keys: imageNumber 1-5
    nuByImage = [2/5, 4/5, 2, 4, NaN];

    %% Locate automatic/ subfolder
    dataDir = fullfile(dataRoot, 'automatic');
    if ~isfolder(dataDir)
        error('importDB_pilot:noDir', '%s', ...
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
        fprintf('[importDB_pilot] Found %d participant folders', nSubjTotal);
        if opts.MaxSubjects < nSubjTotal
            fprintf(' (capped at %d)', opts.MaxSubjects);
        end
        fprintf('\n  Shapes: [%s]  Visits: [%s]  Sessions: [%s]\n', ...
            num2str(opts.Shapes), num2str(opts.Visits), num2str(opts.Sessions));
    end

    %% Collect trials
    trialCell  = {};
    nSubjects  = 0;
    nRejected  = 0;
    rejectLog  = {};

    for iSubj = 1:numel(idFolders)
        folderName = idFolders{iSubj};
        subjNum    = str2double(folderName(3:end));   % strip 'ID' prefix
        subjID     = sprintf('P%d', subjNum);
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
        % Sessions 1-3 → Visit 1 (chronologically first visit)
        % Sessions 4-6 → Visit 2 (second visit, ~1 week later)
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
                    error('importDB_pilot:badShape', '%s', ...
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
                    % touch_index is a per-trial monotonic counter that resets
                    % to 1 at the start of each trial (NOT a session counter).
                    % Empirically confirmed: 0/50 trials have intra-trial
                    % duplicate touch_index values (verified 2026-04-28).
                    % time_seconds has iOS truncation artefacts (~3 s.f.) giving
                    % ~2-19 colliding timestamps per trial -- ignore for time axis.
                    % Uniform 240 Hz assumed; t = (0:N-1)'/240 built below.
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
                    trial.database        = 'pilot';
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
        warning('importDB_pilot:noTrials', ...
            'No valid trials extracted for Shapes=[%s] Visits=[%s].', ...
            num2str(opts.Shapes), num2str(opts.Visits));
        trials = struct([]);
        return
    end

    trials = [trialCell{:}]';

    %% Summary
    if opts.Verbose
        fprintf('\n[importDB_pilot] %d trials from %d participants', ...
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

function dr = defaultPilotDataRoot()
% DEFAULTPILOTDATAROOT  Resolve default path relative to this file.
% Layout: src/functions/importDB_pilot_v001.m -> ../../data/pilot
    thisDir  = fileparts(mfilename('fullpath'));   % .../src/functions
    srcDir   = fileparts(thisDir);                 % .../src
    repoRoot = fileparts(srcDir);                  % .../PowerLawSimulationPreReg
    dr       = fullfile(repoRoot, 'data', 'pilot');
end
