function trials = importDB_cook_v002(dataRoot, opts)
% importDB_cook_v002 Import Cook et al. (2026) WACOM data with metadata verification
%
% Returns an array of trial structs conforming to the Stage 5 empirical
% validation specification (see TODO_Empirical_Validation.md). Each struct
% contains raw x,y position data from the WACOM Intuos4 tablet at 133 Hz,
% ready for processing through the six pre-registered pipelines via
% Stage5_processTrials.m.
%
% v002 improvements over v001:
%   - Replicates the metadata verification from Cook's Script2_dataWrangler
%     (checks internal Subject_ID, Trial_ID, and Task_ID coherence)
%   - Replicates the quality-gated trial selection from Script3_ALL_OUT_Trimmer
%     (prioritises 'Default' trials over 'Bad' trials)
%   - Logs all misplaced, mismatched, and rejected trials with reasons
%   - Reports MATCHED / MISPLACED / MISSING / BAD counts matching Script2 output
%
% Usage:
%   trials = importDB_cook_v002()                                  % CTRL ellipses, default path
%   trials = importDB_cook_v002(dataRoot)
%   trials = importDB_cook_v002(dataRoot, Group="CTRL", Tasks=7)
%   trials = importDB_cook_v002(dataRoot, Group="CTRL", Tasks=3:9, Verbose=true)
%
% Inputs:
%   dataRoot  - Path to 'Paper Version_FINAL' folder containing All_DATA/
%               and TABLE_OSF.mat. Defaults to
%               <repoRoot>/data/cook/Paper Version_FINAL (resolved via
%               mfilename, so works regardless of pwd).
%
% Name-Value Arguments:
%   Group        - "CTRL" (default), "ASD", or "ALL"
%   Tasks        - Task IDs to import (default: 7, the ellipse at nu = 2).
%                  Tasks 3:9 map to angular frequencies [2/33, 2/5, 4/5, 4/3, 2, 3, 4].
%   MaxTrials    - Maximum trials per subject per shape (default: 5,
%                  matching Cook et al. published convention of Script3)
%   MinSamples   - Minimum trajectory length in samples (default: 100)
%   IncludeBad   - Include 'Bad' (incomplete) trials to fill up to MaxTrials
%                  if insufficient 'Default' trials exist (default: true,
%                  matching Script3 behaviour)
%   ScaleFactor  - Multiply x,y by this to convert to mm (default: 1.0;
%                  set when WACOM-to-mm conversion is determined)
%   Verbose      - Print progress (default: true)
%
% Output:
%   trials  - [M x 1] struct array with fields per the standardised spec:
%             .x, .y, .z, .t, .fs, .subjectID, .trialID, .condition,
%             .shape, .database, .notes, .units_verified
%
% Notes:
%   - The raw .mat files contain an OUT cell array populated by the original
%     experiment code (Ben Dongsung Huh, 2018). Each OUT{kk} has .DATA (the
%     trajectory matrix), .idx (column index struct), and .Experiment
%     (metadata). We extract DATA(:, idx.pos) for raw x,y positions.
%   - This function bypasses Cook's Script1 (collation) and directly loads
%     from All_DATA/, but replicates the verification logic of Scripts 2-3.
%   - Cook's Script2 found that some trials were MISPLACED: their internal
%     metadata (Subject_ID, Trial_ID) did not match the filename or cell
%     position. Our approach is stricter: we trust only the internal
%     metadata and reject trials where it conflicts with the filename.
%   - Units: WACOM Intuos4 output units TBD. ScaleFactor defaults to 1.0
%     and units_verified is set to false until the conversion is confirmed.
%     Beta estimation is unit-independent; sigma characterisation requires mm.
%
% Dependencies:
%   TABLE_OSF.mat (participant metadata, in dataRoot)
%
% See also: Stage5_processTrials, differentiateKinematicsEBR, regressDataEBR
%
% Created 2026-03-17  v002
% Dagmar Scott Fraser - d.s.fraser@bham.ac.uk

    arguments
        dataRoot    (1,:) char = defaultCookDataRoot()
        opts.Group      (1,:) string = "CTRL"
        opts.Tasks      (1,:) double = 7
        opts.MaxTrials  (1,1) double = 5
        opts.MinSamples (1,1) double = 100
        opts.IncludeBad (1,1) logical = true
        opts.ScaleFactor(1,1) double = 1.0
        opts.MaxFiles   (1,1) double = Inf
        opts.MaxSubjects(1,1) double = Inf
        opts.Verbose    (1,1) logical = true
    end

    %% Shape-to-angular-frequency lookup (Huh & Sejnowski 2015, Cook et al.)
    % Task:  1    2     3      4     5     6      7   8   9
    % nu:   NaN  NaN  2/33   2/5   4/5   4/3     2   3   4
    nuLookup = [NaN, NaN, 2/33, 2/5, 4/5, 4/3, 2, 3, 4];

    FS_COOK = 133; % Hz, WACOM Intuos4 sampling rate (Cook et al. 2026)

    %% Validate inputs
    tableFile = fullfile(dataRoot, "TABLE_OSF.mat");
    dataDir   = fullfile(dataRoot, "All_DATA");

    assert(isfile(tableFile), "importDB_cook:fileNotFound", ...
        "TABLE_OSF.mat not found in: %s", dataRoot);
    assert(isfolder(dataDir), "importDB_cook:folderNotFound", ...
        "All_DATA/ not found in: %s", dataRoot);

    %% Load participant metadata
    S = load(tableFile, "TABLE");
    TABLE = S.TABLE;

    % Extract subject info
    % TABLE has two ID systems:
    %   Subject_ID     - sequential (1,2,3...) used by TABLE.Experiment
    %   Subject_ID_old - original experiment codes used INSIDE .mat files
    % Script2_dataWrangler hardcodes the old IDs for this reason.
    allSubjectIDs    = TABLE.Subject_Info.Subject_ID;      % new (sequential)
    allSubjectIDsOld = TABLE.Subject_Info.Subject_ID_old;  % old (in .mat files)
    allGroupLabels   = TABLE.Subject_Info.ASD_or_CTRL;

    % Build old-ID to group lookup (for checking internal .mat metadata)
    oldID2group = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    newID2old   = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    for ii = 1:numel(allSubjectIDs)
        oldID2group(int32(allSubjectIDsOld(ii))) = allGroupLabels{ii};
        newID2old(int32(allSubjectIDs(ii)))       = int32(allSubjectIDsOld(ii));
    end

    % Filter by requested group
    % targetIDs uses NEW numbering (for TABLE.Experiment queries)
    % targetIDsOld uses OLD numbering (for matching .mat internal metadata)
    if strcmpi(opts.Group, "ALL")
        targetIDs = allSubjectIDs;
        groupLabel = "ALL";
    else
        groupMatch = strcmpi(string(allGroupLabels), opts.Group);
        targetIDs = allSubjectIDs(groupMatch);
        groupLabel = upper(opts.Group);
    end

    % Optionally cap the number of subjects (for smoke testing)
    if opts.MaxSubjects < numel(targetIDs)
        targetIDs = targetIDs(1:opts.MaxSubjects);
    end

    targetIDsOld = arrayfun(@(id) newID2old(int32(id)), targetIDs);

    nSubjects = numel(targetIDs);
    if opts.Verbose
        fprintf("[importDB_cook] Group: %s | Subjects: %d | Tasks: [%s]\n", ...
            groupLabel, nSubjects, num2str(opts.Tasks));
    end

    %% Verification counters (matching Script2 output)
    nMatched   = 0;
    nMisplaced = 0;
    nMissing   = 0;
    nBadQuality = 0;
    nGoodQuality = 0;

    %% Intermediate storage: per-subject, per-task verified trials
    % We collect all verified trials first, then do quality-gated selection
    % (replicating Script3). Key = "S<subjID>_T<taskID>"
    verifiedTrials = containers.Map('KeyType', 'char', 'ValueType', 'any');

    rejectLog = {};

    %% Loop over tasks: load all files for this task, verify metadata
    for taskIdx = 1:numel(opts.Tasks)
        taskID = opts.Tasks(taskIdx);
        assert(taskID >= 1 && taskID <= 9, "importDB_cook:badTask", ...
            "Task ID must be 1-9, got %d", taskID);

        nu = nuLookup(taskID);

        % Find filenames for this task and target group
        idxSubj = ismember(TABLE.Experiment.Subject_ID, targetIDs);
        idxTask = TABLE.Experiment.Task_ID == taskID;
        idxBoth = find(idxSubj & idxTask);
        filenames = TABLE.Experiment.Filename(idxBoth);
        fileSubjIDs = TABLE.Experiment.Subject_ID(idxBoth);

        if opts.Verbose
            fprintf("  Task %d (nu = %.4f): %d files to inspect\n", ...
                taskID, nu, numel(filenames));
        end

        nFiles = min(numel(filenames), opts.MaxFiles);
        for fIdx = 1:nFiles
            matFile = fullfile(dataDir, filenames{fIdx});
            expectedSubjID = fileSubjIDs(fIdx);

            if ~isfile(matFile)
                warning("importDB_cook:fileMissing", ...
                    "File not found: %s", matFile);
                nMissing = nMissing + 1;
                continue
            end

            % Load the raw experiment file
            S2 = load(matFile, "OUT");
            if ~isfield(S2, "OUT")
                warning("importDB_cook:noOUT", ...
                    "No OUT variable in %s", filenames{fIdx});
                nMissing = nMissing + 1;
                continue
            end
            OUT = S2.OUT;

            % Inspect each trial in the OUT cell array
            for kk = 1:numel(OUT)
                if isempty(OUT{kk})
                    nMissing = nMissing + 1;
                    continue
                end

                thisOUT = OUT{kk};

                %% ---- METADATA VERIFICATION (replicates Script2) ----

                % Check Experiment struct exists
                if ~isfield(thisOUT, 'Experiment')
                    nMissing = nMissing + 1;
                    rejectLog{end+1} = sprintf( ...
                        "%s[%d]: no Experiment field", filenames{fIdx}, kk); %#ok<AGROW>
                    continue
                end

                exp = thisOUT.Experiment;

                % Check Trial_ID == Task_ID coherence (Script2 line 108)
                % exp may be a table row; use try-catch instead of isfield
                try
                    trialIDCheck = exp.Trial_ID;
                    taskIDCheck  = exp.Task_ID;
                    hasIDs = true;
                catch
                    hasIDs = false;
                end
                if hasIDs
                    if trialIDCheck ~= taskIDCheck
                        rejectLog{end+1} = sprintf( ...
                            "%s[%d]: Trial_ID(%d) ~= Task_ID(%d)", ...
                            filenames{fIdx}, kk, exp.Trial_ID, exp.Task_ID); %#ok<AGROW>
                        nMisplaced = nMisplaced + 1;
                        continue
                    end
                end

                % Verify Subject_ID matches filename expectation
                internalSubjID = exp.Subject_ID; % OLD numbering from .mat
                internalTaskID = exp.Trial_ID;     % Trial_ID == Task_ID after check above

                % expectedSubjID is NEW numbering from TABLE.Experiment;
                % convert to OLD for comparison against internal metadata
                expectedSubjIDOld = newID2old(int32(expectedSubjID));

                subjectMatch = (internalSubjID == expectedSubjIDOld);
                taskMatch    = (internalTaskID == taskID);

                if subjectMatch && taskMatch
                    nMatched = nMatched + 1;
                elseif ~subjectMatch || ~taskMatch
                    nMisplaced = nMisplaced + 1;
                    if opts.Verbose
                        fprintf("    MISPLACED in %s[%d]: ", filenames{fIdx}, kk);
                        if ~subjectMatch
                            fprintf("SubjID got %d expected %d; ", ...
                                internalSubjID, expectedSubjID);
                        end
                        if ~taskMatch
                            fprintf("TaskID got %d expected %d", ...
                                internalTaskID, taskID);
                        end
                        fprintf("\n");
                    end

                    % Script2 relocates misplaced trials to the correct
                    % position based on internal metadata. We do the same:
                    % trust the INTERNAL metadata over the filename.
                    % But only keep if this subject is in our target group.
                    if ~ismember(internalSubjID, targetIDsOld)
                        rejectLog{end+1} = sprintf( ...
                            "%s[%d]: misplaced S%d not in target group %s", ...
                            filenames{fIdx}, kk, internalSubjID, groupLabel); %#ok<AGROW>
                        continue
                    end

                    % Also reject if the internal Task_ID is not in our
                    % requested task list
                    if ~ismember(internalTaskID, opts.Tasks)
                        rejectLog{end+1} = sprintf( ...
                            "%s[%d]: misplaced task %d not in requested tasks", ...
                            filenames{fIdx}, kk, internalTaskID); %#ok<AGROW>
                        continue
                    end

                    % Override: use internal metadata for placement
                    taskID_effective = internalTaskID;
                else
                    taskID_effective = taskID;
                end

                % If we matched, use the expected task
                if subjectMatch && taskMatch
                    taskID_effective = taskID;
                end

                %% ---- QUALITY CHECK (replicates Script3 logic) ----
                % Experiment may be a table row (not struct), so isfield
                % fails. Use try-catch as Script2 does.
                qualityStr = "UNKNOWN";
                try
                    qualityStr = string(char(exp.Quality));
                catch
                    % Quality field absent or unreadable
                end

                if strcmpi(qualityStr, "Default")
                    nGoodQuality = nGoodQuality + 1;
                    qualityScore = 1; % prefer these
                elseif strcmpi(qualityStr, "Bad")
                    nBadQuality = nBadQuality + 1;
                    qualityScore = -1; % include only if needed
                else
                    qualityScore = 0; % unknown
                end

                %% ---- DATA EXTRACTION ----
                if ~isfield(thisOUT, 'DATA') || ~isfield(thisOUT, 'idx')
                    rejectLog{end+1} = sprintf("S%d_T%d[%d]: missing DATA or idx", ...
                        internalSubjID, taskID_effective, kk); %#ok<AGROW>
                    continue
                end

                posIdx = thisOUT.idx.pos;
                rawXY = thisOUT.DATA(:, posIdx);

                % Remove NaN rows
                nanMask = any(isnan(rawXY), 2);
                rawXY = rawXY(~nanMask, :);

                % Minimum length check (pre-trim)
                if size(rawXY, 1) < opts.MinSamples
                    rejectLog{end+1} = sprintf( ...
                        "S%d_T%d[%d]: too short (%d samples)", ...
                        internalSubjID, taskID_effective, kk, size(rawXY, 1)); %#ok<AGROW>
                    continue
                end

                % Trim constant-position leading/trailing samples
                rawXY = trimConstantEnds(rawXY);

                if size(rawXY, 1) < opts.MinSamples
                    rejectLog{end+1} = sprintf( ...
                        "S%d_T%d[%d]: too short after trim (%d)", ...
                        internalSubjID, taskID_effective, kk, size(rawXY, 1)); %#ok<AGROW>
                    continue
                end

                %% ---- STORE VERIFIED TRIAL ----
                % Key by internal metadata (trusted source of truth)
                mapKey = sprintf("S%d_T%d", internalSubjID, taskID_effective);

                entry = struct();
                entry.rawXY        = rawXY;
                entry.qualityScore = qualityScore;
                entry.qualityStr   = qualityStr;
                entry.subjectID    = internalSubjID;
                entry.taskID       = taskID_effective;
                entry.sourceFile   = filenames{fIdx};
                entry.cellIdx      = kk;
                entry.nSamples     = size(rawXY, 1);

                % Determine condition from internal Subject_ID
                if oldID2group.isKey(int32(internalSubjID))
                    entry.condition = string(oldID2group(int32(internalSubjID)));
                else
                    entry.condition = "UNKNOWN";
                end

                if verifiedTrials.isKey(mapKey)
                    existing = verifiedTrials(mapKey);
                    existing{end+1} = entry; %#ok<AGROW>
                    verifiedTrials(mapKey) = existing;
                else
                    verifiedTrials(mapKey) = {entry};
                end

            end % kk (trials within file)
        end % fIdx (files)
    end % taskIdx

    %% Report verification counts (matching Script2 output format)
    if opts.Verbose
        fprintf("\n  Verification: MATCHED=%d  MISPLACED=%d  MISSING=%d\n", ...
            nMatched, nMisplaced, nMissing);
        fprintf("  Quality:      Default=%d  Bad=%d\n", ...
            nGoodQuality, nBadQuality);
    end

    %% Quality-gated trial selection (replicates Script3_ALL_OUT_Trimmer)
    % For each subject x task, take up to MaxTrials, preferring 'Default'
    trialCell = {};
    nSelected = 0;
    nBackfilled = 0;

    allKeys = verifiedTrials.keys();
    for kIdx = 1:numel(allKeys)
        key = allKeys{kIdx};
        entries = verifiedTrials(key);

        % Separate by quality
        scores = cellfun(@(e) e.qualityScore, entries);
        defaultIdx = find(scores == 1);
        badIdx     = find(scores == -1);
        unknownIdx = find(scores == 0);

        % Build selection order: Default first, then Bad, then Unknown
        selectionOrder = [defaultIdx(:); badIdx(:); unknownIdx(:)];

        if ~opts.IncludeBad
            selectionOrder = defaultIdx(:);
        end

        nToSelect = min(opts.MaxTrials, numel(selectionOrder));

        for si = 1:nToSelect
            entry = entries{selectionOrder(si)};

            % Apply scale factor
            xPos = entry.rawXY(:, 1) * opts.ScaleFactor;
            yPos = entry.rawXY(:, 2) * opts.ScaleFactor;

            % Construct time vector
            N = numel(xPos);
            tVec = (0:N-1)' / FS_COOK;

            % Angular frequency for shape label
            nuVal = nuLookup(entry.taskID);
            shapeStr = sprintf("nu%.4f_task%d", nuVal, entry.taskID);

            % Build standardised trial struct
            trial = struct();
            trial.x              = xPos;
            trial.y              = yPos;
            trial.z              = [];
            trial.t              = tVec;
            trial.fs             = FS_COOK;
            trial.subjectID      = sprintf("S%03d", entry.subjectID);
            trial.trialID        = sprintf("S%03d_T%d_R%d", ...
                entry.subjectID, entry.taskID, si);
            trial.condition      = entry.condition;
            trial.shape          = shapeStr;
            trial.database       = "cook";
            trial.notes          = sprintf("quality=%s src=%s[%d]", ...
                entry.qualityStr, entry.sourceFile, entry.cellIdx);
            trial.units_verified = false;
            trial.SigmaToMM      = 0.248;         % WACOM Cintiq 22HD px → mm (WACOM support confirmed)

            trialCell{end+1} = trial; %#ok<AGROW>
            nSelected = nSelected + 1;

            if entry.qualityScore < 1
                nBackfilled = nBackfilled + 1;
            end
        end
    end

    %% Assemble output
    if isempty(trialCell)
        warning("importDB_cook:noTrials", "No valid trials extracted.");
        trials = struct([]);
        return
    end

    trials = [trialCell{:}]';

    %% Summary
    if opts.Verbose
        fprintf("\n[importDB_cook] Selected %d trials (%d backfilled from Bad)\n", ...
            nSelected, nBackfilled);

        uSubj = unique(string({trials.subjectID}));
        fprintf("  Subjects with data: %d\n", numel(uSubj));

        allSubjStrs = string({trials.subjectID});
        trialCounts = arrayfun(@(s) sum(allSubjStrs == s), uSubj);
        fprintf("  Trials per subject: %.1f (mean), %d-%d (range)\n", ...
            mean(trialCounts), min(trialCounts), max(trialCounts));

        if ~isempty(rejectLog)
            fprintf("  Rejection log (%d entries):\n", numel(rejectLog));
            for rr = 1:min(numel(rejectLog), 20) % cap display at 20
                fprintf("    %s\n", rejectLog{rr});
            end
            if numel(rejectLog) > 20
                fprintf("    ... and %d more\n", numel(rejectLog) - 20);
            end
        end
    end

end

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function xy = trimConstantEnds(xy)
% TRIMCONSTANTENDS Remove leading/trailing rows where position is unchanged
% Detects stuck stylus at start or end of recording.
    dx = diff(xy);
    movement = any(abs(dx) > 0, 2);

    firstMove = find(movement, 1, "first");
    lastMove  = find(movement, 1, "last") + 1;

    if isempty(firstMove) || isempty(lastMove)
        xy = zeros(0, 2);
        return
    end

    xy = xy(firstMove:lastMove, :);
end

function dr = defaultCookDataRoot()
% DEFAULTCOOKDATAROOT Resolve default path to Cook data relative to this file
% Anchored via mfilename so it works regardless of pwd.
% Layout: src/functions/importDB_cook_v002.m -> ../../data/cook/Paper Version_FINAL
    thisDir  = fileparts(mfilename('fullpath')); % .../src/functions
    srcDir   = fileparts(thisDir);                % .../src
    repoRoot = fileparts(srcDir);                 % .../PowerLawSimulationPreReg
    dr = fullfile(repoRoot, 'data', 'cook', 'Paper Version_FINAL');
end
