function trials = importDB_hickman_v002(opts)
% importDB_hickman_v002 Import Hickman et al. (2024) from raw sources.
%
% Study 1 (PD ON/OFF): reads raw Pavlovia CSVs (~38 Hz irregular),
%   resamples to 60 Hz via PCHIP on raw timestamps.
% Study 2 (Haloperidol/Placebo): reads per-trial .mat files, extracts
%   DATA0 rows 2-3 at the WACOM native rate of 133 Hz.
%
% USAGE:
%   trials = importDB_hickman_v002()
%   trials = importDB_hickman_v002(Study=1, Group="PDM0", Shapes=3)
%   trials = importDB_hickman_v002(Study=2, Group="HALO", Shapes=2)
%
% INPUTS (name-value):
%   Study     - 1 (PD ON/OFF) or 2 (Haloperidol/Placebo). Default: 1.
%   Group     - Study 1: "PDM0"|"PDM1"|"all". Study 2: "HALO"|"PLAC"|"all".
%               Default: "all".
%   Shapes    - Shape column indices. Default: 3 (ellipse, angfreq 2).
%               [1=4/5, 2=4/3, 3=ellipse, 4=angfreq4]
%   Verbose   - Print summary. Default: true.
%   DataRoot  - Path to data/hickman/. Default: auto-detect.
%
% CHANGES FROM v001:
%   - Study 1: reads raw Pavlovia CSVs instead of baseData.mat.
%     Resamples from irregular ~38 Hz to 60 Hz via PCHIP interpolation.
%   - Study 2: reads per-trial .mat files from Raw Data/ subfolders.
%     Extracts DATA0(2:3,:) at 133 Hz. Ignores DATA0(1,:) which is
%     NOT wall-clock time. Handles blinding via pID_P_H_file.mat.
%   - Both studies: positions in experiment-native units (not mm).
%
% PREREQ:
%   Study 1: shapesImport.m on path (in Supporting Scripts/).
%   Study 2: pID_P_H_file.mat in Study 2/ directory.
%
% Fraser, D.S. (2026)

    arguments
        opts.Study    (1,1) double {mustBeMember(opts.Study, [1 2])} = 1
        opts.Group    (1,1) string = "all"
        opts.Shapes   (1,:) double = 3
        opts.Verbose  (1,1) logical = true
        opts.DataRoot (1,1) string = ""
    end

    if opts.DataRoot == ""
        thisFile = mfilename('fullpath');
        if ~isempty(thisFile)
            opts.DataRoot = fullfile(fileparts(fileparts(thisFile)), ...
                '..', 'data', 'hickman');
        else
            error('importDB_hickman:noRoot', ...
                'Cannot auto-detect DataRoot. Provide DataRoot explicitly.');
        end
    end

    if opts.Study == 1
        trials = importStudy1(opts);
    else
        trials = importStudy2(opts);
    end
end

%% ========================================================================
%  STUDY 1: Pavlovia CSVs -> 60 Hz via PCHIP
%  ========================================================================
function trials = importStudy1(opts)

    studyDir = fullfile(opts.DataRoot, 'Study 1');
    assert(isfolder(studyDir), 'importDB_hickman:noDir', 'Not found: %s', studyDir);
    addpath(genpath(fullfile(studyDir, 'Supporting Scripts')));

    shapeNames = ["angfreq_4to5", "angfreq_4to3", "ellipse", "angfreq_4"];
    shapeIDs   = ["images/angfreq_4to5.png", "images/angfreq_4to3.png", ...
                  "images/angfreq_2.png",    "images/angfreq_4.png"];
    targetFs   = 60;
    minSamples = 100;

    condDirs = struct('folder', {"PD_OFF","PD_ON"}, ...
                      'prefix', {"PDM0","PDM1"}, ...
                      'label',  {"PD_OFF","PD_ON"});

    trials = [];

    for ci = 1:numel(condDirs)
        if opts.Group ~= "all" && opts.Group ~= string(condDirs(ci).prefix)
            continue
        end

        rawDir = fullfile(studyDir, 'Raw Data', condDirs(ci).folder);
        if ~isfolder(rawDir)
            warning('importDB_hickman:noRawDir', 'Not found: %s', rawDir);
            continue
        end
        csvFiles = dir(fullfile(rawDir, '*.csv'));
        if opts.Verbose
            fprintf('Study 1 %s: %d CSVs in %s\n', ...
                condDirs(ci).prefix, numel(csvFiles), rawDir);
        end

        for fi = 1:numel(csvFiles)
            csvPath = fullfile(rawDir, csvFiles(fi).name);
            pid = regexp(csvFiles(fi).name, '^\d+', 'match', 'once');
            shapes = shapesImport(csvPath);

            for row = 1:size(shapes, 1)
                sid = shapes{row, 8};
                if ~ischar(sid), continue; end

                shapeCol = find(shapeIDs == string(sid));
                if isempty(shapeCol) || ~ismember(shapeCol, opts.Shapes)
                    continue
                end
                if shapes{row, 7} == 0, continue; end

                try
                    tRaw = eval(shapes{row, 6});
                    xRaw = eval(shapes{row, 1}) * 1200;
                    yRaw = eval(shapes{row, 2}) * 1200;
                catch
                    continue
                end
                if numel(tRaw) < minSamples, continue; end

                % PCHIP resample to 60 Hz on raw timestamps
                tUniform = (tRaw(1):1/targetFs:tRaw(end))';
                xUniform = interp1(tRaw, xRaw, tUniform, 'pchip');
                yUniform = interp1(tRaw, yRaw, tUniform, 'pchip');
                if numel(tUniform) < minSamples, continue; end

                trialIdx = countTrialForPidShape(trials, pid, ...
                    condDirs(ci).prefix, shapeCol);

                tr = buildTrial(xUniform, yUniform, tUniform, targetFs, ...
                    pid, condDirs(ci).prefix, condDirs(ci).label, ...
                    shapeNames(shapeCol), shapeCol, trialIdx, 1, ...
                    sprintf('Raw CSV, PCHIP to %d Hz from ~%.0f Hz', ...
                    targetFs, 1/median(diff(tRaw))), 0.1133);
                trials = appendTrial(trials, tr);
            end
        end
    end
    if opts.Verbose, printSummary(trials, 1); end
end

%% ========================================================================
%  STUDY 2: per-trial .mat, DATA0 at 133 Hz
%  ========================================================================
function trials = importStudy2(opts)

    studyDir = fullfile(opts.DataRoot, 'Study 2');
    assert(isfolder(studyDir), 'importDB_hickman:noDir', 'Not found: %s', studyDir);

    shapeNames   = ["angfreq_4to5","angfreq_4to3","ellipse","angfreq_4"];
    shapeFreqStr = ["4/5","4/3","2","4"];
    nativeFs     = 133;
    minSamples   = 100;

    B = load(fullfile(studyDir, 'pID_P_H_file.mat'), 'pID_P_H');
    pID_P_H = B.pID_P_H;
    nSubj = size(pID_P_H, 1);

    trials = [];

    for k = 1:nSubj
        pidStr = pID_P_H{k, 1};
        for dayCol = 2:3
            intervention = pID_P_H{k, dayCol};
            if intervention == -1, continue; end

            if intervention == 0
                condLabel = "Placebo"; condPrefix = "PLAC";
            else
                condLabel = "Haloperidol"; condPrefix = "HALO";
            end
            if opts.Group ~= "all" && opts.Group ~= condPrefix
                continue
            end

            if dayCol == 2
                subjDir = fullfile(studyDir, 'Raw Data', pidStr);
            else
                subjDir = fullfile(studyDir, 'Raw Data', pidStr, 'DAY2');
            end
            if ~isfolder(subjDir), continue; end

            matFiles = dir(fullfile(subjDir, '*.mat'));
            for fi = 1:numel(matFiles)
                S = load(fullfile(subjDir, matFiles(fi).name), 'OUT');
                if ~isfield(S, 'OUT'), continue; end

                for ci = 1:numel(S.OUT)
                    if isempty(S.OUT{ci}), continue; end
                    out = S.OUT{ci};
                    if ~isfield(out, 'DATA0') || ~isfield(out, 'Experiment')
                        continue
                    end

                    try
                        freqStr = string(char(out.Experiment.Frequency));
                    catch
                        continue
                    end
                    shapeCol = find(shapeFreqStr == freqStr);
                    if isempty(shapeCol) || ~ismember(shapeCol, opts.Shapes)
                        continue
                    end

                    D0 = out.DATA0;
                    if size(D0, 1) < 3, continue; end
                    x = D0(2,:)';
                    y = D0(3,:)';
                    N = numel(x);
                    if N < minSamples, continue; end

                    t = (0:N-1)' / nativeFs;

                    try
                        qual = string(char(out.Experiment.Quality));
                    catch
                        qual = "Unknown";
                    end

                    trialIdx = countTrialForPidShape(trials, pidStr, ...
                        condPrefix, shapeCol);

                    tr = buildTrial(x, y, t, nativeFs, ...
                        pidStr, condPrefix, condLabel, ...
                        shapeNames(shapeCol), shapeCol, trialIdx, 2, ...
                        sprintf('DATA0 at %d Hz, Quality=%s, Day%d', ...
                        nativeFs, qual, dayCol - 1), 0.248);
                    trials = appendTrial(trials, tr);
                end
            end
        end
    end
    if opts.Verbose, printSummary(trials, 2); end
end

%% ========================================================================
function tr = buildTrial(x, y, t, fs, pid, condPrefix, condLabel, ...
        shapeName, shapeCol, trialIdx, study, noteStr, sigmaToMM)
    tr.x              = x(:);
    tr.y              = y(:);
    tr.z              = [];
    tr.t              = t(:);
    tr.fs             = fs;
    tr.subjectID      = sprintf('S%s', pid);
    tr.trialID        = sprintf('S%s_%s_S%d_T%d', ...
        pid, condPrefix, shapeCol, trialIdx);
    tr.condition      = condLabel;
    tr.shape          = shapeName;
    tr.database       = "hickman";
    tr.notes          = noteStr;
    tr.units_verified = false;
    tr.SigmaToMM      = sigmaToMM;
    tr.group          = condPrefix;
    tr.study          = study;
end

function trials = appendTrial(trials, tr)
    if isempty(trials)
        trials = tr;
    else
        trials(end+1) = tr; %#ok<AGROW>
    end
end

function idx = countTrialForPidShape(trials, pid, condPrefix, shapeCol)
    idx = 1;
    if isempty(trials), return; end
    pattern = sprintf('S%s_%s_S%d_T', pid, condPrefix, shapeCol);
    matches = sum(startsWith(string({trials.trialID}), pattern));
    idx = matches + 1;
end

function printSummary(trials, study)
    if isempty(trials)
        fprintf('\n=== importDB_hickman_v002 (Study %d): 0 trials ===\n', study);
        return
    end
    fprintf('\n=== importDB_hickman_v002 (Study %d) ===\n', study);
    fprintf('Trials loaded: %d\n', numel(trials));
    fprintf('Sampling rate: %d Hz\n', trials(1).fs);
    subjs = unique(string({trials.subjectID}));
    conds = unique(string({trials.condition}));
    durs = arrayfun(@(tr) tr.t(end) - tr.t(1), trials);
    nSamp = arrayfun(@(tr) numel(tr.x), trials);
    fprintf('Unique subjects: %d\n', numel(subjs));
    fprintf('Conditions: %s\n', strjoin(conds, ', '));
    fprintf('Duration: %.1f +/- %.1f s (range [%.1f, %.1f])\n', ...
        mean(durs), std(durs), min(durs), max(durs));
    fprintf('Samples: %.0f +/- %.0f (range [%d, %d])\n', ...
        mean(nSamp), std(nSamp), min(nSamp), max(nSamp));
    fprintf('Position units: experiment-native (NOT mm)\n');
end
