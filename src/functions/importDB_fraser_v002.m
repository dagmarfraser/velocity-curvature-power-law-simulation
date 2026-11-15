function [trials, rejects, info] = importDB_fraser_v002(dataRoot, opts)
% importDB_fraser_v002  Import the Fraser N=48 Ellipse corpus (Analysis_2026).
%
% v002 replaces v001 for the N=48 rebuild (Step 1 audit, 2026-11-03). v001 is
% frozen: it produced the v015 corpus and must stay reproducible. Trial
% struct fields are the same as v001 plus the provenance fields listed under
% OUTPUT, so downstream code needs no change.
%
% WHAT CHANGED FROM v001, AND WHY
%   1. Units. x,y are UIKit POINTS, not pixels: 5.209 points per mm on the
%      M4 13" iPad Pro (264.615 ppi / 25.4 = 10.418 px/mm, at 2x scale). v001
%      used 10.41793, about 2x too large a divisor, so every mm quantity
%      (sigma, VGF, velocity) it implied was about 2x too small. Beta and alpha
%      are scale-invariant. The constant is returned in .ptPerMm and
%      .SigmaToMM; callers must read it from the trial, not hard-code it.
%   2. Cohort. ID101, ID102 and ID214 are excluded by default (incomplete,
%      replaced by ID401-403; data/fraser/automatic/README.md). v001 imported
%      every folder, so the v015 "48 participants" were 45 complete plus
%      these three. Every included participant must have exactly 6 sessions
%      in 2 visits (visit boundary = gap > 90 min between session starts, as
%      the README specifies); anything else is an error, not a skip.
%   3. Direct touches. Only touch_type == 2 (Apple Pencil) samples are kept.
%      touch_type 0 (finger or palm) samples were interleaved with the Pencil
%      stream in 159 cohort trials, giving negative time steps and 50-120 mm
%      jumps. The dropped count is recorded per trial (.nDirectTouchDropped).
%   4. Pen-lift gaps (Dagmar, 2026-11-03: split and keep the longest piece).
%      After the edge clip, any step in time_seconds > GapThresh (20 ms,
%      about 5 nominal 4.2 ms steps) splits the trial; the longest gap-free
%      run is kept, and EdgeClip samples are trimmed from each run end that
%      adjoins a gap (the pen lifting and landing). Trials with no gap are
%      processed exactly as v001 processed them (bit-identical x, y, t;
%      checked 2026-11-03 on 2676 trials).
%      Threshold choice (2026-11-03): at 50 ms, 42 trials split and 14 more
%      with 20-50 ms gaps failed the clock check (item 5); at 20 ms, 51
%      split, 2 fail the clock check (slow drift, 12 ms), kept runs >= 718
%      samples. 20 ms keeps those 12 trials without loosening the check.
%   5. Clock check. Within the kept run the clock stays uniform,
%      t = (0:N-1)'/240, because time_seconds is rounded to 1 ms (steps of
%      4 and 5 ms) and differentiateKinematicsEBR's filters assume uniform
%      fs. v002 verifies it: the run's time_seconds must be strictly
%      increasing and stay within ClockTol of the uniform clock, otherwise
%      the trial is rejected with a reason.
%   6. Fail loud. Every dropped trial is a row in the second output
%      (rejects), with a reason code. Nothing is skipped silently, whatever
%      the Verbose setting. A missing file is an error.
%
% SHAPE CODES (imageNumber): 1 phi=0.4, 2 phi=0.8, 3 Ellipse, 4 Square,
% 5 Cursive. Only Ellipse (3) is supported: the gap rule is defined for a
% single continuous stroke, and Cursive has pen lifts by design.
%
% USAGE
%   [trials, rejects, info] = importDB_fraser_v002();
%   [trials, rejects, info] = importDB_fraser_v002(dataRoot, Verbose=false);
%
% OUTPUT
%   trials   [M x 1] struct, v001 fields plus:
%            .ptPerMm, .nDirectTouchDropped, .nGaps, .maxGap_s,
%            .keptFrac (kept run / clipped trial), .clockMaxDev_s, .maxStep_mm
%   rejects  table: participant, stem, trial, reason, detail
%   info     struct: constants, counts by reason, participant list
%
% See also: importDB_fraser_v001 (frozen), plotFraserEllipseGaps_v001,
%           differentiateKinematicsEBR, regressDataEBR
%
% Fraser, D.S. (2026)  v002

    arguments
        dataRoot         (1,:) char   = defaultFraserDataRoot()
        opts.Exclude     (1,:) string = ["ID101","ID102","ID214"]
        opts.Visits      (1,:) double = [1 2]
        opts.Sessions    (1,:) double = 1:3
        opts.MinSamples  (1,1) double = 200
        opts.EdgeClip    (1,1) double = 50
        opts.GapThresh   (1,1) double = 0.02    % s (about 5 nominal steps); see header item 4
        opts.ClockTol    (1,1) double = 0.010   % s, max |time_seconds - uniform clock|
        opts.VisitGapMin (1,1) double = 90      % min between session starts
        opts.MaxSubjects (1,1) double = Inf
        opts.Verbose     (1,1) logical = true
    end

    %% Constants
    FS         = 240;
    PT_PER_MM  = 5.209;          % UIKit points per mm, M4 13" iPad Pro -- see header
    SHAPE_CODE = 3;              % Ellipse only
    NU         = 2;              % Huh & Sejnowski (2015) nu for the ellipse
    PENCIL     = 2;              % UITouch.TouchType.pencil

    dataDir = fullfile(dataRoot, 'automatic');
    if ~isfolder(dataDir)
        error('importDB_fraser_v002:noDir', '%s', ['automatic/ not found in: ' dataRoot]);
    end

    %% Participants
    d = dir(fullfile(dataDir, 'ID*'));
    d = d([d.isdir]);
    allIDs = sort(string({d.name}));
    missingExcl = setdiff(opts.Exclude, allIDs);
    if ~isempty(missingExcl)
        error('importDB_fraser_v002:excludeNotFound', '%s', ...
            ['Excluded IDs not present (check the list): ' strjoin(missingExcl, ', ')]);
    end
    ids = allIDs(~ismember(allIDs, opts.Exclude));
    if opts.MaxSubjects < numel(ids)
        ids = ids(1:opts.MaxSubjects);
    end

    trialCell = {};
    rej = cell(0, 5);
    nEllipse = 0;

    for folderName = ids
        subjNum = str2double(extractAfter(folderName, "ID"));
        subjID  = sprintf('F%d', subjNum);
        subjDir = fullfile(dataDir, folderName);

        %% Sessions and visits (by timestamp gap, README rule)
        mf = dir(fullfile(subjDir, 'betaTestv002 * tracingSection.csv'));
        tok = regexp({mf.name}, '^betaTestv002 (\d{4}-\d{2}-\d{2} \d{6}) tracingSection\.csv$', ...
            'tokens', 'once');
        if any(cellfun(@isempty, tok))
            error('importDB_fraser_v002:badName', '%s', ...
                sprintf('%s: unparsable session file name(s)', folderName));
        end
        stems = sort(string(cellfun(@(c) c{1}, tok, 'UniformOutput', false)));
        if numel(stems) ~= 6
            error('importDB_fraser_v002:sessionCount', '%s', ...
                sprintf('%s has %d sessions, expected 6. Exclude it or fix the data.', ...
                folderName, numel(stems)));
        end
        tStart = datetime(stems, 'InputFormat', 'yyyy-MM-dd HHmmss');
        newVisit = [true; minutes(diff(tStart(:))) > opts.VisitGapMin];
        visitNum = cumsum(newVisit);
        if max(visitNum) ~= 2 || any(accumarray(visitNum, 1) ~= 3)
            error('importDB_fraser_v002:visitSplit', '%s', ...
                sprintf('%s: sessions do not split 3+3 into 2 visits by the %d-min rule', ...
                folderName, opts.VisitGapMin));
        end
        sessionInVisit = zeros(6, 1);
        for v = 1:2
            sessionInVisit(visitNum == v) = 1:3;
        end

        for iSess = 1:6
            vNum = visitNum(iSess); sivNum = sessionInVisit(iSess);
            if ~ismember(vNum, opts.Visits) || ~ismember(sivNum, opts.Sessions)
                continue
            end
            stem = stems(iSess);
            metaFile = fullfile(subjDir, "betaTestv002 " + stem + " tracingSection.csv");
            trajFile = fullfile(subjDir, "betaTestv002 " + stem + " tracingSection-enhanced.csv");
            if ~isfile(trajFile)
                error('importDB_fraser_v002:noEnhanced', '%s', ['Missing: ' char(trajFile)]);
            end
            metaTbl = readtable(metaFile);
            trajTbl = readtable(trajFile);
            needCols = ["response_index","touch_index","time_seconds","x","y","touch_type"];
            if ~all(ismember(needCols, string(trajTbl.Properties.VariableNames)))
                error('importDB_fraser_v002:columns', '%s', ...
                    sprintf('%s: enhanced CSV lacks a required column', trajFile));
            end

            trialNums = metaTbl.trial(metaTbl.tracingScene_object1_imageNumber == SHAPE_CODE);

            for trNum = trialNums'
                nEllipse = nEllipse + 1;
                addRej = @(reason, detail) {folderName, stem, trNum, string(reason), string(detail)};

                T = trajTbl(trajTbl.response_index == trNum, :);
                if isempty(T)
                    rej(end+1, :) = addRej("noRows", "no samples in enhanced CSV"); %#ok<AGROW>
                    continue
                end

                % Direct touches out, then order by the Pencil sample counter
                nDirect = sum(T.touch_type ~= PENCIL);
                T = sortrows(T(T.touch_type == PENCIL, :), 'touch_index');
                N = height(T);

                % Edge clip, as v001
                clip = opts.EdgeClip;
                if N <= 2*clip
                    rej(end+1, :) = addRej("tooShortBeforeClip", sprintf('%d samples', N)); %#ok<AGROW>
                    continue
                end
                T = T(clip+1 : N-clip, :);
                nClipped = height(T);

                % Pen-lift gaps: split, keep the longest run
                dt = diff(T.time_seconds);
                iGap = find(dt > opts.GapThresh);
                runStart = [1; iGap + 1];
                runEnd   = [iGap; nClipped];
                [~, j] = max(runEnd - runStart + 1);
                a = runStart(j); b = runEnd(j);
                if j > 1,              a = a + clip; end   % run begins at a landing
                if j < numel(runStart), b = b - clip; end  % run ends at a lift
                if b - a + 1 < opts.MinSamples
                    rej(end+1, :) = addRej("tooShortAfterClip", ...
                        sprintf('%d samples kept (%d gaps)', max(b-a+1,0), numel(iGap))); %#ok<AGROW>
                    continue
                end
                T = T(a:b, :);
                N = height(T);

                % Clock check against the uniform 240 Hz clock
                ts = T.time_seconds;
                if any(diff(ts) < 0)
                    rej(end+1, :) = addRej("timeNotMonotonic", "negative step after Pencil-only filter"); %#ok<AGROW>
                    continue
                end
                tVec = (0:N-1)' / FS;
                dev = (ts - ts(1)) - tVec;
                clockDev = max(abs(dev - median(dev)));
                if clockDev > opts.ClockTol
                    rej(end+1, :) = addRej("clockDeviation", ...
                        sprintf('max %.4f s from uniform 240 Hz', clockDev)); %#ok<AGROW>
                    continue
                end

                x = double(T.x);
                y = -double(T.y);       % iOS y down -> y up
                maxStep = max(hypot(diff(x), diff(y))) / PT_PER_MM;

                trial = struct();
                trial.x              = x;
                trial.y              = y;
                trial.z              = [];
                trial.t              = tVec;
                trial.fs             = FS;
                trial.subjectID      = subjID;
                trial.trialID        = sprintf('%s_V%d_S%d_T%03d_sh%d', subjID, vNum, sivNum, trNum, SHAPE_CODE);
                trial.condition      = sprintf('VISIT%d', vNum);
                trial.shape          = sprintf('nu%.4f_img%d', NU, SHAPE_CODE);
                trial.database       = 'fraser';
                trial.notes          = sprintf('stem=%s visit=%d sessionInVisit=%d sessionOverall=%d importer=v002', ...
                                           stem, vNum, sivNum, iSess);
                trial.units_verified = true;
                trial.SigmaToMM      = 1 / PT_PER_MM;
                trial.ptPerMm        = PT_PER_MM;
                trial.visit          = vNum;
                trial.sessionInVisit = sivNum;
                trial.sessionOverall = iSess;
                trial.shapeCode      = SHAPE_CODE;
                trial.nu             = NU;
                trial.edgeClip       = clip;
                trial.nDirectTouchDropped = nDirect;
                trial.nGaps          = numel(iGap);
                trial.maxGap_s       = max([dt; 0]);
                trial.keptFrac       = N / nClipped;
                trial.clockMaxDev_s  = clockDev;
                trial.maxStep_mm     = maxStep;
                trialCell{end+1} = trial; %#ok<AGROW>
            end
        end
    end

    rejects = cell2table(rej, 'VariableNames', {'participant','stem','trial','reason','detail'});
    if isempty(trialCell)
        error('importDB_fraser_v002:noTrials', '%s', 'No trials imported.');
    end
    trials = [trialCell{:}]';

    info = struct();
    info.ptPerMm = PT_PER_MM; info.fs = FS; info.excluded = opts.Exclude;
    info.participants = ids; info.nEllipse = nEllipse; info.nKept = numel(trials);
    info.nRejected = height(rejects);
    info.nSplit = sum([trials.nGaps] > 0);
    info.nWithDirectTouch = sum([trials.nDirectTouchDropped] > 0);
    if height(rejects) > 0
        [g, r] = findgroups(rejects.reason);
        info.rejectsByReason = table(r, splitapply(@numel, rejects.trial, g), ...
            'VariableNames', {'reason','n'});
    else
        info.rejectsByReason = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'reason','n'});
    end

    if opts.Verbose
        fprintf('[importDB_fraser_v002] %d participants, %d Ellipse trials: %d kept, %d rejected.\n', ...
            numel(ids), nEllipse, numel(trials), height(rejects));
        fprintf('  %.3f points/mm. Split at a gap: %d. Direct-touch samples dropped in: %d trials.\n', ...
            PT_PER_MM, info.nSplit, info.nWithDirectTouch);
        if height(rejects) > 0, disp(info.rejectsByReason); end
    end
    if height(rejects) > 0
        warning('importDB_fraser_v002:rejected', '%d of %d Ellipse trials rejected; see the rejects output.', ...
            height(rejects), nEllipse);
    end
end

function dr = defaultFraserDataRoot()
    thisDir = fileparts(mfilename('fullpath'));
    dr = fullfile(fileparts(fileparts(thisDir)), 'data', 'fraser');
end
