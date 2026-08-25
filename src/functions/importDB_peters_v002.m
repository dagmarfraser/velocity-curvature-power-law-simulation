function trials = importDB_peters_v002(opts)
% importDB_peters_v002  Import Peters et al. (2025) visuomotor incongruence data.
%
% PsychoPy ellipse drawing task (One by Wacom CTL-672, nominal 75 Hz,
% PsychoPy height units). Template-guided counter-clockwise ellipses.
% Three conditions: congruent, temporal delay, spatial offset.
%
% RESAMPLING NOTE (v002):
%   Raw data have miniscule frame-to-frame jitter (~1-2 ms) around the
%   true 75 Hz sample grid. Uniform output (x, y) uses LINEAR interpolation
%   only. PCHIP and resample() are deliberately avoided:
%     - PCHIP adds cubic smoothing atop an already pre-filtered signal.
%     - resample() inserts a hidden FIR anti-aliasing filter.
%   Linear interpolation over sub-ms jitter is near-trivial (< 0.2 mm
%   position error at typical pen speed ~180 mm/s). Raw irregular data are
%   retained in xRaw/yRaw/tRaw for non-uniform SG differentiation (Case 8).
%
% MOUSE SMOOTHING CAVEAT:
%   Jakub Limanowski confirmed the tablet was used as a Windows mouse via
%   PsychoPy. This exposes data to the OS "Enhance Pointer Precision" (EPP)
%   filter (non-linear, speed-dependent) BEFORE PsychoPy records positions.
%   This constitutes an unknown upstream pre-filter. If EPP was active,
%   alpha_pos estimates will be elevated (artificially reddened) relative
%   to WACOM-as-digitiser datasets (Cook 133 Hz, Hickman 133 Hz). Verify
%   with Peng Wang whether EPP was disabled in the experiment setup.
%
% EXCLUSION NOTE:
%   Peters et al. excluded stopped, wrong-direction, and back-and-forth
%   trials, plus trials > 10% ellipse deviation (17.7% rejection overall).
%   These exclusions are NOT yet applied here. CSVs may contain training
%   trials (~150 rows per participant beyond the 144 experimental trials).
%   See README_PetersDataset_v001.md for full discussion.
%
% USAGE:
%   trials = importDB_peters_v002()
%   trials = importDB_peters_v002(Conditions="congruent")
%   trials = importDB_peters_v002(Conditions=["congruent","delay"], Verbose=true)
%
% OUTPUT FIELDS (Stage 5 standard interface + Peters extras):
%   x, y          - pen position (height units) at ResampleHz, linear-interp
%   t             - time vector, seconds, zero-referenced
%   fs            - = ResampleHz (default 75)
%   xRaw, yRaw    - original pen positions at irregular timestamps
%   tRaw          - original timestamps, zero-referenced (seconds)
%   subjectID     - "VP###"
%   trialID       - "VP###_T###" (row index within CSV)
%   shape         - "ellipse"
%   database      - "peters"
%   units         - coordinate unit string
%   SigmaToMM     - 336 (HP E27 G5 27-inch screen height in mm)
%   condition     - "congruent" | "delay" | "offset"
%   delay_ms      - delay in ms; 0 for congruent/offset
%   offset_height - offset in height units; 0 for congruent/delay
%   rating        - participant's classification response string
%   notes         - human-readable trial description
%
% Peters, F., Wang, P., & Limanowski, J. (2025). doi:10.1101/2025.06.18.660374
% Fraser, D.S. (2026)

    arguments
        opts.DataRoot    (1,1) string  = defaultPetersRoot()
        opts.Verbose     (1,1) logical = true
        opts.MinSamples  (1,1) double  = 75       % 1 s at 75 Hz
        opts.Conditions  (1,:) string  = ["congruent", "delay", "offset"]
        opts.ResampleHz  (1,1) double  = 75
    end

    FS_TARGET    = opts.ResampleHz;
    SCREEN_HT_MM = 336;     % HP E27 G5 FHD 27-inch: physical screen height ~336 mm
    DATABASE     = "peters";

    % Locate movement-part CSVs; skip playback_part files
    listing = dir(fullfile(opts.DataRoot, 'vp*_STSoD_movement_part_*.csv'));
    if isempty(listing)
        error('importDB_peters:noFiles', ...
            'No vp*_STSoD_movement_part_*.csv found in: %s', char(opts.DataRoot));
    end
    [~, order] = sort({listing.name});
    listing    = listing(order);

    trialCell = {};

    for fIdx = 1:numel(listing)
        fpath = fullfile(opts.DataRoot, listing(fIdx).name);

        % Read all columns as strings to preserve PsychoPy list-array cells
        rdOpts = detectImportOptions(fpath, 'TextType', 'string');
        rdOpts = setvartype(rdOpts, rdOpts.VariableNames, 'string');
        rdOpts.VariableNamingRule = 'preserve';
        T = readtable(fpath, rdOpts);

        tok  = regexp(listing(fIdx).name, 'vp(\d+)', 'tokens', 'once');
        vpID = sprintf('VP%03d', str2double(tok{1}));

        for rIdx = 1:height(T)
            row = T(rIdx, :);

            %% -- Decode condition --
            sc = row.('stim_cond');
            so = str2double(row.('stim_overall'));
            mm = str2double(row.('ms_mm'));

            if ismissing(sc) || sc == ""
                cond = "congruent";   % lone unclassified row; treat as congruent
                so   = 0;
                mm   = 0;
            elseif sc == "spat_off" && so == 0
                cond = "congruent";
            elseif sc == "spat_off" && so > 0
                cond = "offset";
            elseif sc == "temp_delay"
                cond = "delay";
            else
                continue   % unrecognised row; skip
            end

            if ~any(opts.Conditions == cond), continue; end

            %% -- Parse pen-position list-string arrays --
            % Silently skip break/incomplete rows where position cells are absent
            if ismissing(row.('mouse.x')) || ismissing(row.('mouse.y')) || ...
                    ismissing(row.('mouse.time'))
                continue
            end
            try
                xRaw = parseListStr(row.('mouse.x'){1});
                yRaw = parseListStr(row.('mouse.y'){1});
                tRaw = parseListStr(row.('mouse.time'){1});
            catch ME
                if opts.Verbose
                    fprintf('  SKIP %s row %d (parse error): %s\n', ...
                        vpID, rIdx, ME.message);
                end
                continue
            end

            if numel(xRaw) ~= numel(tRaw) || numel(yRaw) ~= numel(tRaw) ...
                    || numel(tRaw) < 2
                if opts.Verbose
                    fprintf('  SKIP %s row %d: array length mismatch\n', vpID, rIdx);
                end
                continue
            end

            %% -- Resample to uniform FS_TARGET via LINEAR interpolation --
            % Jitter is sub-millisecond; linear is near-identity over this scale.
            % PCHIP and resample() deliberately avoided (see header note).
            tVec     = tRaw - tRaw(1);    % zero-reference
            tUniform = (0 : 1/FS_TARGET : tVec(end))';
            if numel(tUniform) < opts.MinSamples, continue; end

            x = interp1(tVec, xRaw, tUniform, 'linear', 'extrap');
            y = interp1(tVec, yRaw, tUniform, 'linear', 'extrap');

            %% -- Build trial struct (Stage 5 interface) --
            tr               = struct();
            tr.x             = x(:);
            tr.y             = y(:);
            tr.t             = tUniform(:);
            tr.fs            = FS_TARGET;
            tr.xRaw          = xRaw(:);
            tr.yRaw          = yRaw(:);
            tr.tRaw          = tVec(:);
            tr.subjectID     = vpID;
            tr.trialID       = sprintf('%s_T%03d', vpID, rIdx);
            tr.shape         = "ellipse";
            tr.database      = DATABASE;
            tr.units         = "height (PsychoPy; 1 unit = ~336 mm screen height)";
            tr.SigmaToMM     = SCREEN_HT_MM;
            tr.condition     = cond;
            tr.delay_ms      = mm;
            tr.offset_height = so;
            try
                tr.rating    = string(row.('rating.mov'){1});
            catch
                tr.rating    = missing;
            end
            tr.notes = sprintf('%s | delay=%.0f ms | offset=%.4f hgt | %s row %d', ...
                cond, mm, so, vpID, rIdx);

            trialCell{end+1} = tr; %#ok<AGROW>
        end
    end

    if isempty(trialCell)
        warning('importDB_peters:noTrials', '%s', 'No valid trials extracted.');
        trials = struct([]);
        return
    end

    trials = [trialCell{:}]';

    if opts.Verbose
        printSummary(trials);
    end
end

%% -------------------------------------------------------------------------
function x = parseListStr(s)
% parseListStr  Convert PsychoPy "[a, b, c]" string to double column vector.
    s = strtrim(char(s));
    s = regexprep(s, '^\[|\]$', '');   % strip outer brackets
    x = str2double(strsplit(strtrim(s), ','))';
end

%% -------------------------------------------------------------------------
function printSummary(trials)
    fprintf('\n=== importDB_peters_v002 ===\n');
    fprintf('Trials      : %d\n', numel(trials));
    fprintf('fs          : %d Hz (linear-interpolated to uniform grid)\n', trials(1).fs);
    uSubj = unique(string({trials.subjectID}));
    fprintf('Participants: %d\n', numel(uSubj));
    for c = ["congruent", "delay", "offset"]
        nc = sum(string({trials.condition}) == c);
        fprintf('  %-12s: %d trials\n', c, nc);
    end
    nSamp = arrayfun(@(tr) numel(tr.x), trials);
    durs  = arrayfun(@(tr) tr.t(end),   trials);
    fprintf('Duration    : %.1f +/- %.1f s (range [%.1f, %.1f])\n', ...
        mean(durs), std(durs), min(durs), max(durs));
    fprintf('Samples     : %.0f +/- %.0f (range [%d, %d])\n', ...
        mean(nSamp), std(nSamp), min(nSamp), max(nSamp));
    fprintf('SigmaToMM   : 336 (HP E27 G5 27-inch screen height)\n');
    fprintf('CAUTION     : Mouse Smoothing status unknown; EPP may pre-filter data.\n');
    fprintf('CAUTION     : Exclusion criteria not applied; may include training trials.\n');
end

%% -------------------------------------------------------------------------
function dr = defaultPetersRoot()
% defaultPetersRoot  Resolve data/peters/STSoD_rawdata/ relative to src/functions/.
    thisDir  = fileparts(mfilename('fullpath'));   % src/functions/
    srcDir   = fileparts(thisDir);                % src/
    repoRoot = fileparts(srcDir);                 % repo root
    dr = fullfile(repoRoot, 'data', 'peters', 'STSoD_rawdata');
end
