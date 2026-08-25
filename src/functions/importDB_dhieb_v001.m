function trials = importDB_dhieb_v001(opts)
% importDB_dhieb_v001 Import Dhieb et al. (2022) RETED ellipse data.
%
% Returns trial struct array conforming to the Stage 5 standard interface.
% Each .xls file contains one 60-second recording of continuous ellipse
% tracing on a GENIUS MousePen i608X tablet at 100 Hz.
%
% CRITICAL CAVEAT: Dhieb's data has been pre-filtered by Windows
% "Enhance Mouse Precision" (pointer precision smoothing, confirmed by
% author correspondence; the Chebyshev type 2 processing described in
% the paper refers to a separate analysis step, not the raw acquisition
% filter). The effective lowpass cutoff is approximately 10 Hz. The
% top-octave instrument noise estimate from characteriseNoise_v002 will
% reflect this rolloff, NOT the true noise floor. Only Layer 2 (biological
% template subtraction) gives meaningful noise estimates for this dataset.
%
% UNIT CONVERSION (empirically derived 2026-04-18):
%   Physical ellipse (Dhieb et al. 2022, p.787): long axis = 60 mm,
%   small axis = 36 mm, inclination = 30 deg.
%   Semi-axes: a = 30 mm, b = 18 mm.
%   Signal amplitude in data units measured from 90 trials:
%     x-channel median half-range = 186 du -> 27.5 mm expected -> 0.1478 mm/du
%     y-channel median half-range = 168 du -> 21.6 mm expected -> 0.1291 mm/du
%   Adopted conversion (x-channel, long-axis more robustly estimated):
%     SigmaToMM = 0.1478 mm/du
%   Previous placeholder (25.4/2000 = 0.0127 mm/du for generic 2000 LPI)
%   underestimated by ~11x. Corrected sigma ~7 mm (consistent with Cook
%   and Hickman at ~7-8 mm).
%
% Data: x,y coordinates only (no time column). Columns are unnamed.
% Leading/trailing rows of (0,0) are pen-up markers and are trimmed.
% Units: display/application pixels from GENIUS MousePen i608X driver
%   (NOT raw tablet LPI counts). Physical conversion: 0.1478 mm/du.
%
% Dhieb et al. (2022) DOI: 10.1080/08839514.2022.2111484
% Fraser, D.S. (2025)
% Unit conversion added: d.s.fraser@bham.ac.uk, 2026-04-18

    arguments
        opts.DataRoot   (1,1) string = defaultDhiebRoot()
        opts.Verbose    (1,1) logical = true
        opts.MinSamples (1,1) double = 200
    end

    FS          = 100;    % from paper: "sampled at 100 Hz"
    DATABASE    = "dhieb";
    SIGMA_TO_MM = 0.1478; % mm per data unit (derived from known ellipse dimensions)

    % Also load metadata for age/gender/handedness
    metaFile = fullfile(fileparts(opts.DataRoot), ...
        "Characteristics of volunteers who participated in RETED acquisition.xlsx");
    hasMeta = isfile(metaFile);
    if hasMeta
        meta = readtable(metaFile);
        metaIDs = meta.Number;
    end

    files = dir(fullfile(opts.DataRoot, "*A.xls"));
    if isempty(files)
        error("importDB_dhieb:noFiles", "No .xls files found in %s", opts.DataRoot);
    end

    trialCell = {};

    for fIdx = 1:numel(files)
        fname = files(fIdx).name;
        fpath = fullfile(opts.DataRoot, fname);

        % Extract subject number from filename (e.g. "10A.xls" -> 10)
        tokens = regexp(fname, "^(\d+)A\.xls$", "tokens");
        if isempty(tokens), continue; end
        subjNum = str2double(tokens{1}{1});

        % Read data
        T = readtable(fpath, "Sheet", 1);
        if width(T) < 2
            if opts.Verbose, fprintf("  SKIP %s: < 2 columns\n", fname); end
            continue;
        end

        x = T{:,1};
        y = T{:,2};

        % Trim leading/trailing pen-up (0,0) rows
        penDown = ~(x == 0 & y == 0);
        first   = find(penDown, 1, 'first');
        last    = find(penDown, 1, 'last');
        if isempty(first) || (last - first + 1) < opts.MinSamples
            if opts.Verbose
                fprintf("  SKIP %s: too few valid samples\n", fname);
            end
            continue;
        end
        x = x(first:last);
        y = y(first:last);

        % Build trial struct
        trial = struct();
        trial.x          = x(:);
        trial.y          = y(:);
        trial.fs         = FS;
        trial.subjectID  = sprintf("S%03d", subjNum);
        trial.trialID    = sprintf("S%03d_ellipse", subjNum);
        trial.shape      = "ellipse";
        trial.database   = DATABASE;
        trial.units      = "display_pixels";
        trial.SigmaToMM  = SIGMA_TO_MM;  % multiply sigma (data units) by this -> mm

        % Attach metadata if available
        if hasMeta
            mIdx = find(metaIDs == subjNum, 1);
            if ~isempty(mIdx)
                trial.age    = meta.Age(mIdx);
                trial.gender = string(meta.Gender{mIdx});
                trial.hand   = string(meta.DominantHand{mIdx});
                trial.notes  = sprintf( ...
                    "age=%d, %s, %s, prefiltered ~10Hz (Windows Enhance Mouse Precision), SigmaToMM=%.4f", ...
                    meta.Age(mIdx), meta.Gender{mIdx}, ...
                    meta.DominantHand{mIdx}, SIGMA_TO_MM);
            else
                trial.notes = sprintf( ...
                    "no metadata match; prefiltered ~10Hz (Windows Enhance Mouse Precision), SigmaToMM=%.4f", ...
                    SIGMA_TO_MM);
            end
        else
            trial.notes = sprintf( ...
                "prefiltered ~10Hz (Windows Enhance Mouse Precision), SigmaToMM=%.4f", SIGMA_TO_MM);
        end

        trialCell{end+1} = trial; %#ok<AGROW>
    end

    if isempty(trialCell)
        warning("importDB_dhieb:noTrials", "No valid trials extracted.");
        trials = struct([]);
        return;
    end

    trials = [trialCell{:}]';

    if opts.Verbose
        fprintf("[importDB_dhieb] %d trials from %d subjects (fs=%d Hz)\n", ...
            numel(trials), numel(trials), FS);
        fprintf("  WARNING: Data pre-filtered ~10 Hz (Windows Enhance Mouse Precision, confirmed by author correspondence).\n");
        fprintf("  Not a Chebyshev filter; that refers to a separate analysis step described in the paper.\n");
        fprintf("  Top-octave instrument noise estimate will be unreliable.\n");
        fprintf("  Unit conversion: %.4f mm/du (from known ellipse dimensions)\n", ...
            SIGMA_TO_MM);
        if hasMeta
            ages = [trials.age];
            fprintf("  Age range: %d-%d years (mean %.0f)\n", ...
                min(ages), max(ages), mean(ages));
        end
    end

end

function dr = defaultDhiebRoot()
    thisDir  = fileparts(mfilename("fullpath"));
    srcDir   = fileparts(thisDir);
    repoRoot = fileparts(srcDir);
    dr = fullfile(repoRoot, "data", "dhieb_RETED", "data", "data");
end
