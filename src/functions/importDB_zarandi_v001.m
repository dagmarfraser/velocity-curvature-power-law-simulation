function trials = importDB_zarandi_v001(opts)
% importDB_zarandi_v001 Import Zarandi et al. (2023) WACOM ellipse data.
%
% Returns trial struct array conforming to the Stage 5 standard interface
% (same fields as importDB_cook_v002). Extracts dominant-hand natural-tempo
% ellipse trials following the parsing logic in PowerLawEmpirical.m.
%
% Zarandi et al. (2023) DOI: 10.1038/s41598-023-34861-x
% WACOM Intuos Pro, 100 Hz. Units: cm (pending confirmation).
%
% USAGE:
%   trials = importDB_zarandi_v001()
%   trials = importDB_zarandi_v001(Verbose=true)
%
% Fraser, D.S. (2025)

    arguments
        opts.DataRoot  (1,1) string = defaultZarandiRoot()
        opts.Verbose   (1,1) logical = true
        opts.MinSamples (1,1) double = 100
    end

    FS = 100;  % 100 Hz sampling
    DATABASE = "zarandi";

    subjectFiles = arrayfun(@(s) sprintf("subject%d.txt", s), 0:13);

    trialCell = {};

    for sIdx = 1:numel(subjectFiles)
        fpath = fullfile(opts.DataRoot, subjectFiles(sIdx));
        if ~isfile(fpath)
            if opts.Verbose, fprintf("  SKIP: %s not found\n", subjectFiles(sIdx)); end
            continue;
        end

        [subjID, dominant_hand, ~, trialVec, ~, x, y, ~, ~, ~, ~] = ...
            importfileZarandi(char(fpath));

        domIdx = find(dominant_hand);
        if isempty(domIdx), continue; end

        trialDom = trialVec(domIdx);
        xDom = x(domIdx);
        yDom = y(domIdx);

        % Segment: three speed blocks separated by trial==0 runs.
        % Block 2 (between 2nd and 3rd zero-run start) = natural tempo.
        zeroStarts = strfind([0, trialDom'==0], [0 1]);
        if numel(zeroStarts) < 3
            if opts.Verbose, fprintf("  SKIP subject%d: < 3 trial blocks\n", sIdx-1); end
            continue;
        end
        block2End = zeroStarts(3);

        for trialNum = 0:9
            tStart = strfind([0, trialDom'==trialNum], [0 1]);
            tEnd   = strfind([trialDom'==trialNum, 0], [1 0]);

            if length(tStart) < 2 || tStart(2) >= block2End
                continue;
            end

            xTrial = xDom(tStart(2):tEnd(2));
            yTrial = yDom(tStart(2):tEnd(2));

            if numel(xTrial) < opts.MinSamples, continue; end

            trial = struct();
            trial.x         = xTrial(:);
            trial.y         = yTrial(:);
            trial.fs        = FS;
            trial.subjectID = sprintf("S%02d", sIdx - 1);
            trial.trialID   = sprintf("S%02d_T%d", sIdx - 1, trialNum);
            trial.shape      = "ellipse";
            trial.database   = DATABASE;
            trial.notes      = sprintf("dominant hand, natural tempo, trial %d", trialNum);
            trial.units      = "cm (unverified)";
            trial.SigmaToMM  = 10.0;              % cm → mm (Zarandi native units)

            trialCell{end+1} = trial; %#ok<AGROW>
        end
    end

    if isempty(trialCell)
        warning("importDB_zarandi:noTrials", "No valid trials extracted.");
        trials = struct([]);
        return;
    end

    trials = [trialCell{:}]';

    if opts.Verbose
        uSubj = unique(string({trials.subjectID}));
        fprintf("[importDB_zarandi] %d trials from %d subjects (fs=%d Hz)\n", ...
            numel(trials), numel(uSubj), FS);
    end

end

function dr = defaultZarandiRoot()
    thisDir  = fileparts(mfilename("fullpath"));
    srcDir   = fileparts(thisDir);
    repoRoot = fileparts(srcDir);
    dr = fullfile(repoRoot, "data", "zarandi", ...
        "IJCNN2022_Ellipses-main", "IJCNN2022_Ellipses-main", "data");
end
