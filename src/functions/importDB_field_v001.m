function trials = importDB_field_v001(dataRoot, opts)
% importDB_field_v001  Import ShapeMap phone CSVs (data/field). Successor-paper
% use only: NOT part of the prereg corpus and not validated against its grid.
%
% One CSV = one trial. Returns the MOVEMENT SEGMENT only (first to last change
% in position), with no edge clipping: the caller differentiates first and
% clips afterwards, so filter edge transients do not enter the fit.
%
% TIME / fs   Taken from the numeric 'timestamp' column (seconds), never assumed.
%             fs = 1/median(dt). Errors if dt is non-positive or any interval
%             departs from the median by more than MaxDtDev (no resampling).
% UNITS       x,y returned in PIXELS (Fraser convention); SigmaToMM = mm/px from
%             Screen_size_mm / Pixels in the file. Errors if the two axes disagree.
% COORDINATES y negated so y increases upward.
% FINGER LOST A trial with finger_lost ~= 0 inside the movement is rejected, logged.
%
% USAGE   trials = importDB_field_v001();   trials = importDB_field_v001(dataRoot)
% OPTIONS MinSamples (200), MaxDtDev (0.20), Verbose (true)
%
% Fraser, D.S. (2026)  v001

    arguments
        dataRoot         (1,:) char   = defaultFieldDataRoot()
        opts.MinSamples  (1,1) double = 200
        opts.MaxDtDev    (1,1) double = 0.20
        opts.Verbose     (1,1) logical = true
    end

    need  = {'Participant_ID','Phone_model','iOS','Screen_size_mm','Pixels', ...
             'Shape_ID','Trial','Frame','mouse_x','mouse_y','finger_lost','timestamp'};
    files = dir(fullfile(dataRoot, '*.csv'));
    if isempty(files)
        error('importDB_field:noFiles', 'No CSVs in %s', dataRoot);
    end

    trialCell = {};
    rejectLog = {};
    for i = 1:numel(files)
        fn = files(i).name;
        T  = readtable(fullfile(files(i).folder, fn), 'TextType', 'string');
        miss = setdiff(need, T.Properties.VariableNames);
        if ~isempty(miss)
            error('importDB_field:badColumns', '%s lacks: %s', fn, strjoin(miss, ', '));
        end

        %% scale
        dimMM = str2double(regexp(char(T.Screen_size_mm(1)), '[\d.]+', 'match'));
        dimPx = str2double(regexp(char(T.Pixels(1)),         '[\d.]+', 'match'));
        if numel(dimMM) ~= 2 || numel(dimPx) ~= 2 || any(isnan([dimMM dimPx]))
            error('importDB_field:badScale', '%s: cannot parse screen size/pixels', fn);
        end
        mmx = dimMM(1) / dimPx(1);
        mmy = dimMM(2) / dimPx(2);
        if abs(mmx - mmy) / mmx > 0.01
            error('importDB_field:anisotropic', '%s: mm/px differs between axes (%.5f vs %.5f)', ...
                fn, mmx, mmy);
        end

        %% movement segment
        xy  = [T.mouse_x, T.mouse_y];
        chg = find(any(diff(xy) ~= 0, 2));
        if isempty(chg)
            rejectLog{end+1} = sprintf('%s: never moves', fn); %#ok<AGROW>
            continue
        end
        seg = chg(1):(chg(end) + 1);
        if numel(seg) < opts.MinSamples
            rejectLog{end+1} = sprintf('%s: %d samples in movement', fn, numel(seg)); %#ok<AGROW>
            continue
        end
        if any(T.finger_lost(seg) ~= 0)
            rejectLog{end+1} = sprintf('%s: finger_lost inside movement', fn); %#ok<AGROW>
            continue
        end

        %% time
        t  = T.timestamp(seg) - T.timestamp(seg(1));
        dt = diff(t);
        if any(dt <= 0)
            error('importDB_field:nonMonotonic', '%s: timestamp not strictly increasing', fn);
        end
        dev = max(abs(dt - median(dt))) / median(dt);
        if dev > opts.MaxDtDev
            error('importDB_field:irregular', ...
                '%s: max interval deviation %.1f%% exceeds %.0f%%', fn, 100*dev, 100*opts.MaxDtDev);
        end

        shapeName = strtrim(string(T.Shape_ID(1)));
        pid       = strtrim(string(T.Participant_ID(1)));
        phone     = strtrim(string(T.Phone_model(1)));
        iosVer    = strtrim(string(T.iOS(1)));
        tr = struct();
        tr.x            = double(T.mouse_x(seg));
        tr.y            = -double(T.mouse_y(seg));
        tr.z            = [];
        tr.t            = t(:);
        tr.fs           = 1 / median(dt);
        tr.subjectID    = char(pid);
        tr.trialID      = sprintf('%s_%s_T%d', tr.subjectID, shapeName, T.Trial(1));
        tr.condition    = 'FIELD';
        tr.shape        = char(shapeName);
        tr.database     = 'field';
        tr.notes        = sprintf('%s iOS %s; file %s', phone, iosVer, fn);
        tr.units_verified = false;
        tr.SigmaToMM    = mmx;
        tr.maxDtDev     = dev;
        trialCell{end+1} = tr; %#ok<AGROW>
    end

    if isempty(trialCell)
        error('importDB_field:noTrials', 'No valid trials. Rejections: %s', ...
            strjoin(rejectLog, '; '));
    end
    trials = [trialCell{:}]';

    if opts.Verbose
        fprintf('[importDB_field] %d trials from %d files (%d rejected)\n', ...
            numel(trials), numel(files), numel(rejectLog));
        for r = 1:numel(rejectLog)
            fprintf('  rejected: %s\n', rejectLog{r});
        end
    end
end

function dr = defaultFieldDataRoot()
% src/functions/importDB_field_v001.m -> ../../data/field
    dr = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'field');
end
