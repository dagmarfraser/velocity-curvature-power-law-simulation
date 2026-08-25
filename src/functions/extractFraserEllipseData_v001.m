function extractFraserEllipseData_v001()
% EXTRACTFRASERELLIPSEDATA_V001  Copy Ellipse-only trials from the
% Analysis_2026 raw automatic/ID*/ folders into
% PowerLawSimulationPreReg/data/fraser/automatic/.
%
% Filters tracingSection.csv and tracingSection-enhanced.csv to
% imageNumber==3 (Ellipse) rows only, preserving the exact two-file-per-
% session schema used by data/pilot/automatic/ (see importDB_pilot_v001.m).
% Source retains all 5 shapes and the full raw file set (.txt,
% -trajectories.csv) untouched; only these two CSV types, Ellipse rows
% only, are copied. All 48 raw participant folders are copied (including
% the 3 incomplete ones, ID101/ID102/ID214 per Analysis_2026/claude.md) --
% completeness filtering is an importer-time decision, not a copy-time one,
% matching how data/pilot/automatic/ was handled.
%
% Fraser, D.S. (2026) v001
%
% See also: importDB_pilot_v001

    SRC_ROOT = ['/Users/dsfraser/Dropbox/Brain2Bee/StimuliApp/MSc Analysis/' ...
                'Analysis_2026/data/raw/automatic'];
    DST_ROOT = ['/Users/dsfraser/Dropbox/Brain2Bee/PowerLawSimulationPreReg/' ...
                'data/fraser/automatic'];
    SHAPE_CODE = 3;   % Ellipse (Huh & Sejnowski nu=2), see importDB_pilot_v001

    if ~isfolder(SRC_ROOT)
        error('extractFraserEllipse:noSrc', '%s', ...
            ['Source folder not found: ' SRC_ROOT]);
    end
    if ~isfolder(DST_ROOT)
        mkdir(DST_ROOT);
    end

    idDirs = dir(fullfile(SRC_ROOT, 'ID*'));
    idDirs = idDirs([idDirs.isdir]);
    if isempty(idDirs)
        error('extractFraserEllipse:noParticipants', '%s', ...
            ['No ID* folders found under: ' SRC_ROOT]);
    end

    log = struct('id', {}, 'nSessions', {}, 'nMetaRows', {}, 'nEnhRows', {});
    totalSrcBytes = 0;
    totalDstBytes = 0;

    for i = 1:numel(idDirs)
        idName = idDirs(i).name;
        srcDir = fullfile(SRC_ROOT, idName);
        dstDir = fullfile(DST_ROOT, idName);

        metaFiles = dir(fullfile(srcDir, 'betaTestv002 * tracingSection.csv'));
        if isempty(metaFiles)
            continue
        end
        if ~isfolder(dstDir)
            mkdir(dstDir);
        end

        totalMetaRows = 0;
        totalEnhRows  = 0;

        for f = 1:numel(metaFiles)
            metaName = metaFiles(f).name;
            stem     = extractBefore(metaName, ' tracingSection.csv');
            enhName  = [stem ' tracingSection-enhanced.csv'];

            metaPath = fullfile(srcDir, metaName);
            enhPath  = fullfile(srcDir, enhName);

            if ~isfile(enhPath)
                error('extractFraserEllipse:missingEnhanced', '%s', ...
                    sprintf('Missing enhanced.csv for %s / %s', idName, stem));
            end

            meta = readtable(metaPath);
            ellipseTrials = meta.trial(meta.tracingScene_object1_imageNumber == SHAPE_CODE);
            if isempty(ellipseTrials)
                error('extractFraserEllipse:noEllipseTrials', '%s', ...
                    sprintf('No Ellipse trials found for %s / %s -- unexpected, check imageNumber coding', ...
                    idName, stem));
            end

            metaOut = meta(ismember(meta.trial, ellipseTrials), :);

            enh    = readtable(enhPath, 'TextType', 'string');
            enhOut = enh(ismember(enh.response_index, ellipseTrials), :);

            writetable(metaOut, fullfile(dstDir, metaName));
            writetable(enhOut,  fullfile(dstDir, enhName));

            totalMetaRows = totalMetaRows + height(metaOut);
            totalEnhRows  = totalEnhRows  + height(enhOut);

            srcInfo = dir(metaPath); srcBytes = srcInfo.bytes;
            srcInfo2 = dir(enhPath); srcBytes = srcBytes + srcInfo2.bytes;
            dstInfo = dir(fullfile(dstDir, metaName)); dstBytes = dstInfo.bytes;
            dstInfo2 = dir(fullfile(dstDir, enhName)); dstBytes = dstBytes + dstInfo2.bytes;
            totalSrcBytes = totalSrcBytes + srcBytes;
            totalDstBytes = totalDstBytes + dstBytes;
        end

        log(end+1) = struct('id', idName, 'nSessions', numel(metaFiles), ...
            'nMetaRows', totalMetaRows, 'nEnhRows', totalEnhRows); %#ok<AGROW>

        fprintf('%s: %d sessions, %d Ellipse trials, %d trajectory rows\n', ...
            idName, numel(metaFiles), totalMetaRows, totalEnhRows);
    end

    fprintf('\nDone. %d participant folders written to %s\n', numel(log), DST_ROOT);
    fprintf('Source (all 5 shapes) meta+enhanced bytes touched: %.1f MB\n', totalSrcBytes/1e6);
    fprintf('Dest   (Ellipse only)  meta+enhanced bytes written: %.1f MB\n', totalDstBytes/1e6);
    fprintf('Total sessions: %d, total Ellipse trials: %d, total trajectory rows: %d\n', ...
        sum([log.nSessions]), sum([log.nMetaRows]), sum([log.nEnhRows]));

end
