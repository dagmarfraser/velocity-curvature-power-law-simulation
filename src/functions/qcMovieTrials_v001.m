function outFile = qcMovieTrials_v001(trials, outFile, opts)
% qcMovieTrials_v001 Generate MP4 QC movie from Stage 5 trial structs
%
% One frame per trial. VLC frame-advance with 'e' to scrub.
%
% Usage:
%   qcMovieTrials_v001(trials)
%   qcMovieTrials_v001(trials, 'myQC.mp4')
%   qcMovieTrials_v001(trials, [], FrameRate=4)
%
% Created 2026-03-17  v001
% Dagmar Scott Fraser - d.s.fraser@bham.ac.uk

    arguments
        trials     (1,:) struct
        outFile    (1,:) string = ""
        opts.FrameRate   (1,1) double = 4
        opts.FixedAxes   (1,1) logical = true
        opts.FigSize     (1,2) double = [800 600]
    end

    nTrials = numel(trials);
    if nTrials == 0
        warning("qcMovieTrials:empty", "No trials to render.");
        return
    end

    %% Default output filename
    if outFile == ""
        db = trials(1).database;
        ts = datestr(now, 'yyyymmdd_HHMMSS');
        outFile = sprintf("qc_%s_%s.mp4", db, ts);
    end

    %% Pre-compute annotation strings
    titleStrs = cell(nTrials, 1);
    for ii = 1:nTrials
        tr = trials(ii);
        nuStr = extractNuStr(tr.shape);
        qualStr = extractQualStr(tr.notes);
        titleStrs{ii} = sprintf('[%d/%d]  %s | %s | nu=%s | %s | N=%d (%.1fs)', ...
            ii, nTrials, char(tr.subjectID), char(tr.condition), ...
            nuStr, qualStr, numel(tr.x), tr.t(end));
    end

    %% Compute global axis limits
    if opts.FixedAxes
        xMins = arrayfun(@(t) min(t.x), trials);
        xMaxs = arrayfun(@(t) max(t.x), trials);
        yMins = arrayfun(@(t) min(t.y), trials);
        yMaxs = arrayfun(@(t) max(t.y), trials);
        xLim = [min(xMins) max(xMaxs)];
        yLim = [min(yMins) max(yMaxs)];
        pad = max(diff(xLim), diff(yLim)) * 0.05;
        xLim = xLim + [-pad pad];
        yLim = yLim + [-pad pad];
    end

    %% Setup figure — visible but positioned offscreen for reliable getframe
    fig = figure('Position', [2000 2000 opts.FigSize], ...
        'Color', 'w', 'Renderer', 'opengl');
    ax = axes(fig);

    %% Setup video writer
    vw = VideoWriter(char(outFile), 'MPEG-4');
    vw.FrameRate = opts.FrameRate;
    vw.Quality = 70;
    open(vw);

    fprintf('[qcMovie] Writing %d frames to %s ...\n', nTrials, outFile);
    tStart = tic;

    for ii = 1:nTrials
        cla(ax, 'reset');

        tr = trials(ii);
        x = tr.x;
        y = tr.y;

        % Simple line plot — fast and reliable
        plot(ax, x, y, '-', 'Color', [0.3 0.5 0.8], 'LineWidth', 0.8);
        hold(ax, 'on');
        plot(ax, x(1), y(1), 'go', 'MarkerSize', 8, 'LineWidth', 2);
        plot(ax, x(end), y(end), 'rs', 'MarkerSize', 8, 'LineWidth', 2);
        hold(ax, 'off');

        axis(ax, 'equal');
        if opts.FixedAxes
            xlim(ax, xLim);
            ylim(ax, yLim);
        end

        title(ax, titleStrs{ii}, 'FontSize', 10, 'Interpreter', 'none');

        drawnow;
        frame = getframe(fig);
        writeVideo(vw, frame);

        if mod(ii, 100) == 0 || ii == nTrials
            elapsed = toc(tStart);
            fps = ii / elapsed;
            eta = (nTrials - ii) / fps;
            fprintf('  %d/%d frames (%.1f fps, ETA %.0fs)\n', ii, nTrials, fps, eta);
        end
    end

    close(vw);
    close(fig);
    fprintf('[qcMovie] Done: %s (%.1fs)\n', outFile, toc(tStart));

end

function s = extractNuStr(shapeStr)
    shapeStr = string(shapeStr);
    tok = extractBetween(shapeStr, "nu", "_");
    if ~isempty(tok); s = char(tok{1}); else; s = '?'; end
end

function s = extractQualStr(notesStr)
    notesStr = string(notesStr);
    tok = extractBetween(notesStr, "quality=", " ");
    if ~isempty(tok); s = char(tok{1}); else; s = '?'; end
end
