function plotVerdictHeatmapNvsN_v001()
% PLOTVERDICTHEATMAPNVSN_V001  N=20 vs N=200 verdict heatmap, continuous shading.
%
% Companion to plotEmpiricalCoverageHeatmap_v001.m, addressing a specific
% complaint about it: the 3-colour discrete PASS/CONDITIONAL/FAIL fill
% hides real variation *within* a verdict, and FAIL is by far the biggest
% group (32/42 cells) -- Dhieb/BWFD-LMLS at 5.7% and Hickman PLAC/SG-LMLS
% at 84.9% are both "FAIL" but are not remotely the same result.
%
% This version colours every cell by its actual coverage fraction using
% a continuous gradient WITHIN each gate band (deep-to-pale red across
% [0, 0.90), deep-to-pale amber across [0.90, 0.95), pale-to-deep green
% across [0.95, 1.0]) -- so the PASS/90%/CONDITIONAL/95%/FAIL boundaries
% remain visible as colour-family jumps, but the nuance inside each band
% is no longer thrown away.
%
% Reuses compare42CellVerdict_v001.m's own output table directly (single
% source of truth for the N=20 vs N=200 numbers) rather than reloading or
% re-deriving anything -- that function already does the loading/
% aggregation and is already validated (Working Notes item 10).
%
% Panels: (1) N_REPS=20 (v012, Finding #160's citable table), (2)
% N_REPS=200 (v013/v014 confirmatory corpus), with the four cells whose
% verdict changes outlined in gold, (3) a compact text summary of exactly
% those transitions.
%
% Fraser, D.S. (2026)

    srcDir = fileparts(mfilename('fullpath'));
    cd(srcDir);

    fprintf('Calling compare42CellVerdict_v001()...\n');
    T = compare42CellVerdict_v001();   % 42-row table, already validated

    PIPELINES_NATIVE  = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    PIPELINES_DISPLAY = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
    [foundMask, dispIdx] = ismember(PIPELINES_DISPLAY, PIPELINES_NATIVE);
    if ~all(foundMask)
        error('plotVerdictHeatmapNvsN:LabelMismatch', ...
            'PIPELINES_DISPLAY entries not found in PIPELINES_NATIVE.');
    end
    if height(T) ~= 42
        error('plotVerdictHeatmapNvsN:UnexpectedRowCount', ...
            'Expected 42 rows from compare42CellVerdict_v001, got %d.', height(T));
    end

    datasetNames = T.dataset(1:6:end);          % one label per 6-row block
    if numel(datasetNames) ~= 7
        error('plotVerdictHeatmapNvsN:UnexpectedDatasetCount', ...
            'Expected 7 datasets, got %d.', numel(datasetNames));
    end

    cov20  = reshape(T.covN20,  6, 7)';   % rows = dataset, cols = pipeline (native order)
    cov200 = reshape(T.covN200, 6, 7)';
    changed = reshape(T.changed, 6, 7)';

    cov20  = cov20(:,  dispIdx);          % reorder to paper's display convention
    cov200 = cov200(:, dispIdx);
    changed = changed(:, dispIdx);

    nDS = numel(datasetNames);
    nPP = numel(PIPELINES_DISPLAY);

    %% Figure ---------------------------------------------------------------
    fig = figure('Name', 'Verdict heatmap: N=20 vs N=200', 'Position', [60 60 1400 620]);

    axL = axes(fig, 'Position', [0.10 0.26 0.37 0.62]);
    drawGradientHeatmap_local(axL, cov20, datasetNames, PIPELINES_DISPLAY, false(nDS,nPP));
    title(axL, 'N_{REPS}=20 (v012, Finding #160)', 'FontSize', 11);

    axR = axes(fig, 'Position', [0.52 0.26 0.37 0.62]);
    drawGradientHeatmap_local(axR, cov200, datasetNames, PIPELINES_DISPLAY, changed);
    title(axR, 'N_{REPS}=200 (v013/v014 confirmatory) -- gold outline = verdict changed', 'FontSize', 11);

    % Shared gradient legend strip along the bottom
    axLeg = axes(fig, 'Position', [0.10 0.055 0.79 0.045]);
    drawGradientLegend_local(axLeg);

    % Compact transition summary, right-hand margin
    axSum = axes(fig, 'Position', [0.91 0.26 0.08 0.62]);
    axis(axSum, 'off'); hold(axSum, 'on'); xlim(axSum,[0 1]); ylim(axSum,[0 1]);
    text(axSum, 0, 1, 'Changed cells', 'FontWeight', 'bold', 'FontSize', 9, ...
        'VerticalAlignment', 'top');

    changeRows = find(T.changed);
    yPos = 0.85;
    for r = changeRows(:)'
        txt = sprintf('%s\n%s\n%.1f%%->%.1f%%\n%s->%s', ...
            T.dataset(r), T.pipeline(r), 100*T.covN20(r), 100*T.covN200(r), ...
            T.verdictN20(r), T.verdictN200(r));
        text(axSum, 0, yPos, txt, 'FontSize', 7.5, 'VerticalAlignment', 'top', ...
            'Interpreter', 'none');
        yPos = yPos - 0.22;
    end
    if isempty(changeRows)
        text(axSum, 0, 0.85, '(none)', 'FontSize', 8, 'VerticalAlignment', 'top');
    end

    %% Save
    outFile = fullfile('..', 'figures', 'verdictHeatmapNvsN_v001.png');
    saveas(fig, outFile);
    fprintf('\nSaved: %s\n', outFile);
    fprintf('Cells with changed verdict: %d\n', numel(changeRows));

end

%% ========================  local functions  ============================

function rgb = coverageToColor_local(cov)
% Continuous colour WITHIN each gate band; colour family jumps at the
% 90% and 95% gate thresholds themselves, so the bands stay legible.
    if isnan(cov)
        rgb = [0.82 0.82 0.82];
        return
    end
    cov = max(0, min(1, cov));
    if cov < 0.90
        t  = cov / 0.90;
        c0 = [0.45 0.05 0.05];    % deep red,  cov=0
        c1 = [0.92 0.58 0.40];    % pale red,  cov->0.90-
    elseif cov < 0.95
        t  = (cov - 0.90) / 0.05;
        c0 = [0.88 0.58 0.08];    % amber,      cov=0.90
        c1 = [0.97 0.86 0.40];    % pale amber, cov->0.95-
    else
        t  = (cov - 0.95) / 0.05;
        c0 = [0.62 0.84 0.46];    % pale green, cov=0.95
        c1 = [0.05 0.42 0.14];    % deep green, cov=1.0
    end
    rgb = c0 + t * (c1 - c0);
end

function drawGradientHeatmap_local(ax, covGrid, rowLabels, colLabels, highlightMask)
    [nDS, nPP] = size(covGrid);
    axis(ax, 'ij'); hold(ax, 'on');
    xlim(ax, [0.5 nPP+0.5]); ylim(ax, [0.5 nDS+0.5]);
    ax.XTick = 1:nPP; ax.XTickLabel = colLabels; ax.XTickLabelRotation = 30;
    ax.YTick = 1:nDS; ax.YTickLabel = rowLabels;
    box(ax, 'on');

    for d = 1:nDS
        for p = 1:nPP
            c = coverageToColor_local(covGrid(d,p));
            rectangle(ax, 'Position', [p-0.5, d-0.5, 1, 1], ...
                'FaceColor', c, 'EdgeColor', [1 1 1], 'LineWidth', 0.5);
            if isnan(covGrid(d,p))
                lbl = 'n/a';
            else
                lbl = sprintf('%.1f%%', 100*covGrid(d,p));
            end
            % Light text on dark cells, dark text on pale cells (mean
            % luminance threshold), so labels stay legible across the
            % whole gradient rather than being fixed white or black.
            lum = 0.299*c(1) + 0.587*c(2) + 0.114*c(3);
            if lum < 0.55
                txtCol = [1 1 1];
            else
                txtCol = [0.1 0.1 0.1];
            end
            text(ax, p, d, lbl, 'HorizontalAlignment', 'center', ...
                'Color', txtCol, 'FontSize', 8, 'FontWeight', 'bold');
            if highlightMask(d,p)
                rectangle(ax, 'Position', [p-0.5, d-0.5, 1, 1], ...
                    'FaceColor', 'none', 'EdgeColor', [0.85 0.65 0.05], ...
                    'LineWidth', 3);
            end
        end
    end
end

function drawGradientLegend_local(ax)
    axis(ax, 'off'); hold(ax, 'on');
    xlim(ax, [0 1]); ylim(ax, [0 1]);
    nSeg = 200;
    xs = linspace(0, 1, nSeg+1);
    for i = 1:nSeg
        c = coverageToColor_local((xs(i)+xs(i+1))/2);
        rectangle(ax, 'Position', [xs(i) 0 (xs(i+1)-xs(i)) 1], ...
            'FaceColor', c, 'EdgeColor', 'none');
    end
    for thresh = [0 0.90 0.95 1.0]
        plot(ax, [thresh thresh], [0 1], 'k-', 'LineWidth', 1);
        text(ax, thresh, -0.6, sprintf('%.0f%%', 100*thresh), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    title(ax, 'rising-branch coverage (continuous within FAIL / CONDITIONAL / PASS bands)', ...
        'FontSize', 8.5, 'FontWeight', 'normal');
end
