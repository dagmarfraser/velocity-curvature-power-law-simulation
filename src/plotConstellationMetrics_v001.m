function plotConstellationMetrics_v001(opts)
% plotConstellationMetrics_v001  Produce Fig 7 (shaped_xu validation figures).
%
% Reads constellationMetrics_v002.mat (no re-simulation) and produces three
% figures matching the skeleton FIGURES table filenames:
%
%   Fig 7a  loopClosure_shaped_xu_v007_perDataset_SGIRLS.png
%           Validity gradient: cccMed (PRIMARY) and loopCCC (SECONDARY) per
%           dataset for shaped_xu. Shows that per-trial contrast CCC is the
%           defensible absolute metric; loopCCC is the noise-model-selection
%           secondary. Grade thresholds at 0.50, 0.70, 0.90.
%
%   Fig 7b  loopClosure_shaped_xu_v007_allPipelines.png
%           Three-model comparison: cccMed for xu, shaped_xu, bootstrap.
%           Visualises the two-factor picture: xu fails completely; shaped_xu
%           and bootstrap both adequate (max gap 0.038).
%
%   Fig 7c  loopClosure_shaped_xu_v007_scatter_obsVsGenStar.png
%           Predicted vs observed inter-pipeline delta-beta scatter (shaped_xu;
%           valid trials pooled; identity line; coloured by dataset).
%           Zarandi excluded (degenerate by construction; MAE~0 noted in caption).
%
% USAGE
%   plotConstellationMetrics_v001
%   plotConstellationMetrics_v001('SaveFigs', false)   % display only
%
% Fraser, D.S. (2026)  v001

    arguments
        opts.SaveFigs     (1,1) logical = true
        opts.TiledScatter (1,1) logical = false
        opts.FigDir       (1,1) string  = ""   % default: project figures/ dir
    end

    %% -- Paths --------------------------------------------------------------
    srcDir = fileparts(mfilename("fullpath"));
    if opts.FigDir == ""
        opts.FigDir = fullfile(fileparts(srcDir), "figures");
    end
    if ~isfolder(opts.FigDir) && opts.SaveFigs
        mkdir(opts.FigDir);
    end
    matFile = fullfile(srcDir, "constellationMetrics_v002.mat");
    if ~isfile(matFile)
        error("plotConstellationMetrics_v001:missingMat", "%s", ...
            "constellationMetrics_v002.mat not found. Run constellationMetrics_v002 first.");
    end

    %% -- Load ---------------------------------------------------------------
    fprintf("Loading constellationMetrics_v002.mat ...\n");
    R = load(matFile, "results");
    r = R.results;

    sxPD  = r.byNoiseModel.shaped_xu.perDataset;
    xuPD  = r.byNoiseModel.xu.perDataset;
    boPD  = r.byNoiseModel.bootstrap.perDataset;

    % Dataset display order for bars (excluding Zarandi; Dhieb last with caveat)
    orderedDs = ["Pilot","Cook CTRL","Cook ASD","Hickman HALO","Hickman PLAC","Dhieb"];
    nDs = numel(orderedDs);

    % Helper: extract field for ordered datasets from a perDataset struct array
    getField = @(pd, fname) arrayfun(@(ds) extractField(pd, ds, fname), orderedDs);

    %% =====================================================================
    %% Fig 7a: Validity gradient (cccMed PRIMARY + loopCCC SECONDARY)
    %% =====================================================================
    fprintf("Building Fig 7a: validity gradient ...\n");

    cccMedSx  = getField(sxPD,  "cccMed");
    loopCCCSx = getField(sxPD,  "loopCCCstored");
    nValidSx  = getField(sxPD,  "Nvalid");

    % Dhieb gets a dagger marker (EPP caveat)
    labels = string(orderedDs);
    labels(end) = labels(end) + " †";

    fig1 = figure("Name", "Fig7a_ValidityGradient", "Position", [50 50 820 420]);
    ax = axes(fig1);

    xPos  = 1:nDs;
    bw    = 0.35;
    colPrimary   = [0.18 0.45 0.73];   % blue: cccMed (primary)
    colSecondary = [0.70 0.80 0.90];   % light blue: loopCCC (secondary)
    colDhieb     = [0.85 0.40 0.20];   % orange: EPP-caveated

    for i = 1:nDs
        cPri = colPrimary;
        cSec = colSecondary;
        if i == nDs   % Dhieb: muted colours
            cPri = colDhieb;
            cSec = colDhieb + 0.25;
            cSec(cSec > 1) = 1;
        end
        b1 = bar(ax, xPos(i) - bw/2, cccMedSx(i),  bw, "FaceColor", cPri, "EdgeColor", "none");
        hold(ax, "on");
        b2 = bar(ax, xPos(i) + bw/2, loopCCCSx(i), bw, "FaceColor", cSec, "EdgeColor", "k", ...
            "LineWidth", 0.8);
        if i == 1
            b1Ref = b1; b2Ref = b2;
        end
    end

    % Grade thresholds
    ylines = [0.50, 0.70, 0.90];
    for yl = ylines
        xline(ax, NaN);   % dummy; use yline-equivalent
        plot(ax, [0.3, nDs+0.7], [yl yl], "k:", "LineWidth", 0.8);
        text(ax, nDs + 0.75, yl, sprintf("%.2f", yl), ...
            "FontSize", 8, "VerticalAlignment", "middle", "Color", [0.4 0.4 0.4]);
    end

    % Zarandi annotation
    zarSx = extractField(sxPD, "Zarandi", "loopCCCstored");
    text(ax, 0.45, 0.09, sprintf("Zarandi: Degenerate\n(Scenario A; MAE≈0)\nloopCCC=%.3f", zarSx), ...
        "FontSize", 7.5, "Color", [0.5 0.5 0.5], "HorizontalAlignment", "left", ...
        "BackgroundColor", [1 1 1 0.75], "EdgeColor", [0.7 0.7 0.7]);

    % N labels above bars
    for i = 1:nDs
        text(ax, xPos(i), cccMedSx(i) + 0.02, sprintf("N=%d", nValidSx(i)), ...
            "FontSize", 7, "HorizontalAlignment", "center", "Color", [0.3 0.3 0.3]);
    end

    set(ax, "XTick", xPos, "XTickLabel", labels, "XTickLabelRotation", 15, ...
        "YLim", [0 1.05], "YTick", 0:0.2:1, "FontSize", 10, "Box", "off");
    xlabel(ax, "Dataset", "FontSize", 11);
    ylabel(ax, "CCC", "FontSize", 11);
    title(ax, "Fig 7a: Constellation CCC validity gradient (shaped\_xu)", ...
        "FontSize", 11, "FontWeight", "bold");
    legend(ax, [b1Ref, b2Ref], ...
        {"cccMed (primary: per-trial contrast CCC)", ...
         "loopCCC (secondary: pooled-across-pipelines)"}, ...
        "Location", "northwest", "FontSize", 9);
    text(ax, nDs, -0.12, "† EPP-suspected (R7 FAIL; held out of pool)", ...
        "FontSize", 8, "Color", colDhieb, "HorizontalAlignment", "right");

    saveOrShow(fig1, opts, fullfile(opts.FigDir, "loopClosure_shaped_xu_v007_perDataset_SGIRLS.png"));

    %% =====================================================================
    %% Fig 7b: Three-model comparison (cccMed for xu, shaped_xu, bootstrap)
    %% =====================================================================
    fprintf("Building Fig 7b: three-model comparison ...\n");

    cccMedXu = getField(xuPD,  "cccMed");
    cccMedBo = getField(boPD,  "cccMed");

    fig2 = figure("Name", "Fig7b_ThreeModelComparison", "Position", [50 50 820 420]);
    ax2  = axes(fig2);
    disableDefaultInteractivity(ax2);

    bwM   = 0.26;
    colXu  = [0.80 0.35 0.35];   % red: xu (fails)
    colSx  = [0.18 0.45 0.73];   % blue: shaped_xu (primary)
    colBo  = [0.25 0.65 0.35];   % green: bootstrap (robustness)

    hold(ax2, "on");
    for i = 1:nDs
        bar(ax2, xPos(i) - bwM,   cccMedXu(i),  bwM, "FaceColor", colXu,  "EdgeColor", "none");
        bar(ax2, xPos(i),          cccMedSx(i),  bwM, "FaceColor", colSx,  "EdgeColor", "none");
        bBo = bar(ax2, xPos(i) + bwM, cccMedBo(i), bwM, "FaceColor", colBo,  "EdgeColor", "none");
    end
    bXuRef = bar(ax2, 0, NaN, 0.1, "FaceColor", colXu,  "EdgeColor", "none", "Visible", "off");
    bSxRef = bar(ax2, 0, NaN, 0.1, "FaceColor", colSx,  "EdgeColor", "none", "Visible", "off");
    bBoRef = bar(ax2, 0, NaN, 0.1, "FaceColor", colBo,  "EdgeColor", "none", "Visible", "off");

    plot(ax2, [0.3, nDs+0.7], [0 0], "k-", "LineWidth", 0.5);  % zero line

    set(ax2, "XTick", xPos, "XTickLabel", labels, "XTickLabelRotation", 15, ...
        "YLim", [-0.15 1.05], "YTick", -0.1:0.2:1, "FontSize", 10, "Box", "off");
    xlabel(ax2, "Dataset", "FontSize", 11);
    ylabel(ax2, "cccMed (per-trial contrast CCC)", "FontSize", 11);
    title(ax2, "Fig 7b: Three-model comparison (xu, shaped\_xu, bootstrap)", ...
        "FontSize", 11, "FontWeight", "bold");
    legend(ax2, [bXuRef, bSxRef, bBoRef], ...
        {"xu (alpha-consistent phase, no amplitude structure; fails)", ...
         "shaped\_xu (primary: semi-parametric, alpha-auditable)", ...
         "bootstrap (robustness check: non-parametric)"}, ...
        "Location", "northwest", "FontSize", 9);

    saveOrShow(fig2, opts, fullfile(opts.FigDir, "loopClosure_shaped_xu_v007_allPipelines.png"));

    %% =====================================================================
    %% Fig 7c: Predicted vs observed delta-beta scatter (shaped_xu)
    %% =====================================================================
    fprintf("Building Fig 7c: predicted vs observed delta-beta scatter ...\n");

    % Pool valid-trial contrast pairs across all non-degenerate datasets
    dsColours = lines(nDs);  % one colour per dataset
    fig3 = figure("Name", "Fig7c_PredVsObsScatter", "Position", [50 50 500 500]);
    ax3  = axes(fig3);
    hold(ax3, "on");

    for i = 1:nDs
        pd = getPerDataset(sxPD, orderedDs(i));
        if isempty(pd) || pd.Nvalid < 2; continue; end
        dbo  = pd.deltaBetaObs;
        dbp  = pd.deltaBetaPred;
        valid = isfinite(dbo(:,1)) & isfinite(dbp(:,1));
        if sum(valid) < 2; continue; end
        x_all = dbo(valid, :);
        y_all = dbp(valid, :);
        scatter(ax3, x_all(:), y_all(:), 2, dsColours(i,:), ...
            "filled", "MarkerFaceAlpha", 0.18, "DisplayName", char(labels(i)));
    end

    % Identity line
    axLim = 0.25;
    plot(ax3, [-axLim axLim], [-axLim axLim], "k-", "LineWidth", 1.2);
    plot(ax3, [-axLim axLim], [0 0], "k:", "LineWidth", 0.6);
    plot(ax3, [0 0], [-axLim axLim], "k:", "LineWidth", 0.6);

    set(ax3, "XLim", [-axLim axLim], "YLim", [-axLim axLim], ...
        "FontSize", 10, "Box", "off", "DataAspectRatio", [1 1 1]);
    xlabel(ax3, "Observed \Delta\beta (inter-pipeline)", "FontSize", 11);
    ylabel(ax3, "Predicted \Delta\beta (shaped\_xu)", "FontSize", 11);
    title(ax3, "Fig 7c: Predicted vs observed constellation contrasts", ...
        "FontSize", 11, "FontWeight", "bold");

    legend(ax3, "show", "Location", "northwest", "FontSize", 8);
    text(ax3, axLim, -axLim + 0.01, "Zarandi excluded (degenerate; MAE\approx0)", ...
        "FontSize", 8, "HorizontalAlignment", "right", "Color", [0.5 0.5 0.5]);

    saveOrShow(fig3, opts, fullfile(opts.FigDir, "loopClosure_shaped_xu_v007_scatter_obsVsGenStar.png"));

    %% =====================================================================
    %% Fig 7d (optional): per-dataset tiled scatter — contrast-type coloured
    %% =====================================================================
    if opts.TiledScatter
        fprintf("Building Fig 7d: tiled per-dataset scatter (contrast-type colour) ...\n");

        % Contrast type for each of the 15 lower-triangle pairs (k=1..15).
        % Pipeline order: BWFD-OLS=1 SG-OLS=2 BWFD-LMLS=3 SG-LMLS=4 BWFD-IRLS=5 SG-IRLS=6
        % Derivation contrast: same regression method, different derivation (BWFD vs SG)
        % Regression contrast: same derivation, different regression (OLS / LMLS / IRLS)
        % Mixed: both differ
        % k: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        contrastType = [1, 2, 3, 3, 2, 1, 2, 3, 2, 3, 3, 2, 3, 2, 1];
        ctColors = [0.18 0.45 0.73;   % 1-Derivation: blue
                    0.85 0.33 0.10;   % 2-Regression: orange
                    0.45 0.67 0.19];  % 3-Mixed:      green
        ctNames  = {"Derivation (BWFD\leftrightarrowSG, same regression; k=1,6,15)", ...
                    "Regression (same derivation, OLS/LMLS/IRLS; k=2,5,7,9,12,14)", ...
                    "Mixed (diff derivation AND regression; k=3,4,8,10,11,13)"};

        fig4 = figure("Name", "Fig7d_TiledScatter", "Position", [50 50 1200 780]);
        tl   = tiledlayout(fig4, 2, 3, "TileSpacing", "compact", "Padding", "compact");
        title(tl, "Predicted vs observed \Delta\beta by dataset (shaped\_xu)", ...
            "FontSize", 11, "FontWeight", "bold");

        axLim    = 0.25;
        legHands = gobjects(3, 1);   % filled for all 3 types in panel 1 — no placeholders

        for i = 1:nDs
            pd  = getPerDataset(sxPD, orderedDs(i));
            axT = nexttile(tl);
            hold(axT, "on");

            if ~isempty(pd) && pd.Nvalid >= 2
                dbo   = pd.deltaBetaObs;
                dbp   = pd.deltaBetaPred;
                valid = isfinite(dbo(:,1)) & isfinite(dbp(:,1));
                nV    = sum(valid);

                if nV >= 2
                    xAll = dbo(valid, :);   % Nvalid × 15
                    yAll = dbp(valid, :);

                    % Scatter per contrast type (3 calls, 3 colours)
                    for ct = 1:3
                        mask = contrastType == ct;
                        xCT  = xAll(:, mask);
                        yCT  = yAll(:, mask);
                        h = scatter(axT, xCT(:), yCT(:), 4, ctColors(ct,:), ...
                            "filled", "MarkerFaceAlpha", 0.30);
                        if i == 1; legHands(ct) = h; end
                    end

                    % OLS regression line across all pairs
                    xv = xAll(:);  yv = yAll(:);
                    fin = isfinite(xv) & isfinite(yv);
                    p   = polyfit(xv(fin), yv(fin), 1);
                    plot(axT, [-axLim axLim], polyval(p, [-axLim axLim]), ...
                        "r--", "LineWidth", 1.4);
                    text(axT, -axLim*0.95, axLim*0.88, ...
                        sprintf("slope=%.2f", p(1)), "FontSize", 8, ...
                        "Color", [0.75 0 0], "HorizontalAlignment", "left");
                end

                Cb    = pd.cccMed / max(abs(pd.pearsonMed), 0.001);
                dsLab = char(orderedDs(i));
                if i == nDs; dsLab = [dsLab " \dagger"]; end %#ok<AGROW>
                tStr  = sprintf("%s\nccc=%.3f  r=%.3f  Cb=%.3f\nMAE=%.4f  N=%d", ...
                    dsLab, pd.cccMed, pd.pearsonMed, Cb, pd.maeMed, pd.Nvalid);
            else
                tStr = sprintf("%s\nDegenerate", char(orderedDs(i)));
            end

            % Identity line and crosshairs
            plot(axT, [-axLim axLim], [-axLim axLim], "k-",  "LineWidth", 1.2);
            plot(axT, [-axLim axLim], [0 0],          "k:",  "LineWidth", 0.5);
            plot(axT, [0 0],          [-axLim axLim], "k:",  "LineWidth", 0.5);

            set(axT, "XLim", [-axLim axLim], "YLim", [-axLim axLim], ...
                "FontSize", 9, "Box", "off", "DataAspectRatio", [1 1 1]);
            title(axT, tStr, "FontSize", 8.5, "Color", dsColours(i,:));

            % Legend on first panel only
            if i == 1 && all(isvalid(legHands))
                legend(axT, legHands, ctNames, "Location", "northwest", ...
                    "FontSize", 7, "Box", "off");
            end

            if mod(i-1, 3) == 0
                ylabel(axT, "Predicted \Delta\beta", "FontSize", 9);
            end
            if i > 3
                xlabel(axT, "Observed \Delta\beta", "FontSize", 9);
            end
        end

        saveOrShow(fig4, opts, ...
            fullfile(opts.FigDir, "loopClosure_shaped_xu_v007_scatter_tiled.png"));
    end
        fprintf("Building Fig 7d: tiled per-dataset scatter ...\n");

        fig4 = figure("Name", "Fig7d_TiledScatter", "Position", [50 50 1200 780]);
        tl   = tiledlayout(fig4, 2, 3, "TileSpacing", "compact", "Padding", "compact");
        title(tl, "Predicted vs observed \Delta\beta by dataset (shaped\_xu)", ...
            "FontSize", 11, "FontWeight", "bold");

        axLim = 0.25;

        for i = 1:nDs
            pd  = getPerDataset(sxPD, orderedDs(i));
            axT = nexttile(tl);
            hold(axT, "on");

            if ~isempty(pd) && pd.Nvalid >= 2
                dbo   = pd.deltaBetaObs;
                dbp   = pd.deltaBetaPred;
                valid = isfinite(dbo(:,1)) & isfinite(dbp(:,1));
                nV    = sum(valid);

                if nV >= 2
                    xAll = dbo(valid,:);
                    yAll = dbp(valid,:);
                    cAll = abs(xAll);       % colour by |observed delta-beta|

                    scatter(axT, xAll(:), yAll(:), 4, cAll(:), ...
                        "filled", "MarkerFaceAlpha", 0.30);
                    colormap(axT, "hot");
                    clim(axT, [0, axLim]);

                    % OLS regression line
                    xv = xAll(:);  yv = yAll(:);
                    p  = polyfit(xv(isfinite(xv) & isfinite(yv)), ...
                                 yv(isfinite(xv) & isfinite(yv)), 1);
                    xFit = [-axLim, axLim];
                    plot(axT, xFit, polyval(p, xFit), "r--", "LineWidth", 1.4);
                    text(axT, -axLim*0.95, axLim*0.88, ...
                        sprintf("slope=%.2f", p(1)), ...
                        "FontSize", 8, "Color", [0.75 0 0], ...
                        "HorizontalAlignment", "left");
                end

                % Dhieb: EPP caveat marker
                dsLabel = char(orderedDs(i));
                if i == nDs; dsLabel = [dsLabel " \dagger"]; end %#ok<AGROW>

                tStr = sprintf("%s\ncccMed=%.3f  r=%.3f  MAE=%.4f\nN_{trials}=%d", ...
                    dsLabel, pd.cccMed, pd.pearsonMed, pd.maeMed, pd.Nvalid);
            else
                tStr = sprintf("%s\nDegenerate", char(orderedDs(i)));
            end

            % Identity line and zero crosshairs
            plot(axT, [-axLim axLim], [-axLim axLim], "k-",  "LineWidth", 1.2);
            plot(axT, [-axLim axLim], [0 0],          "k:", "LineWidth", 0.5);
            plot(axT, [0 0],          [-axLim axLim], "k:", "LineWidth", 0.5);

            set(axT, "XLim", [-axLim axLim], "YLim", [-axLim axLim], ...
                "FontSize", 9, "Box", "off", "DataAspectRatio", [1 1 1]);
            title(axT, tStr, "FontSize", 8.5, "Color", dsColours(i,:));

            if mod(i-1, 3) == 0
                ylabel(axT, "Predicted \Delta\beta", "FontSize", 9);
            end
            if i > 3
                xlabel(axT, "Observed \Delta\beta", "FontSize", 9);
            end
        end

    fprintf("\nplotConstellationMetrics_v001 done.\n");
end

%% ============================ LOCAL HELPERS ================================

function val = extractField(pd, dsName, fname)
% Extract a scalar field from a perDataset struct for a named dataset.
% Returns NaN if the dataset is not found.
    val = NaN;
    for i = 1:numel(pd)
        if pd(i).dataset == dsName
            val = pd(i).(fname);
            return
        end
    end
end

function pd = getPerDataset(pdArr, dsName)
% Return the perDataset struct for a named dataset, or [] if not found.
    pd = [];
    for i = 1:numel(pdArr)
        if pdArr(i).dataset == dsName
            pd = pdArr(i);
            return
        end
    end
end

function saveOrShow(fig, opts, fpath)
% Save figure as PNG (300 dpi) or display if SaveFigs is false.
    if opts.SaveFigs
        exportgraphics(fig, fpath, "Resolution", 300);
        fprintf("  Saved: %s\n", fpath);
        close(fig);
    else
        drawnow;
    end
end
