function plotPaperRoadmapSchematic_v001()
% PLOTPAPERROADMAPSCHEMATIC_V001  Two-panel "Figure 0" for the paper skeleton.
%
% Panel A: schematic flow of one trial through the paper's own logic (no
%   data dependency) -- raw trajectory -> pipeline -> beta_obs, forking
%   into the precision question (SEM, Sec 2.3) and the identifiability
%   question (forward-map local invertibility, Sec 4/5), recombining into
%   loop closure -> beta_gen*.
% Panel B: the "central paradox" as a literal scatter, not a paragraph --
%   all 42 empirical (dataset x pipeline) cells plotted by SEM adequacy
%   (x-axis) against global forward-map invertibility verdict (y-axis).
%   Sourced from real project outputs, not illustrative numbers:
%     - Precision axis: perCoordinateSEM_v2_001.mat's own `sem` column,
%       snapped to each dataset's own empirical (alpha,sigma,fs) centroid
%       (Sec 2.3 values), evaluated at the betaGen grid node nearest 1/3.
%     - Identifiability axis: each dataset's own canonical v012
%       loop-closure corpus (loopClosureResults_<dataset>_all_shaped_xu_
%       v012.mat), aggregated via the SAME method compare42CellVerdict_
%       v001.m uses (rise-coverage against results.invertStatus, PASS
%       >=95%, CONDITIONAL >=90%, FAIL<90%, "no_beta_obs" trials excluded
%       from the denominator) -- Finding #160's own methodology, not
%       re-derived independently.
%
% ASSUMPTION CARRIED FROM compare42CellVerdict_v001.m, not independently
% verified here: invertStatus(pp) column order (pp=1..6) is assumed to
% match PIPELINES below (BWFD-OLS, SG-OLS, BWFD-LMLS, SG-LMLS, BWFD-IRLS,
% SG-IRLS), the same assumption that script makes. Each .mat file also
% stores its own `pipelineLabels` variable, which could be read and
% cross-checked against PIPELINES for an extra guard if this script is
% promoted beyond a schematic-figure role.
%
% Suggested placement: immediately after the Key Terms glossary, before
% Section 0's prose -- the reader's map before the EIV/identifiability
% formalism.
%
% Fraser, D.S. (2026)

    srcDir = fileparts(mfilename('fullpath'));
    cd(srcDir);

    %% CONFIG -------------------------------------------------------------
    PIPELINES = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];

    REGISTRY = struct( ...
        'name',   {'Zarandi','Cook CTRL','Cook ASD','Dhieb','Hickman PLAC','Hickman HALO','Fraser'}, ...
        'matTag', {'Zarandi','Cook_CTRL','Cook_ASD','Dhieb','Hickman_PLAC','Hickman_HALO','Fraser'}, ...
        'alpha',  {3.18,      4.77,       5.06,      2.50,   5.34,          5.42,          4.289}, ...
        'sigma',  {4.77,      8.15,       7.84,      7.50,   7.17,          7.42,          2.009}, ...
        'fs',     {100,       133,        133,       100,    133,           133,           240});
    nDS = numel(REGISTRY);
    nPP = numel(PIPELINES);

    MDC          = 0.03;
    SEM_ADEQUATE = MDC / 2.77;      % 0.0108, Sec 2.3

    %% Panel B data (1/2): precision axis, real SEM ------------------------
    fprintf('Loading perCoordinateSEM_v2_001.mat...\n');
    if ~isfile('perCoordinateSEM_v2_001.mat')
        error('plotPaperRoadmapSchematic:NoSEMFile', ...
            'FAILED PATH: perCoordinateSEM_v2_001.mat not found in %s.', srcDir);
    end
    Ssem = load('perCoordinateSEM_v2_001.mat', 'coordTable');
    T = Ssem.coordTable;
    allAlpha = sort(unique(T.alpha));
    allSigma = sort(unique(T.sigma));
    allFs    = sort(unique(T.fs));
    allBeta  = sort(unique(T.betaGen));
    [~, bgIdx] = min(abs(allBeta - 1/3));
    snapBeta = allBeta(bgIdx);      % evaluate SEM at the canonical-value grid node

    semGrid = NaN(nDS, nPP);
    for d = 1:nDS
        [~,ai] = min(abs(allAlpha - REGISTRY(d).alpha));
        [~,si] = min(abs(allSigma - REGISTRY(d).sigma));
        [~,fi] = min(abs(allFs    - REGISTRY(d).fs));
        snapA = allAlpha(ai); snapS = allSigma(si); snapF = allFs(fi);
        for p = 1:nPP
            row = T(T.alpha==snapA & T.sigma==snapS & T.fs==snapF & ...
                    T.betaGen==snapBeta & T.pipeline==PIPELINES(p), :);
            if isempty(row)
                error('plotPaperRoadmapSchematic:NoSEMRow', '%s', sprintf( ...
                    'No perCoordinateSEM_v2_001 row for %s / %s at snapped ', ...
                    '(alpha=%.3f, sigma=%.2f, fs=%d, betaGen=%.4f).', ...
                    REGISTRY(d).name, PIPELINES(p), snapA, snapS, snapF, snapBeta));
            end
            semGrid(d, p) = row.sem(1);
        end
    end
    semAdequate = semGrid < SEM_ADEQUATE;

    %% Panel B data (2/2): identifiability axis, real v012 gate -----------
    fprintf('Loading v012 loop-closure corpora (7 datasets)...\n');
    verdictGrid  = strings(nDS, nPP);
    coverageGrid = NaN(nDS, nPP);
    for d = 1:nDS
        matFile = sprintf('loopClosureResults_%s_all_shaped_xu_v012.mat', REGISTRY(d).matTag);
        if ~isfile(matFile)
            error('plotPaperRoadmapSchematic:NoLoopClosureFile', ...
                'FAILED PATH: %s not found in %s.', matFile, srcDir);
        end
        Dset = load(matFile, 'results');
        for p = 1:nPP
            statuses = arrayfun(@(r) r.invertStatus(p), Dset.results);
            valid  = statuses ~= "no_beta_obs";
            nValid = sum(valid);
            if nValid == 0
                verdictGrid(d,p) = "no-data";
                continue
            end
            cov = sum(statuses(valid) == "rise") / nValid;
            coverageGrid(d,p) = cov;
            if cov >= 0.95
                verdictGrid(d,p) = "PASS";
            elseif cov >= 0.90
                verdictGrid(d,p) = "CONDITIONAL";
            else
                verdictGrid(d,p) = "FAIL";
            end
        end
    end

    % Fail loud, not silent: cross-check against Finding #160's own split.
    nPass = sum(verdictGrid(:) == "PASS");
    nCond = sum(verdictGrid(:) == "CONDITIONAL");
    nFail = sum(verdictGrid(:) == "FAIL");
    fprintf('\nRecomputed 42-cell split: %d PASS + %d CONDITIONAL + %d FAIL  (Finding #160: 4+6+32)\n', ...
        nPass, nCond, nFail);
    if ~isequal([nPass nCond nFail], [4 6 32])
        warning('plotPaperRoadmapSchematic:VerdictMismatch', '%s', sprintf([ ...
            'Recomputed split does not match Finding #160''s reported 4+6+32. ' ...
            'Do not silently trust Panel B -- check dataset/pipeline order and ' ...
            'mat versions against compare42CellVerdict_v001.m before using this figure.']));
    end

    %% Figure ---------------------------------------------------------------
    fig = figure('Name', 'Paper roadmap schematic', 'Position', [60 60 1500 640]);

    % ---- Panel A: schematic flow (left) ----
    axA = axes(fig, 'Position', [0.02 0.06 0.55 0.86]);
    axis(axA, 'off'); hold(axA, 'on'); xlim(axA, [0 10]); ylim(axA, [0 10]);

    grey  = [0.90 0.90 0.90];
    blue  = [0.80 0.90 0.98];
    green = [0.80 0.98 0.85];
    red   = [0.97 0.85 0.85];

    boxes = struct( ...
        'id',    {'traj','pipe','obs','sem','ident','stop','faillLoud','loop','final'}, ...
        'pos',   {[0.2 7.6],[2.3 7.6],[4.6 7.6],[1.0 4.6],[5.3 4.6],[0.2 1.6],[7.6 1.6],[3.3 1.6],[3.0 0.0]}, ...
        'w',     {2.0, 2.0, 1.9, 2.1, 2.3, 2.4, 2.3, 3.0, 3.8}, ...
        'h',     {1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.4}, ...
        'lines', { ...
            {'Raw trajectory','(x,y,t)'}, ...
            {'Pipeline','\{BWFD,SG\}x\{OLS,LMLS,IRLS\}'}, ...
            {'\beta_{obs}, VGF_{obs}'}, ...
            {'PRECISION','SEM at this coord','vs MDC/2.77 (§2.3)'}, ...
            {'IDENTIFIABILITY','locally monotonic','here? (§4/§5)'}, ...
            {'report \beta_{obs}','+ SEM caveat, stop'}, ...
            {'FAIL LOUD:','no correction','possible'}, ...
            {'loop closure:','shaped\_xu per-trial','inversion (§5)'}, ...
            {'\beta_{gen}^{*}, four uncertainty layers (§2.5)'} ...
        }, ...
        'col',   {grey, grey, grey, blue, green, red, red, green, green});

    for b = 1:numel(boxes)
        boxCentre_local(axA, boxes(b));
    end

    arrow_local(axA, edgePt_local(boxes(1),'r'), edgePt_local(boxes(2),'l'));
    arrow_local(axA, edgePt_local(boxes(2),'r'), edgePt_local(boxes(3),'l'));
    arrow_local(axA, edgePt_local(boxes(3),'b'), edgePt_local(boxes(4),'t'));
    arrow_local(axA, edgePt_local(boxes(3),'b'), edgePt_local(boxes(5),'t'));
    arrow_local(axA, edgePt_local(boxes(4),'b'), edgePt_local(boxes(6),'t'));
    arrow_local(axA, edgePt_local(boxes(4),'b'), edgePt_local(boxes(8),'t'));
    arrow_local(axA, edgePt_local(boxes(5),'b'), edgePt_local(boxes(7),'t'));
    arrow_local(axA, edgePt_local(boxes(5),'b'), edgePt_local(boxes(8),'t'));
    arrow_local(axA, edgePt_local(boxes(8),'b'), edgePt_local(boxes(9),'t'));

    title(axA, 'Panel A -- one trial''s path through the paper''s own logic', 'FontSize', 11);

    % ---- Panel B: orthogonality scatter (right) ----
    axB = axes(fig, 'Position', [0.63 0.12 0.34 0.76]); hold(axB, 'on');

    vCode = zeros(nDS, nPP);
    vCode(verdictGrid == "FAIL")        = 1;
    vCode(verdictGrid == "CONDITIONAL") = 2;
    vCode(verdictGrid == "PASS")        = 3;

    rng(20260817, 'twister');   % reproducible jitter only, not analysis
    xJit = double(semAdequate) + (rand(nDS,nPP) - 0.5) * 0.14;
    yJit = vCode + (rand(nDS,nPP) - 0.5) * 0.22;

    cols = [0.75 0.20 0.20; 0.85 0.60 0.10; 0.15 0.55 0.20];  % FAIL/CONDITIONAL/PASS
    for d = 1:nDS
        for p = 1:nPP
            c = cols(vCode(d,p), :);
            plot(axB, xJit(d,p), yJit(d,p), 'o', 'MarkerFaceColor', c, ...
                 'MarkerEdgeColor', c*0.6, 'MarkerSize', 7);
        end
    end

    xlim(axB, [-0.4 1.4]); ylim(axB, [0.5 3.5]);
    xticks(axB, [0 1]); xticklabels(axB, {'SEM inadequate','SEM adequate'});
    yticks(axB, [1 2 3]); yticklabels(axB, {'FAIL','CONDITIONAL','PASS'});
    ylabel(axB, 'Global invertibility verdict (§4, Finding #160)');
    grid(axB, 'on'); box(axB, 'on');

    nTopRight = sum(semAdequate(:) & vCode(:) == 3);
    title(axB, {'Panel B -- the central paradox, plotted', ...
        sprintf('all 42 dataset x pipeline cells; top-right = %d', nTopRight)}, ...
        'FontSize', 11);

    %% Save
    outFile = fullfile('..', 'figures', 'paperRoadmap_v001.png');
    saveas(fig, outFile);
    fprintf('\nSaved: %s\n', outFile);
    fprintf('Top-right (SEM adequate AND global PASS) cell count: %d / %d\n', ...
        nTopRight, nDS*nPP);

end

%% ========================  local functions  ============================

function boxCentre_local(ax, b)
    rectangle(ax, 'Position', [b.pos b.w b.h], 'Curvature', 0.08, ...
        'FaceColor', b.col, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.1);
    text(ax, b.pos(1) + b.w/2, b.pos(2) + b.h/2, strjoin(b.lines, '\n'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 8.3, 'Interpreter', 'tex');
end

function pt = edgePt_local(b, side)
% Midpoint of the requested edge ('l','r','t','b') of box b, in axes data coords.
    switch side
        case 'l', pt = [b.pos(1),         b.pos(2) + b.h/2];
        case 'r', pt = [b.pos(1) + b.w,    b.pos(2) + b.h/2];
        case 't', pt = [b.pos(1) + b.w/2,  b.pos(2) + b.h];
        case 'b', pt = [b.pos(1) + b.w/2,  b.pos(2)];
        otherwise
            error('edgePt_local:BadSide', 'Unknown side ''%s''.', side);
    end
end

function arrow_local(ax, xy1, xy2)
    quiver(ax, xy1(1), xy1(2), xy2(1)-xy1(1), xy2(2)-xy1(2), 0, ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 1.3, 'MaxHeadSize', 1.4);
end
