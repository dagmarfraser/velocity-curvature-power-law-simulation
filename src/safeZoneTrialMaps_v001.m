%% safeZoneTrialMaps_v001.m
% P-2 (docs/TODO_v003_Rewrite_v001.md): safe zones and clearance margins at each trial's own
% tempo. v002's safe zones (Part 2; plotForwardMapAtCentroid_v001 §4) were read from one
% fixed-VGF grid cell per dataset, whose fold is a tempo effect (#224); every upper bound
% came out at the fold, 0.50. Here the SAME rule is applied to every v015 per-trial forward
% map (each built at the trial's own f0, sigma, alpha and geometry, shaped_xu noise):
%   LB   first beta_gen (from the 2nd node) with beta_rec > 0.20 and the next 3 steps rising
%   UB   beta_gen at the map's peak (argmax beta_rec); "at top" if the peak is the last node
%   zone per trial: [max over pipelines of LB, min over pipelines of UB]  (conservative, as v002)
% Per dataset: median per-trial bounds; where beta_gen* sits; clearance beta_gen* - LB
% (SG-IRLS and conservative), against MDC 0.03. v002's grid-cell bounds printed alongside.
% Zarandi is reported but outside the validity domain (D20); Dhieb has one trial per subject.
% Read-only; no simulation. Run from src/:  safeZoneTrialMaps_v001
% Inputs: src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Outputs: results/safeZoneTrialMaps_v001.mat, figures/safeZoneTrialMaps_v001.png
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
PL       = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  % runner order (v012 L118)
FOCUS    = "SG-IRLS";
LB_FLOOR = 0.20;  LB_RUN = 3;  SKIP = 1;          % plotForwardMapAtCentroid_v001 L~205-208
MDC      = 0.03;
V002     = table(["Fraser"; "Cook_CTRL"; "Cook_ASD"; "Hickman_PLAC"; "Hickman_HALO"; "Zarandi"], ...
                 [0.27; 0.23; 0.23; 0.23; 0.23; 0.37], [0.50; 0.50; 0.50; 0.50; 0.50; 0.53], ...
                 'VariableNames', ["dataset" "v002LB" "v002UB"]);   % v002 Part 2, grid cells
OUT_MAT  = fullfile(ROOT, "results", "safeZoneTrialMaps_v001.mat");
OUT_PNG  = fullfile(ROOT, "figures", "safeZoneTrialMaps_v001.png");
for f = [fileparts(OUT_MAT) fileparts(OUT_PNG)]
    if ~isfolder(f), error("szTrial:outDir", "%s", "Missing folder: " + f); end
end
iF = find(PL == FOCUS);

%% Per trial x pipeline bounds
rows = {};  curves = struct();
for d = DATASETS
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("szTrial:input", "%s", "Missing: " + f); end
    L = load(f, "results", "betaGenVec");  Rr = L.results;  bgv = L.betaGenVec(:)';
    C = NaN(numel(Rr), numel(bgv));
    for i = 1:numel(Rr)
        Mc = Rr(i).betaRecCurveMed;
        if isempty(Mc), continue; end
        if size(Mc, 1) ~= numel(PL) || size(Mc, 2) ~= numel(bgv)
            error("szTrial:shape", "%s", sprintf("%s trial %d: betaRecCurveMed is %s, expected %dx%d", ...
                d, i, mat2str(size(Mc)), numel(PL), numel(bgv)));
        end
        C(i, :) = Mc(iF, :);
        [lb, ub, atTop] = deal(NaN(1, numel(PL)));
        for p = 1:numel(PL)
            [lb(p), ub(p), atTop(p)] = bounds_local(Mc(p, :), bgv, LB_FLOOR, LB_RUN, SKIP);
        end
        bs = Rr(i).betaGenStar(:)';
        rows(end+1, :) = {d, Rr(i).f0, bs(iF), string(Rr(i).invertStatus(iF)), lb(iF), ub(iF), atTop(iF), ...
            max(lb), min(ub), all(isfinite(lb))}; %#ok<SAGROW>
    end
    curves.(d) = struct("bgv", bgv, "C", C);
end
T = cell2table(rows, 'VariableNames', ["dataset" "f0" "bgs" "status" "lb" "ub" "atTop" "lbCons" "ubCons" "lbAll"]);
T.inZone     = T.bgs >= T.lb & T.bgs <= T.ub;
T.clear      = T.bgs - T.lb;
T.clearCons  = T.bgs - T.lbCons;
fprintf("Trials with maps: %d\n", height(T));

%% Per dataset summary
S = table();
for d = DATASETS
    x = T(T.dataset == d, :);  inv = x(x.status == "rise" & isfinite(x.bgs), :);
    S = [S; table(d, height(x), 100*mean(isfinite(x.lb)), median(x.lb, "omitnan"), median(x.ub, "omitnan"), ...
        100*mean(x.atTop == 1), median(x.lbCons, "omitnan"), median(x.ubCons, "omitnan"), height(inv), ...
        median(inv.bgs), 100*mean(inv.inZone), median(inv.clear, "omitnan"), 100*mean(inv.clear >= MDC), ...
        median(inv.clearCons, "omitnan"), 100*mean(inv.clearCons >= MDC), ...
        'VariableNames', ["dataset" "nTrials" "pctLB" "LB" "UB" "pctPeakAtTop" "LBcons" "UBcons" ...
        "nInv" "bgsMed" "pctInZone" "clearMed" "pctClearMDC" "clearConsMed" "pctClearConsMDC"])]; %#ok<AGROW>
end
S = outerjoin(S, V002, "Keys", "dataset", "MergeKeys", true, "Type", "left");
S = S(match_local(DATASETS, S.dataset), :);
fprintf("\n%s bounds per trial (median over trials); conservative = across all six pipelines per trial.\n", FOCUS);
fprintf("Clearance = beta_gen* - LB on invertible trials; MDC = %.2f. v002 columns: one grid cell per dataset.\n\n", MDC);
disp(S)
fprintf("Peak at the top of the sweep (no fold within [%.2f, %.2f]): %s\n", curves.(DATASETS(1)).bgv([1 end]), ...
    strjoin(compose("%s %.0f%%", S.dataset, S.pctPeakAtTop), ", "));

%% Figure: per-dataset median SG-IRLS map with IQR, bounds and beta_gen*
fg = figure("Color", "w", "Position", [60 60 1400 640]);  tl = tiledlayout(2, 4, "TileSpacing", "compact");
for k = 1:numel(DATASETS)
    d = DATASETS(k);  c = curves.(d);  s = S(S.dataset == d, :);
    nexttile; hold on;
    q = prctile(c.C, [25 50 75], 1);  ok = all(isfinite(q), 1);
    fill([c.bgv(ok) fliplr(c.bgv(ok))], [q(1, ok) fliplr(q(3, ok))], [0.75 0.85 0.95], "EdgeColor", "none");
    plot(c.bgv, q(2, :), "LineWidth", 1.8, "Color", [0 0.45 0.74]);
    plot([0 max(c.bgv)], [0 max(c.bgv)], "k:");
    xl_local([s.LB s.UB], "--", [0.13 0.55 0.13]);  xl_local(s.bgsMed, "-", [0.85 0.33 0.1]);
    xl_local([s.v002LB s.v002UB], ":", [0.5 0.5 0.5]);
    title(sprintf("%s (n = %d)", strrep(d, "_", " "), s.nTrials), "FontWeight", "normal");
    xlabel("\beta_{gen}");  ylabel("\beta_{rec}");  xlim([0 max(c.bgv)]);  ylim([0 max(c.bgv)]);  box on;
end
nexttile; axis off;
text(0, 0.8, {"Median per-trial " + FOCUS + " map, IQR band", "green dashed: median LB, UB (this script)", ...
    "grey dotted: v002 grid-cell bounds", "orange: median \beta_{gen}* (invertible trials)"}, "FontSize", 9);
title(tl, "Safe zones at each trial's own tempo (v015 per-trial maps)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "T", "S", "LB_FLOOR", "LB_RUN", "SKIP", "MDC", "FOCUS", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function [lb, ub, atTop] = bounds_local(br, bg, floor_, run_, skip_)
% v002 rule on one curve. NaN nodes: bounds are read on the finite part.
    ok = isfinite(br);  lb = NaN;  ub = NaN;  atTop = NaN;
    if sum(ok) < run_ + 2, return; end
    br = br(ok);  bg = bg(ok);
    [~, pk] = max(br);  ub = bg(pk);  atTop = double(pk == numel(br));
    for k = 1 + skip_:pk - run_
        if br(k) > floor_ && all(diff(br(k:k + run_)) > 0), lb = bg(k); break; end
    end
end

function idx = match_local(want, have)
    [ok, idx] = ismember(want, have);
    if ~all(ok), error("szTrial:order", "%s", "Datasets missing from summary"); end
end

function xl_local(x, sty, col)
% xline only where finite (NaN bounds are reported in the table, not drawn)
    for v = x(isfinite(x)), xline(v, sty, "Color", col, "LineWidth", 1.2); end
end
