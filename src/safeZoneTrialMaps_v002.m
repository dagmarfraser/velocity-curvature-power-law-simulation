%% safeZoneTrialMaps_v002.m
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
% v002: v001 (2026-09-28) found the upper bound at the sweep top in 87-97% of in-domain
% trials, but lower bounds far above v002's (Fraser 0.41, Cook CTRL 0.52), with only 41% of
% Fraser's production-invertible trials in zone. Test: is the beta_rec > 0.20 clause (which
% assumes a near-identity map) the cause? LB_FLOORS sweeps it; each row reports agreement of
% "in zone" with the production local check (invertStatus == "rise"). The v001 output stands
% as the record of the rule applied unchanged.
% Read-only; no simulation. Run from src/:  safeZoneTrialMaps_v002
% Inputs: src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Outputs: results/safeZoneTrialMaps_v002.mat, figures/safeZoneTrialMaps_v002.png (floor = LB_FLOORS(end))
% Fraser, D.S. (2026)  v002

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
PL       = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  % runner order (v012 L118)
FOCUS    = "SG-IRLS";
LB_FLOORS = [0.20 0.10 0];  LB_RUN = 3;  SKIP = 1;   % 0.20 = plotForwardMapAtCentroid_v001's clause
MDC      = 0.03;
V002     = table(["Fraser"; "Cook_CTRL"; "Cook_ASD"; "Hickman_PLAC"; "Hickman_HALO"; "Zarandi"], ...
                 [0.27; 0.23; 0.23; 0.23; 0.23; 0.37], [0.50; 0.50; 0.50; 0.50; 0.50; 0.53], ...
                 'VariableNames', ["dataset" "v002LB" "v002UB"]);   % v002 Part 2, grid cells
OUT_MAT  = fullfile(ROOT, "results", "safeZoneTrialMaps_v002.mat");
OUT_PNG  = fullfile(ROOT, "figures", "safeZoneTrialMaps_v002.png");
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
        nF = numel(LB_FLOORS);  [lb, lbC] = deal(NaN(1, nF));  [ub, atTop] = deal(NaN(1, numel(PL)));
        for k = 1:nF
            l = NaN(1, numel(PL));
            for p = 1:numel(PL)
                [l(p), ub(p), atTop(p)] = bounds_local(Mc(p, :), bgv, LB_FLOORS(k), LB_RUN, SKIP);
            end
            lb(k) = l(iF);  lbC(k) = max(l);
        end
        bs = Rr(i).betaGenStar(:)';
        rows(end+1, :) = {d, Rr(i).f0, bs(iF), string(Rr(i).invertStatus(iF)), lb, ub(iF), atTop(iF), ...
            lbC, min(ub)}; %#ok<SAGROW>
    end
    curves.(d) = struct("bgv", bgv, "C", C);
end
T = cell2table(rows, 'VariableNames', ["dataset" "f0" "bgs" "status" "lb" "ub" "atTop" "lbCons" "ubCons"]);
T.inZone     = T.bgs >= T.lb & T.bgs <= T.ub;          % trials x floors
T.clear      = T.bgs - T.lb;
T.clearCons  = T.bgs - T.lbCons;
fprintf("Trials with maps: %d\n", height(T));

%% Per dataset x floor summary; agreement with the production local check
S = table();
for k = 1:numel(LB_FLOORS)
    for d = DATASETS
        x = T(T.dataset == d, :);  inv = x(x.status == "rise" & isfinite(x.bgs), :);
        S = [S; table(LB_FLOORS(k), d, height(x), height(inv), 100*height(inv)/height(x), ...
            median(x.lb(:, k), "omitnan"), median(x.ub, "omitnan"), 100*mean(x.atTop == 1), ...
            median(inv.bgs), 100*mean(inv.inZone(:, k)), median(inv.clear(:, k), "omitnan"), ...
            100*mean(inv.clear(:, k) >= MDC), median(x.lbCons(:, k), "omitnan"), 100*mean(inv.clearCons(:, k) >= MDC), ...
            'VariableNames', ["floor" "dataset" "nTrials" "nInv" "pctInvProd" "LB" "UB" "pctPeakAtTop" ...
            "bgsMed" "pctInvInZone" "clearMed" "pctClearMDC" "LBcons" "pctClearConsMDC"])]; %#ok<AGROW>
    end
end
S = outerjoin(S, V002, "Keys", "dataset", "MergeKeys", true, "Type", "left");
S = sortrows(S, ["floor" "dataset"], ["descend" "ascend"]);
fprintf("\n%s per-trial bounds (median over trials) by LB floor. pctInvProd: production local check (invertStatus rise).\n", FOCUS);
fprintf("pctInvInZone: of those, %% with beta_gen* inside [LB, UB]. Clearance vs MDC %.2f. v002 columns: one grid cell.\n\n", MDC);
disp(S)
fprintf("Peak at the top of the sweep (no fold within [%.2f, %.2f]): %s\n", curves.(DATASETS(1)).bgv([1 end]), ...
    strjoin(compose("%s %.0f%%", S.dataset(S.floor == LB_FLOORS(1)), S.pctPeakAtTop(S.floor == LB_FLOORS(1))), ", "));
SF = S(S.floor == LB_FLOORS(end), :);  SF = SF(match_local(DATASETS, SF.dataset), :);

%% Figure: per-dataset median SG-IRLS map with IQR, bounds and beta_gen*
fg = figure("Color", "w", "Position", [60 60 1400 640]);  tl = tiledlayout(2, 4, "TileSpacing", "compact");
for k = 1:numel(DATASETS)
    d = DATASETS(k);  c = curves.(d);  s = SF(SF.dataset == d, :);
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
title(tl, sprintf("Safe zones at each trial's own tempo (v015 per-trial maps; LB floor %.2f)", LB_FLOORS(end)));
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "T", "S", "LB_FLOORS", "LB_RUN", "SKIP", "MDC", "FOCUS", "-v7.3");
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
