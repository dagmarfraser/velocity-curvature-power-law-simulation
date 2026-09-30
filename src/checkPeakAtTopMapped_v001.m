%% checkPeakAtTopMapped_v001.m
% safeZoneTrialMaps_v002 (Finding #239) reports the SG-IRLS map peaking at the top of the sweep in
% 87-97% of in-domain trials, over "trials with maps: 3,945". Its test for a missing map is an
% empty betaRecCurveMed (L50), but the 84 trials the runner could not model carry an all-NaN one
% (checkEmptyMaps_v001, checkEmptyMapMechanism_v001), so they stay in the denominator as maps that
% do not peak at the top (bounds_local returns atTop = NaN, L123-125). Its other per-trial figures
% are over invertible trials or omit NaN, so they are unaffected. This recomputes the share:
%   V0 AS PUBLISHED   all rows of safeZoneTrialMaps_v002's T. Anchor: its S.pctPeakAtTop.
%   V1 WITH A MAP     the trials without a map removed.
% Reads:  results/safeZoneTrialMaps_v002.mat, results/checkEmptyMaps_v001.mat
% Writes: results/checkPeakAtTopMapped_v001.mat
% USAGE:  from the project root: checkPeakAtTopMapped_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
SZ_MAT  = fullfile(ROOT, "results", "safeZoneTrialMaps_v002.mat");
EM_MAT  = fullfile(ROOT, "results", "checkEmptyMaps_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "checkPeakAtTopMapped_v001.mat");
IN_DOM  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];   % #242
for f = [SZ_MAT EM_MAT], if ~isfile(f), error("peakTop:input", "%s", "FAILED PATH: " + f); end, end

%% Per dataset
Z = load(SZ_MAT, "T", "S", "LB_FLOORS");  E = load(EM_MAT, "T").T;
P = table();
for d = unique(Z.T.dataset, "stable")'
    z = Z.T(Z.T.dataset == d, :);  e = E(E.dataset == d, :);
    if height(z) ~= height(e) || ~isequaln(z.f0, e.f0)
        error("peakTop:join", "%s", d + ": safe-zone rows do not align with checkEmptyMaps_v001's trials");
    end
    has = e.cause == "mapped";
    p0 = 100 * mean(z.atTop == 1);
    s0 = Z.S.pctPeakAtTop(Z.S.dataset == d & Z.S.floor == Z.LB_FLOORS(1));
    if abs(p0 - s0) > 1e-9, error("peakTop:anchor", "%s", sprintf("%s: V0 %.4f%%, published %.4f%%", d, p0, s0)); end
    P = [P; table(d, ismember(d, IN_DOM), height(z), nnz(~has), p0, 100 * mean(z.atTop(has) == 1), nnz(isnan(z.atTop(has))), ...
        'VariableNames', ["dataset" "inDomain" "nTrials" "nNoMap" "pctPeakAtTop0" "pctPeakAtTop1" "nMappedAtTopNaN"])]; %#ok<AGROW>
end
fprintf("REGRESSION ANCHOR passed: V0 reproduces safeZoneTrialMaps_v002's peak-at-top share in all %d datasets.\n\n", height(P));
disp(P)
d = P(P.inDomain, :);
fprintf("In domain, SG-IRLS map peaks at the sweep top: %.1f-%.1f%% over all trials (as published), %.1f-%.1f%% over the %d trials with a map.\n", ...
    min(d.pctPeakAtTop0), max(d.pctPeakAtTop0), min(d.pctPeakAtTop1), max(d.pctPeakAtTop1), sum(d.nTrials - d.nNoMap));
fprintf("All seven datasets: %d trials with a map (published denominator %d).\n", sum(P.nTrials - P.nNoMap), sum(P.nTrials));

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "P", "runDate");
fprintf("Saved: %s\n", OUT_MAT);
