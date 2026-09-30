%% checkEmptyMaps_v001.m
% Why do some trials have no forward map, and what does that do to the coverage verdicts?
% checkPart2VerdictBasis_v001 found that 15.7% of the failing trial-cells in the in-domain FAIL
% cells have no usable map for that pipeline; analyseTrialTempoWindow_v002 found 84 in-domain
% trials with no SG-IRLS map (Fraser 55, Hickman HALO 27, PLAC 2).
% processOneTrial_local (runLoopClosureFftnoise_v013) skips the forward simulation, leaving an
% all-NaN map, in two places:
%   SHORT     the template-subtracted length M is below 4 x (EDGE_CLIP + 50) samples (L806;
%             EDGE_CLIP defaults to 20, L128: a floor of 280 samples);
%   IRA-FAIL  a per-axis IRASA exponent (alphaMaj or alphaMin) is not finite (L844-851), so the
%             shaped_xu surrogate cannot be built. iraAlphaSigma_v001 returns NaN when fewer than
%             3 bins in 1-20 Hz carry a finite, positive fractal PSD (L60-65); an error thrown
%             inside it is caught and also leaves NaN (L837-843).
% With no finite node the inversion returns "neither" (L630-633) wherever beta_obs is finite, and
% the coverage rule (compare42CellVerdict_v002 L79-84) counts that as a failure. The runner's own
% loop-closure CCC skips these trials (L734).
% This script:
%   1. classifies every trial without a map as SHORT or IRA-FAIL on both replicate bases, and
%      errors if the runner's own no-map test (all-NaN betaRecSlice, L734) disagrees;
%   2. anchors to checkPart2VerdictBasis_v001 (15.7%, v015) and splits that share by cause,
%      including maps that exist but are all NaN for one pipeline;
%   3. describes the trials without a map: length, orbits, tempo, subjects, which axis failed;
%   4. recomputes the coverage verdicts:
%        V0 PRODUCTION  anchor: 4/6/32 (N_REPS = 20, v012) and 8/2/32 (N_REPS = 200, v015);
%        V1             trials the runner skipped are not evaluable (out of the denominator);
%        V2 (v015)      also trial-cells whose map is all NaN for that pipeline (v012 keeps no
%                       curves, so V2 = V1 there).
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v012.mat and _v015.mat
% Writes: results/checkEmptyMaps_v001.mat
% USAGE:  from the project root: checkEmptyMaps_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
IN_DOM    = [true true true true true false false];                    % #242
FS        = [240 133 133 133 133 100 100];                               % FS_DS, runner v013 L185-218
PL        = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
EDGE_CLIP = 20;  M_FLOOR = 4 * (EDGE_CLIP + 50);                         % runner v013 L128, L806
EXPECT    = struct("v012", [4 6 32], "v015", [8 2 32]);
PUB_NOMAP = 15.7;                                                        % checkPart2VerdictBasis_v001, v015
OUT_MAT   = fullfile(ROOT, "results", "checkEmptyMaps_v001.mat");

%% Per trial
T = table();  V = table();  ids = struct();
for k = 1:numel(DATASETS)
    for ver = ["v012" "v015"]
        f = fullfile(ROOT, "src", "loopClosureResults_" + DATASETS(k) + "_all_shaped_xu_" + ver + ".mat");
        if ~isfile(f), error("emptyMaps:input", "%s", "FAILED PATH: " + f); end
        R = load(f, "results").results;  n = numel(R);
        M    = arrayfun(@(r) double(r.M), R(:));
        aMaj = arrayfun(@(r) double(r.alphaMaj), R(:));  aMin = arrayfun(@(r) double(r.alphaMin), R(:));
        cause = repmat("mapped", n, 1);
        cause(M < M_FLOOR) = "short";
        cause(M >= M_FLOOR & ~(isfinite(aMaj) & isfinite(aMin))) = "iraFail";
        skip = cause ~= "mapped";
        seen = arrayfun(@(r) all(isnan(r.betaRecSlice(:))), R(:));        % the runner's own test (L734)
        if ~isequal(seen, skip)
            error("emptyMaps:rule", "%s", sprintf("%s %s: %d trials without a map, %d explained (short %d, IRA-FAIL %d), %d disagree", ...
                DATASETS(k), ver, nnz(seen), nnz(skip), nnz(cause == "short"), nnz(cause == "iraFail"), nnz(seen ~= skip)));
        end
        C = arrayfun(@(r) string(r.invertStatus(:))', R(:), "UniformOutput", false);  st = vertcat(C{:});
        if size(st, 2) ~= 6, error("emptyMaps:width", "%s", "invertStatus is not 6 wide in " + f); end
        id = string({R.trialID})';  ids.(DATASETS(k) + "_" + ver) = id(skip);
        mapNaN = repmat(skip, 1, 6);
        if ver == "v015"
            for i = 1:n
                c = R(i).betaRecCurveMed;
                if isempty(c), mapNaN(i, :) = true; continue; end
                if size(c, 1) ~= 6, error("emptyMaps:curve", "%s", sprintf("%s trial %d: betaRecCurveMed has %d rows", DATASETS(k), i, size(c, 1))); end
                mapNaN(i, :) = all(isnan(c), 2)';
            end
            if ~isequal(all(mapNaN, 2), skip)
                error("emptyMaps:curveRule", "%s", DATASETS(k) + ": all-NaN betaRecCurveMed disagrees with the skip rule");
            end
            bObs = cell2mat(arrayfun(@(r) double(r.betaObs(:))', R(:), "UniformOutput", false));
            f0 = arrayfun(@(r) double(r.f0), R(:));
            T = [T; table(repmat(DATASETS(k), n, 1), repmat(IN_DOM(k), n, 1), (1:n)', id, string({R.subjectID})', ...
                M, M / FS(k), f0, M / FS(k) .* f0, cause, st, mapNaN, isfinite(bObs), aMaj, aMin, ...
                arrayfun(@(r) double(r.alphaIRA), R(:)), arrayfun(@(r) double(r.sigmaMM), R(:)), ...
                'VariableNames', ["dataset" "inDomain" "trial" "trialID" "subjectID" "M" "seconds" "f0" "orbits" ...
                "cause" "status" "mapNaN" "bObsFinite" "alphaMaj" "alphaMin" "alphaIRA" "sigmaMM"])]; %#ok<AGROW>
        end
        for p = 1:6
            valid = st(:, p) ~= "no_beta_obs";
            cv = [mean(st(valid, p) == "rise"), mean(st(valid & ~skip, p) == "rise"), mean(st(valid & ~mapNaN(:, p), p) == "rise")];
            V = [V; table(DATASETS(k), IN_DOM(k), ver, PL(p), nnz(valid), nnz(valid & skip), nnz(valid & mapNaN(:, p)), ...
                cv, arrayfun(@verdict_local, cv), 'VariableNames', ["dataset" "inDomain" "basis" "pipeline" "nValid" ...
                "nSkipValid" "nMapNaNValid" "cov" "verdict"])]; %#ok<AGROW>
        end
    end
end
fprintf("RULE CHECK passed: on both bases a trial has no forward map exactly when M < %d samples or a per-axis IRASA alpha is not finite.\n", M_FLOOR);
same = arrayfun(@(d) isequal(sort(ids.(d + "_v012")), sort(ids.(d + "_v015"))), DATASETS);
fprintf("Same trials skipped on both bases: %s\n", strjoin(DATASETS + "=" + string(same), ", "));

%% Anchors
cnt = @(x) [nnz(x == "PASS") nnz(x == "CONDITIONAL") nnz(x == "FAIL")];
for ver = ["v012" "v015"]
    got = cnt(V.verdict(V.basis == ver, 1));
    if ~isequal(got, EXPECT.(ver)), error("emptyMaps:anchor", "%s", sprintf("%s: V0 gives %s, expected %s", ver, mat2str(got), mat2str(EXPECT.(ver)))); end
end
fprintf("REGRESSION ANCHOR passed: production rule gives 4/6/32 (v012) and 8/2/32 (v015).\n");
nT = height(T);
TC = table(repelem(T.dataset, 6), repelem(T.inDomain, 6), repmat(PL(:), nT, 1), reshape(T.status', [], 1), repelem(T.f0, 6), ...
    repelem(T.cause, 6), reshape(T.mapNaN', [], 1), reshape(T.bObsFinite', [], 1), ...
    'VariableNames', ["dataset" "inDomain" "pipeline" "status" "f0" "cause" "mapNaN" "bObsFinite"]);
V15 = V(V.basis == "v015", :);  F = V15(V15.inDomain & V15.verdict(:, 1) == "FAIL", :);
inF = ismember(TC.dataset + "|" + TC.pipeline, F.dataset + "|" + F.pipeline);
pf = TC(inF & ~ismember(TC.status, ["no_beta_obs" "rise"]) & isfinite(TC.f0), :);   % as checkPart2VerdictBasis_v001 L97-126
noUse = pf.status == "neither" & (pf.mapNaN | ~pf.bObsFinite);
if abs(100 * mean(noUse) - PUB_NOMAP) > 0.05
    error("emptyMaps:anchor2", "%s", sprintf("no-usable-map share %.2f%%, expected %.1f%%", 100 * mean(noUse), PUB_NOMAP));
end
fprintf("ANCHOR passed: %.1f%% of the %d failing trial-cells in the %d in-domain FAIL cells (v015) have no usable map.\n", ...
    100 * mean(noUse), height(pf), height(F));
fprintf("  by cause: SHORT %.1f%%, IRA-FAIL %.1f%%, mapped but all NaN for that pipeline %.1f%%, beta_obs not finite %.1f%% (of all failing trial-cells)\n", ...
    100 * mean(noUse & pf.cause == "short"), 100 * mean(noUse & pf.cause == "iraFail"), ...
    100 * mean(noUse & pf.cause == "mapped" & pf.mapNaN), 100 * mean(noUse & pf.cause == "mapped" & ~pf.mapNaN));

%% Who has no map (v015)
fprintf("\nTrials without a forward map, v015 (SHORT: M < %d samples; IRA-FAIL: per-axis IRASA alpha not finite):\n", M_FLOOR);
S = table();
for k = 1:numel(DATASETS)
    x = T(T.dataset == DATASETS(k), :);  s = x(x.cause ~= "mapped", :);  m = x(x.cause == "mapped", :);
    top = 0;  if height(s) > 0, top = max(groupcounts(s.subjectID)); end
    S = [S; table(DATASETS(k), IN_DOM(k), height(x), height(s), nnz(s.cause == "short"), nnz(s.cause == "iraFail"), ...
        100 * height(s) / height(x), numel(unique(x.subjectID)), numel(unique(s.subjectID)), top, ...
        med_local(s.M), med_local(m.M), med_local(s.seconds), med_local(m.seconds), med_local(s.orbits), med_local(m.orbits), ...
        med_local(s.f0), med_local(m.f0), nnz(any(m.mapNaN, 2)), 100 * mean(reshape(s.status, [], 1) == "neither"), ...
        'VariableNames', ["dataset" "inDomain" "nTrials" "nNoMap" "nShort" "nIraFail" "pctNoMap" "nSubj" "nSubjNoMap" ...
        "maxPerSubj" "MmedNoMap" "MmedMapped" "secNoMap" "secMapped" "orbitsNoMap" "orbitsMapped" "f0NoMap" "f0Mapped" ...
        "nMappedWithNaNPipe" "pctCellsNeither"])]; %#ok<AGROW>
end
S{:, 7:end} = round(S{:, 7:end}, 3);
disp(S)
iF = T(T.cause == "iraFail", :);
if height(iF) > 0
    fprintf("IRA-FAIL (%d trials): alphaMaj alone NaN %d, alphaMin alone NaN %d, both NaN %d; trial-level alphaIRA finite in %d, sigmaMM finite in %d;\n", ...
        height(iF), nnz(~isfinite(iF.alphaMaj) & isfinite(iF.alphaMin)), nnz(isfinite(iF.alphaMaj) & ~isfinite(iF.alphaMin)), ...
        nnz(~isfinite(iF.alphaMaj) & ~isfinite(iF.alphaMin)), nnz(isfinite(iF.alphaIRA)), nnz(isfinite(iF.sigmaMM)));
    fprintf("  M %d-%d samples; %.1f-%.1f orbits (median %.1f); f0 finite in %d\n", min(iF.M), max(iF.M), min(iF.orbits), max(iF.orbits), ...
        median(iF.orbits, "omitnan"), nnz(isfinite(iF.f0)));
    disp(iF(1:min(8, height(iF)), ["dataset" "trialID" "M" "seconds" "f0" "orbits" "alphaMaj" "alphaMin" "alphaIRA" "sigmaMM"]))
end
mp = TC(TC.cause == "mapped" & TC.mapNaN, :);
if height(mp) > 0
    fprintf("Mapped trials with an all-NaN map for one pipeline: %d trial-cells\n", height(mp));
    disp(groupcounts(mp, ["dataset" "pipeline"]))
end

%% Coverage verdicts
for ver = ["v012" "v015"]
    x = V(V.basis == ver, :);
    fprintf("\n%s (%s) PASS/COND/FAIL: V0 %s, V1 %s, V2 %s; in domain (30): %s, %s, %s\n", ver, ...
        ifelse_local(ver == "v012", "N_REPS = 20", "N_REPS = 200"), mat2str(cnt(x.verdict(:, 1))), mat2str(cnt(x.verdict(:, 2))), ...
        mat2str(cnt(x.verdict(:, 3))), mat2str(cnt(x.verdict(x.inDomain, 1))), mat2str(cnt(x.verdict(x.inDomain, 2))), mat2str(cnt(x.verdict(x.inDomain, 3))));
    chg = x(x.verdict(:, 1) ~= x.verdict(:, 2) | x.verdict(:, 1) ~= x.verdict(:, 3), :);
    for r = 1:height(chg)
        fprintf("  %-13s %-9s V0 %5.2f%% %-11s V1 %5.2f%% %-11s V2 %5.2f%% %-11s (%d skipped, %d all-NaN of %d)\n", chg.dataset(r), chg.pipeline(r), ...
            100 * chg.cov(r, 1), chg.verdict(r, 1), 100 * chg.cov(r, 2), chg.verdict(r, 2), 100 * chg.cov(r, 3), chg.verdict(r, 3), ...
            chg.nSkipValid(r), chg.nMapNaNValid(r), chg.nValid(r));
    end
    if isempty(chg), fprintf("  no verdict changes\n"); end
end

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "T", "V", "S", "ids", "M_FLOOR", "runDate");
fprintf("\nSaved: %s\n", OUT_MAT);

%% =========================================================================
function v = verdict_local(c)
    if c >= 0.95, v = "PASS"; elseif c >= 0.90, v = "CONDITIONAL"; else, v = "FAIL"; end
end

function y = med_local(x)
    if isempty(x), y = NaN; else, y = median(x, "omitnan"); end
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
