%% checkPart2VerdictBasis_v001.m
% Two VERIFY flags in v004 that rest on the per-trial verdicts:
%   (1) Part 2 Method through-line (v004 L1595): "32 of 42 dataset-pipeline combinations
%       (76.2%) still fail the per-trial check, and those failures are dominated by slow
%       trials whose maps are flat (#227, #239)". Checks the count on both N_REPS bases and
%       on the validity domain, then the composition of the failing trials in each in-domain
%       FAIL cell: tempo (below 0.35 Hz, #227's slowest bin; below the dataset's median f0),
%       status (below the map / above it / on a fold branch) and map gain.
%   (2) Reviewer concern (1)(a) (v004 L3100): "99.3-100% of the production-invertible trials
%       sit on their rising branch" (#239). Reads #239's saved summary and reports the floor-0
%       in-zone share per dataset, with the other two #239 quantities (peak at the top of the
%       sweep; clearance of beta_gen* above the rising branch's start).
% Verdict rule identical to compare42CellVerdict_v002 L73-92 (rise coverage over trials with a
% valid invertStatus; PASS >= 95%, CONDITIONAL >= 90%, FAIL below).
% SELF-CHECKS (Fail Loud): verdict counts reproduce 4/6/32 (N_REPS = 20, v012) and 8/2/32
%   (N_REPS = 200, v015) over 42 cells, and 0/6/24 and 4/2/24 on the 30 in-domain cells
%   ("Checked and consistent", coherence pass v003_v001); #239's peak-at-top shares reproduce
%   87-97% in-domain.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v012.mat and _v015.mat,
%         results/safeZoneTrialMaps_v002.mat
% Writes: results/checkPart2VerdictBasis_v001.mat
% USAGE:  from the project root: checkPart2VerdictBasis_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
IN_DOM   = [true true true true true false false];                    % #242
PL       = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];   % runner order
SLOW_HZ  = 0.35;                                                      % #227's slowest bin
EXPECT   = struct("N20", [4 6 32], "N200", [8 2 32], "N20dom", [0 6 24], "N200dom", [4 2 24]);
SZ_MAT   = fullfile(ROOT, "results", "safeZoneTrialMaps_v002.mat");
OUT_MAT  = fullfile(ROOT, "results", "checkPart2VerdictBasis_v001.mat");
if ~isfile(SZ_MAT), error("p2basis:input", "%s", "FAILED PATH: " + SZ_MAT); end

%% (1a) Verdicts on both bases
V = table();  Tr = table();
for k = 1:numel(DATASETS)
    d = DATASETS(k);
    for ver = ["v012" "v015"]
        f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_" + ver + ".mat");
        if ~isfile(f), error("p2basis:corpus", "%s", "FAILED PATH: " + f); end
        L = load(f, "results", "betaGenVec");  R = L.results;  bgv = L.betaGenVec(:)';
        if any(arrayfun(@(r) numel(r.invertStatus), R) ~= 6)
            error("p2basis:width", "%s", "invertStatus is not 6 wide in " + f);
        end
        st = string(cell2mat_local(arrayfun(@(r) string(r.invertStatus(:))', R, "UniformOutput", false)));
        for p = 1:6
            valid = st(:, p) ~= "no_beta_obs";  cov = mean(st(valid, p) == "rise");
            V = [V; table(d, IN_DOM(k), ver, PL(p), nnz(valid), cov, verdict_local(cov), ...
                'VariableNames', ["dataset" "inDomain" "basis" "pipeline" "nValid" "coverage" "verdict"])]; %#ok<AGROW>
        end
        if ver == "v015"                                                % per-trial rows for (1b)
            n = numel(R);  f0 = arrayfun(@(r) double(r.f0), R(:));  g = NaN(n, 6);  pos = repmat("nomap", n, 6);
            for i = 1:n
                Mc = R(i).betaRecCurveMed;
                if isempty(Mc), continue; end
                g(i, :) = ((Mc(:, end) - Mc(:, 1)) / (bgv(end) - bgv(1)))';   % chord slope of each map (NaN ends stay NaN)
                for q = 1:6                                               % beta_obs against the map's range (as trialTempoWindow)
                    c = Mc(q, :);  b = R(i).betaObs(q);
                    if all(isnan(c)) || ~isfinite(b), pos(i, q) = "nan";
                    elseif b < min(c), pos(i, q) = "below";  elseif b > max(c), pos(i, q) = "above";  else, pos(i, q) = "within"; end
                end
            end
            Tr = [Tr; table(repelem(d, 6 * n, 1), repelem(IN_DOM(k), 6 * n, 1), repelem((1:n)', 6), repelem(f0, 6), ...
                repmat(PL(:), n, 1), reshape(st', [], 1), reshape(g', [], 1), reshape(pos', [], 1), ...
                'VariableNames', ["dataset" "inDomain" "trial" "f0" "pipeline" "status" "chordSlope" "obsPos"])]; %#ok<AGROW>
        end
    end
end
cnt = @(x) [nnz(x == "PASS") nnz(x == "CONDITIONAL") nnz(x == "FAIL")];
got = struct("N20", cnt(V.verdict(V.basis == "v012")), "N200", cnt(V.verdict(V.basis == "v015")), ...
    "N20dom", cnt(V.verdict(V.basis == "v012" & V.inDomain)), "N200dom", cnt(V.verdict(V.basis == "v015" & V.inDomain)));
for fn = string(fieldnames(EXPECT))'
    if ~isequal(got.(fn), EXPECT.(fn))
        error("p2basis:anchor", "%s", sprintf("%s: PASS/COND/FAIL %s, expected %s", fn, mat2str(got.(fn)), mat2str(EXPECT.(fn))));
    end
end
fprintf("SELF-CHECK passed: 42 cells %s (N20) and %s (N200); 30 in-domain cells %s and %s.\n", ...
    mat2str(got.N20), mat2str(got.N200), mat2str(got.N20dom), mat2str(got.N200dom));
fprintf("FAIL cells: 32 of 42 (76.2%%) on both bases; in the validity domain 24 of 30 (80.0%%) on both.\n");
chg = V(V.basis == "v012", :);  c2 = V(V.basis == "v015", :);
fprintf("Cells whose verdict differs between bases: %d of 42\n", nnz(chg.verdict ~= c2.verdict));
for r = find(chg.verdict ~= c2.verdict)'
    fprintf("  %s %s: %.2f%% %s -> %.2f%% %s\n", chg.dataset(r), chg.pipeline(r), 100 * chg.coverage(r), chg.verdict(r), ...
        100 * c2.coverage(r), c2.verdict(r));
end
fprintf("Fraser, all six cells (N_REPS = 20 -> 200):\n");
for r = find(chg.dataset == "Fraser")'
    fprintf("  %-10s %.2f%% %-11s -> %.2f%% %s\n", chg.pipeline(r), 100 * chg.coverage(r), chg.verdict(r), 100 * c2.coverage(r), c2.verdict(r));
end

%% (1b) Composition of failing trials in each in-domain FAIL cell (v015)
F = c2(c2.inDomain & c2.verdict == "FAIL", :);
comp = table();
for r = 1:height(F)
    x = Tr(Tr.dataset == F.dataset(r) & Tr.pipeline == F.pipeline(r) & Tr.status ~= "no_beta_obs" & isfinite(Tr.f0), :);
    medF = median(Tr.f0(Tr.dataset == F.dataset(r) & Tr.pipeline == F.pipeline(r) & isfinite(Tr.f0)));
    bad = x(x.status ~= "rise", :);  good = x(x.status == "rise", :);
    comp = [comp; table(F.dataset(r), F.pipeline(r), height(x), height(bad), 100 * F.coverage(r), ...
        100 * mean(bad.f0 < SLOW_HZ), 100 * mean(x.f0 < SLOW_HZ), 100 * mean(bad.f0 < medF), ...
        median(bad.f0), median(good.f0), 100 * mean(bad.status == "neither"), ...
        100 * mean(ismember(bad.status, ["ambiguous" "desc"])), ...
        median(bad.chordSlope, "omitnan"), median(good.chordSlope, "omitnan"), ...
        100 * mean(bad.status == "neither" & ismember(bad.obsPos, ["below" "above"])), ...
        100 * mean(bad.status == "neither" & bad.obsPos == "within"), ...
        'VariableNames', ["dataset" "pipeline" "nValid" "nFail" "coveragePct" "pctFailSlow" "pctAllSlow" ...
        "pctFailBelowMedF0" "f0MedFail" "f0MedRise" "pctFailNeither" "pctFailFold" "slopeFail" "slopeRise" ...
        "pctFailOffMap" "pctFailNeitherWithin"])]; %#ok<AGROW>
end
comp{:, 5:end} = round(comp{:, 5:end}, 3);
fprintf("\n(1b) In-domain FAIL cells (v015): who fails. pctFailSlow = %% of failing trials below %.2f Hz\n", SLOW_HZ);
fprintf("     (pctAllSlow = the cell's base rate); slope = median chord slope of the trial's map over the sweep.\n");
disp(comp)
pooled = Tr(Tr.inDomain & Tr.status ~= "no_beta_obs" & isfinite(Tr.f0), :);
inFail = ismember(pooled.dataset + "|" + pooled.pipeline, F.dataset + "|" + F.pipeline);
pf = pooled(inFail & pooled.status ~= "rise", :);  pa = pooled(inFail, :);
fprintf("Pooled over the %d in-domain FAIL cells: %d failing trial-cells; %.1f%% below %.2f Hz (base rate %.1f%%); ", ...
    height(F), height(pf), 100 * mean(pf.f0 < SLOW_HZ), SLOW_HZ, 100 * mean(pa.f0 < SLOW_HZ));
medT = groupsummary(Tr(isfinite(Tr.f0), :), "dataset", "median", "f0");
[~, im] = ismember(pf.dataset, medT.dataset);
fprintf("%.1f%% below their dataset's median f0; %.1f%% 'neither', %.1f%% on a fold branch; chord slope %.2f (failing) vs %.2f (rising).\n", ...
    100 * mean(pf.f0 < medT.median_f0(im)), 100 * mean(pf.status == "neither"), 100 * mean(ismember(pf.status, ["ambiguous" "desc"])), ...
    median(pf.chordSlope, "omitnan"), median(pa.chordSlope(pa.status == "rise"), "omitnan"));
nb = pf.status == "neither";
fprintf("Of the failing trial-cells: 'neither' %.1f%% = off the map (beta_obs below or above its range) %.1f%% + no usable map for that pipeline (all nodes NaN or no map) %.1f%% + within its range on no monotonic run %.1f%%; on a fold branch %.1f%%.\n", ...
    100 * mean(nb), 100 * mean(nb & ismember(pf.obsPos, ["below" "above"])), 100 * mean(nb & ismember(pf.obsPos, ["nan" "nomap"])), ...
    100 * mean(nb & pf.obsPos == "within"), 100 * mean(ismember(pf.status, ["ambiguous" "desc"])));
fprintf("  (off the map: below %.1f%%, above %.1f%% of failing trial-cells)\n", 100 * mean(nb & pf.obsPos == "below"), 100 * mean(nb & pf.obsPos == "above"));
byDs = groupsummary(pf, "dataset", @(x) 100 * mean(x < SLOW_HZ), "f0");
disp(byDs)

%% (2) #239's per-trial rising-branch statement
S = load(SZ_MAT, "S").S;
S0 = S(S.floor == 0, :);
[ok, loc] = ismember(DATASETS, S0.dataset);
if ~all(ok), error("p2basis:sz", "%s", "safeZoneTrialMaps_v002 summary lacks a dataset"); end
S0 = S0(loc, :);  S0.inDomain = IN_DOM(:);
pk = S0.pctPeakAtTop(S0.inDomain);
if abs(round(min(pk)) - 87) > 0 || abs(round(max(pk)) - 97) > 0
    error("p2basis:szAnchor", "%s", sprintf("peak-at-top in-domain %.1f-%.1f, #239 says 87-97", min(pk), max(pk)));
end
fprintf("\n(2) #239, SG-IRLS, floor 0 (pure monotonicity), production-invertible trials:\n");
disp(S0(:, ["dataset" "inDomain" "nTrials" "nInv" "pctInvProd" "pctPeakAtTop" "pctInvInZone" "clearMed" "pctClearMDC"]))
fprintf("In-domain: in zone %.1f-%.1f%%; peak at the sweep top %.0f-%.0f%%; clearance median %.2f-%.2f, >= MDC in %.0f-%.0f%%.\n", ...
    min(S0.pctInvInZone(S0.inDomain)), max(S0.pctInvInZone(S0.inDomain)), min(pk), max(pk), ...
    min(S0.clearMed(S0.inDomain)), max(S0.clearMed(S0.inDomain)), min(S0.pctClearMDC(S0.inDomain)), max(S0.pctClearMDC(S0.inDomain)));

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "V", "comp", "S0", "SLOW_HZ", "runDate");
fprintf("\nSaved: %s\n", OUT_MAT);

%% =========================================================================
function v = verdict_local(c)
    if c >= 0.95, v = "PASS"; elseif c >= 0.90, v = "CONDITIONAL"; else, v = "FAIL"; end
end

function M = cell2mat_local(C)
    M = vertcat(C{:});
end
