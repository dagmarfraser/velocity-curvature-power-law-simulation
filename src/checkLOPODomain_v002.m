%% checkLOPODomain_v002.m
% Held-out-pipeline (LOPO) validation on the five-dataset validity domain (Finding #242).
% v001 (Finding #231) excluded Zarandi only (D20), so its "validity domain" basis
% (20,592 predictions, 84.5% within MDC, 76.7% dataset-equal-weighted) still held Dhieb.
% Since Session 125, Dhieb is outside the domain on measurement contamination (ground B:
% pointer acceleration on at collection). v002 adds that basis; nothing else changes.
% Reads src/loopClosureLOPOOutOfSample_v001.mat (loopClosureLOPOOutOfSample_v001.m):
%   rows = [residOut, residIn, pipelineIdx, datasetIdx]; criterion |residOut| < MDC.
% Self-checks (errors, not warnings): per-dataset values reproduce the saved perDataset
%   struct; the all-seven and six-dataset bases reproduce Findings #205 and #231.
% Writes results/lopoDomain_v002.mat. Run from src/.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "src", "loopClosureLOPOOutOfSample_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "lopoDomain_v002.mat");
BASES   = struct("name", {"all7", "domain6_F231", "domain5_F242"}, ...
                 "out",  {string.empty, "Zarandi", ["Zarandi" "Dhieb"]});
PUBLISHED = [21390 83.5 74.1; 20592 84.5 76.7];      % n, pooled %, equal-weighted % (#205, #231)
TOL_PCT   = 0.05;                                      % published values are rounded to 0.1

%% Load and self-check per dataset
if ~isfile(IN_MAT), error("lopoDomain:input", "%s", "Missing: " + IN_MAT); end
S = load(IN_MAT, "rows", "opts", "perDataset");
R = S.rows;  MDC = S.opts.MDC;  ds = string(S.opts.Datasets);
if size(R, 2) ~= 4, error("lopoDomain:layout", "rows must be [residOut residIn pipelineIdx datasetIdx]"); end
ao = abs(R(:, 1));  dIdx = R(:, 4);
pct = @(x) 100 * mean(x < MDC);
for k = 1:numel(ds)
    x = ao(dIdx == k);
    if abs(pct(x) - S.perDataset(k).pctMDC) > 1e-9 || numel(x) ~= S.perDataset(k).n
        error("lopoDomain:repro", "%s", "Per-dataset recomputation does not reproduce saved perDataset for " + ds(k));
    end
end
fprintf("Per-dataset values reproduce the saved perDataset struct exactly.\n");

%% Bases
t = table();
for b = BASES
    miss = setdiff(b.out, ds);
    if ~isempty(miss), error("lopoDomain:domain", "%s", "Not among saved datasets: " + strjoin(miss, ", ")); end
    m = ~ismember(dIdx, find(ismember(ds, b.out)));
    k = unique(dIdx(m))';
    eqw = mean(arrayfun(@(j) pct(ao(dIdx == j)), k));
    t = [t; table(string(b.name), numel(k), sum(m), median(ao(m)), pct(ao(m)), eqw, ...
        'VariableNames', ["basis" "nDatasets" "nPredictions" "medianAbsResid" "pctWithinMDC" "pctDatasetEqualWeighted"])]; %#ok<AGROW>
end

%% Self-check against the published bases
for i = 1:2
    ok = t.nPredictions(i) == PUBLISHED(i, 1) && abs(t.pctWithinMDC(i) - PUBLISHED(i, 2)) <= TOL_PCT ...
        && abs(t.pctDatasetEqualWeighted(i) - PUBLISHED(i, 3)) <= TOL_PCT;
    if ~ok, error("lopoDomain:published", "%s", "Basis " + t.basis(i) + " does not reproduce its published values"); end
end
fprintf("Bases all7 and domain6_F231 reproduce Findings #205 and #231.\n");

fprintf("\nHeld-out-pipeline validation by basis:\n");  disp(t)
fprintf("Per dataset (n, %% within MDC):\n");
disp(table(ds(:), [S.perDataset.n]', round([S.perDataset.pctMDC]', 1), 'VariableNames', ["dataset" "n" "pctWithinMDC"]))

if ~isfolder(fileparts(OUT_MAT)), error("lopoDomain:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "t", "MDC", "BASES", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
