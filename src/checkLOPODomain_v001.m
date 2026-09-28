%% checkLOPODomain_v001.m
% The Abstract's held-out-pipeline figure (83.5% of 21,390 predictions within MDC; Finding
% #205) pools all seven datasets, Zarandi included, while its identifiability counts use the
% six-dataset validity domain (D20, Finding #223). Recompute on the same basis.
% Reads src/loopClosureLOPOOutOfSample_v001.mat (loopClosureLOPOOutOfSample_v001.m):
%   rows = [residOut, residIn, pipelineIdx, datasetIdx] (script L65); criterion |residOut| < MDC (L130).
% Self-check: reproduces the saved all-seven n, median and % < MDC before anything else.
% Also reports the dataset-equal-weighted mean (Fraser is ~75% of cells; Finding #205 caveat).
% Writes results/lopoDomain_v001.mat

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "src", "loopClosureLOPOOutOfSample_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "lopoDomain_v001.mat");
OUT_DOMAIN = "Zarandi";                       % D20

%% Load and self-check
if ~isfile(IN_MAT), error("lopoDomain:input", "%s", "Missing: " + IN_MAT); end
S = load(IN_MAT, "rows", "opts", "perDataset");
R = S.rows;  MDC = S.opts.MDC;  ds = string(S.opts.Datasets);
if size(R, 2) ~= 4, error("lopoDomain:layout", "rows must be [residOut residIn pipelineIdx datasetIdx]"); end
ao = abs(R(:, 1));  dIdx = R(:, 4);
pct = @(x) 100 * mean(x < MDC);
fprintf("All %d datasets: n = %d, median |resid| = %.4f, %% < MDC = %.1f  [Finding #205: 21,390, 0.0081, 83.5]\n", ...
    numel(ds), numel(ao), median(ao), pct(ao));
for k = 1:numel(ds)
    x = ao(dIdx == k);
    if abs(pct(x) - S.perDataset(k).pctMDC) > 1e-9 || numel(x) ~= S.perDataset(k).n
        error("lopoDomain:repro", "%s", "Per-dataset recomputation does not reproduce saved perDataset for " + ds(k));
    end
end
fprintf("Per-dataset values reproduce the saved perDataset struct exactly.\n");

%% Within the validity domain
inDom = ~ismember(dIdx, find(ds == OUT_DOMAIN));
if ~any(ds == OUT_DOMAIN), error("lopoDomain:domain", "%s", OUT_DOMAIN + " not among saved datasets"); end
t = table();
for sel = ["all7" "domain6"]
    m = true(size(ao));  if sel == "domain6", m = inDom; end
    k = unique(dIdx(m))';
    eqw = mean(arrayfun(@(j) pct(ao(dIdx == j)), k));
    t = [t; table(sel, numel(k), sum(m), median(ao(m)), pct(ao(m)), eqw, ...
        'VariableNames', ["basis" "nDatasets" "nPredictions" "medianAbsResid" "pctWithinMDC" "pctDatasetEqualWeighted"])]; %#ok<AGROW>
end
fprintf("\nHeld-out-pipeline validation by basis:\n");  disp(t)

if ~isfolder(fileparts(OUT_MAT)), error("lopoDomain:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "t", "MDC", "OUT_DOMAIN", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
