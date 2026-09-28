%% tabulateBetaGenStarRunDrift_v001.m
% Finding #204 extended to all seven datasets: run-to-run drift in
% betaGenStar between the confirmatory corpus (v013, or v014 where v013 is
% absent -- Cook ASD, Dhieb) and the v015 re-run. Same RNG_SEED, N_REPS,
% config; the only difference is run identity (parfor worker streams).
%
% Definitions match Finding #204's own Zarandi table: |diff| in betaGenStar
% over cells finite in BOTH runs, per pipeline. Unlike runLoopClosureAll_v015
% verify_local (L186-193), ciLo/ciHi are NOT pooled in. Trials are aligned by
% trialID, asserted, not assumed from index order.
%
% Adds dsMedShift: median(betaGenStar v015) - median(v013) over jointly
% finite trials, i.e. how far a dataset-level median moves. NOTE this is a
% plain per-pipeline median, NOT the constellation-median estimator the
% pooled beta_gen* is built from (see claude.md "Estimator trap").
%
% Provenance: written 2026-09-18, Session 112.

%% CONFIG
SRC_DIR  = fullfile(fileparts(mfilename("fullpath")));
DATASETS = ["Zarandi","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Dhieb","Fraser"];
MODEL    = "shaped_xu";
MDC      = 0.03;                                 % clinical MDC, prereg s3.2
OUT_MAT  = fullfile(SRC_DIR, "betaGenStarRunDrift_v001.mat");
OUT_TXT  = fullfile(SRC_DIR, "betaGenStarRunDrift_v001.txt");

%% RUN
rows = {};
for ds = DATASETS
    [r13, r15, labels, refVer] = loadPair_local(SRC_DIR, ds, MODEL);
    [B13, B15, S13, S15] = stack_local(r13, r15, numel(labels), ds);
    for p = 0:numel(labels)                      % p = 0 -> all pipelines
        if p == 0, cols = 1:numel(labels); pl = "ALL";
        else,      cols = p;               pl = labels(p); end
        rows(end+1,:) = [{ds, refVer, pl}, ...
            num2cell(driftStats_local(B13(:,cols), B15(:,cols), S13(:,cols), S15(:,cols), MDC))]; %#ok<SAGROW>
    end
end
T = cell2table(rows, "VariableNames", ["dataset","ref","pipeline","nTrials", ...
    "nJoint","median","p90","max","fracAboveMDC","nanFlip","statusFlip","dsMedShift"]);

%% REPORT
disp(T);
save(OUT_MAT, "T", "MDC", "DATASETS", "MODEL");
writetable(T, OUT_TXT, "Delimiter", "\t");
fprintf("Saved %s\n      %s\n", OUT_MAT, OUT_TXT);

%% ---------------------------------------------------------------------------
function [r13, r15, labels, refVer] = loadPair_local(srcDir, ds, model)
f15 = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_v015.mat", ds, model));
refVer = "v013";
f13 = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_v013.mat", ds, model));
if ~isfile(f13)
    refVer = "v014";
    f13 = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_v014.mat", ds, model));
end
if ~isfile(f15), error("drift:missing", "%s", "No v015 corpus: " + f15); end
if ~isfile(f13), error("drift:missing", "%s", "No v013/v014 comparator for " + ds); end
A = load(f13, "results", "pipelineLabels");
B = load(f15, "results", "pipelineLabels");
if ~isequal(string(A.pipelineLabels), string(B.pipelineLabels))
    error("drift:labels", "%s", ds + ": pipelineLabels differ between runs.");
end
labels = string(B.pipelineLabels);
r13 = A.results; r15 = B.results;
end

function [B13, B15, S13, S15] = stack_local(r13, r15, nP, ds)
id13 = string({r13.trialID}); id15 = string({r15.trialID});
if numel(unique(id13)) ~= numel(id13)
    error("drift:dupID", "%s", ds + ": duplicate trialIDs in comparator.");
end
[tf, loc] = ismember(id13, id15);
if ~all(tf) || numel(id13) ~= numel(id15)
    error("drift:align", "%s", sprintf("%s: trial sets differ (%d vs %d, %d unmatched).", ...
        ds, numel(id13), numel(id15), nnz(~tf)));
end
r15 = r15(loc);
n = numel(r13);
B13 = NaN(n, nP); B15 = B13; S13 = strings(n, nP); S15 = S13;
for t = 1:n
    B13(t,:) = r13(t).betaGenStar(:).'; B15(t,:) = r15(t).betaGenStar(:).';
    s13 = string(r13(t).invertStatus); s15 = string(r15(t).invertStatus);
    if numel(s13) ~= nP || numel(s15) ~= nP
        error("drift:status", "%s", sprintf("%s trial %d: invertStatus has %d/%d " + ...
            "elements, expected %d -- field is not per-pipeline.", ds, t, numel(s13), numel(s15), nP));
    end
    S13(t,:) = s13(:).'; S15(t,:) = s15(:).';
end
end

function out = driftStats_local(a, b, sa, sb, mdc)
joint = isfinite(a) & isfinite(b);
d = abs(a(joint) - b(joint));
if isempty(d), q = [NaN NaN NaN NaN];
else,          q = [median(d), quantile(d, 0.90), max(d), mean(d > mdc)]; end
shift = NaN;
if size(a,2) == 1 && any(joint), shift = median(b(joint)) - median(a(joint)); end
out = [size(a,1), nnz(joint), q, nnz(isfinite(a) ~= isfinite(b)), nnz(sa ~= sb), shift];
end
