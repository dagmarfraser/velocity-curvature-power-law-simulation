function Ssum = analyseCheckAdditiveVsRetimed_v001(trialFile, opts)
% analyseCheckAdditiveVsRetimed_v001  Subject-level analysis of a saved
% checkAdditiveVsRetimed trial table. Re-runnable in seconds; never
% re-runs the corpus. Subject is the unit; trials are never pooled.
%
%   1. lambda, lambdaLow, lambdaHigh: per-subject median, Wilcoxon
%      signed-rank vs 0 (additive, the forward map's model) and vs 1
%      (re-timed).
%   2. Band contrast: per-subject median lambdaLow vs lambdaHigh, paired
%      signed-rank, subjects with both finite.
%   3. Tempo: Spearman correlation of per-subject median lambda and f0.
% Tables from checkAdditiveVsRetimed_v001 (no bands, no f0) are accepted:
% the missing parts are skipped with a warning, never filled.
%
% Output: results/checkAdditiveVsRetimed_subjectSummary_<trialFile stem>.mat
% (refuses to overwrite unless Overwrite=true).
%
% USAGE:
%   S = analyseCheckAdditiveVsRetimed_v001( ...
%       "results/checkAdditiveVsRetimed_v002_20260922_205130.mat");
%
% Fraser, D.S. (2026) v001 (analysis moved out of
% runCheckAdditiveVsRetimedFullCorpus_v002, logic unchanged)

arguments
    trialFile (1,1) string {mustBeFile}
    opts.MinSubj   (1,1) double {mustBeInteger, mustBePositive} = 6
    opts.Save      (1,1) logical = true
    opts.Overwrite (1,1) logical = false
end

L = load(trialFile);
if ~isfield(L, "T")
    error("analyseCheckAdditiveVsRetimed:NoTable", "%s", trialFile + " holds no trial table T.");
end
T = L.T;
runOpts = struct(); if isfield(L, "opts"), runOpts = L.opts; end
hasBands = all(ismember(["lambdaLow","lambdaHigh"], T.Properties.VariableNames));
hasF0    = ismember("f0", T.Properties.VariableNames);
if ~hasBands, warning("analyseCheckAdditiveVsRetimed:NoBands", "%s", "No band columns (v001 table): band analyses skipped."); end
if ~hasF0,    warning("analyseCheckAdditiveVsRetimed:NoF0",    "%s", "No f0 column (v001 table): tempo analysis skipped."); end
if isfield(runOpts, "TrialSelection") && runOpts.TrialSelection ~= "all"
    warning("analyseCheckAdditiveVsRetimed:NotFullCorpus", "%s", sprintf( ...
        "Trial table is TrialSelection=""%s"", NSurr=%d: not the full-corpus run.", ...
        runOpts.TrialSelection, runOpts.NSurr));
end

lamVars = "lambda"; if hasBands, lamVars = ["lambda","lambdaLow","lambdaHigh"]; end
ds = unique(T.dataset, "stable")';
rows = {};
for d = ds
    m = T.dataset == d;
    r = table(d, numel(unique(T.subjectID(m))), sum(m), 'VariableNames', ["dataset","nSubj","nTrials"]);
    for v = ["lambda","lambdaLow","lambdaHigh"]
        for pre = ["nSubj_","med_","pVs0_","pVs1_"], r.(pre + v) = NaN; end
        if ~ismember(v, lamVars), continue, end
        lamS = subjMedian_local(T.(v)(m), T.subjectID(m));
        r.("nSubj_" + v) = numel(lamS);
        r.("med_" + v)   = medOrNaN_local(lamS);
        [p0, p1] = srTests_local(lamS, opts.MinSubj, d + " " + v);
        r.("pVs0_" + v) = p0; r.("pVs1_" + v) = p1;
    end
    r.nSubj_band = NaN; r.medLowMinusHigh = NaN; r.pLowVsHigh = NaN;
    if hasBands
        both = m & isfinite(T.lambdaLow) & isfinite(T.lambdaHigh);
        [lo, g] = subjMedian_local(T.lambdaLow(both), T.subjectID(both));
        hi = [];
        if ~isempty(g), hi = splitapply(@median, T.lambdaHigh(both), g); end
        r.nSubj_band = numel(lo); r.medLowMinusHigh = medOrNaN_local(lo - hi);
        if numel(lo) >= opts.MinSubj, r.pLowVsHigh = signrank(lo, hi); end
    end
    r.rhoLamF0 = NaN; r.pLamF0 = NaN; r.medF0_eval = NaN; r.medF0_all = NaN;
    if hasF0
        ok = m & isfinite(T.lambda) & isfinite(T.f0);
        [lamS, g] = subjMedian_local(T.lambda(ok), T.subjectID(ok));
        if numel(lamS) >= opts.MinSubj
            f0S = splitapply(@median, T.f0(ok), g);
            [rs, ps] = corr(lamS, f0S, "Type", "Spearman");
            r.rhoLamF0 = rs; r.pLamF0 = ps;
        end
        r.medF0_eval = medOrNaN_local(T.f0(ok));
        r.medF0_all  = medOrNaN_local(T.f0(m & isfinite(T.f0)));
    end
    rows{end+1} = r; %#ok<AGROW>

    oth = T.status(m & ~ismember(T.status, ["ok","anchorsUnseparated"]));
    if ~isempty(oth)
        [u, ~, gi] = unique(oth);
        fprintf("  %s exclusions: %s\n", d, strjoin(u + "=" + string(accumarray(gi, 1)), ", "));
    end
end
Ssum = vertcat(rows{:});

fprintf("\n%s\n(1) Wilcoxon signed-rank on per-subject median lambda (total / cycle band / high band)\n", trialFile);
disp(Ssum(:, ["dataset","nSubj","nSubj_lambda","med_lambda","pVs0_lambda","pVs1_lambda", ...
    "nSubj_lambdaLow","med_lambdaLow","pVs0_lambdaLow","pVs1_lambdaLow", ...
    "nSubj_lambdaHigh","med_lambdaHigh","pVs0_lambdaHigh","pVs1_lambdaHigh"]))
fprintf("(2) Band contrast (paired signed-rank)  (3) Tempo (Spearman), subject medians\n");
disp(Ssum(:, ["dataset","nSubj_band","medLowMinusHigh","pLowVsHigh","rhoLamF0","pLamF0","medF0_eval","medF0_all"]))

if opts.Save
    [fDir, stem] = fileparts(trialFile);
    outFile = fullfile(fDir, "checkAdditiveVsRetimed_subjectSummary_" + ...
        erase(stem, "checkAdditiveVsRetimed_") + ".mat");
    if isfile(outFile) && ~opts.Overwrite
        error("analyseCheckAdditiveVsRetimed:Exists", "%s", ...
            outFile + " exists; pass Overwrite=true to replace it.");
    end
    MinSubj = opts.MinSubj; %#ok<NASGU>
    save(outFile, "Ssum", "trialFile", "runOpts", "MinSubj");
    fprintf("Saved: %s\n", outFile);
end
end

% -------------------------------------------------------------------------
function [med, g] = subjMedian_local(vals, subj)
% Per-subject median over finite values; g indexes vals(isfinite(vals)).
    ok = isfinite(vals); med = NaN(0,1); g = [];
    if ~any(ok), return, end
    g = findgroups(subj(ok));
    med = splitapply(@median, vals(ok), g);
end

function v = medOrNaN_local(x)
    if isempty(x), v = NaN; else, v = median(x); end
end

function [p0, p1] = srTests_local(x, minN, label)
    p0 = NaN; p1 = NaN;
    if numel(x) < minN
        warning("analyseCheckAdditiveVsRetimed:FewSubjects", "%s", ...
            sprintf("%s: %d evaluable subjects (< %d), tests not run.", label, numel(x), minN));
        return
    end
    p0 = signrank(x, 0); p1 = signrank(x, 1);
end
