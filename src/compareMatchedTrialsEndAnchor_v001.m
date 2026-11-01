function out = compareMatchedTrialsEndAnchor_v001(opts)
%COMPAREMATCHEDTRIALSENDANCHOR_V001  Template effect vs population effect in the v015 -> v019 cccMed shift (Addendum 73).
%   v019 has 5-10 more valid trials per dataset (trials the zero-residual tail had invalidated). This script splits the shift:
%     all     cccMed over every valid trial in each version (anchored to constellationMetrics_v006: must reproduce Table E1)
%     matched cccMed(v015) and cccMed(v019) on the trials valid in BOTH, their paired difference (median, 95% bootstrap CI) = the TEMPLATE effect on the same trials
%     new     cccMed(v019) on trials valid in v019 only (and v015-only, if any) = the POPULATION that the tail used to exclude
%   Same validity and metric as constellationMetrics_v006: finite betaObs (6), finite median betaRecSlice (6), finite betaGenStarMed; trials with max|dObs| < 0.005 (degenerate) or all |dObs| < 0.03
%   (scenario A) excluded; per-trial Lin's CCC (linCCC_v001) of the 15 predicted vs 15 observed inter-pipeline contrasts.
%   out = compareMatchedTrialsEndAnchor_v001()    run on the RDS src (five v015 and five v019 files); seconds.
%   Options: Dir, Datasets, VerA ("v015"), VerB ("v019"), NBoot (2000). Writes results/compareMatchedTrialsEndAnchor_v001.mat and .txt.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Fraser"]
    opts.VerA (1,1) string = "v015"
    opts.VerB (1,1) string = "v019"
    opts.NBoot (1,1) double {mustBeInteger, mustBePositive} = 2000
end
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));
dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
if ~isfolder(resDir), mkdir(resDir); end
diary(fullfile(resDir, 'compareMatchedTrialsEndAnchor_v001.txt'));  cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('compareMatchedTrialsEndAnchor_v001 | %s | %s vs %s | data %s\n', string(datetime('now')), opts.VerA, opts.VerB, dataDir);
% anchors: constellationMetrics_v006 cccMed (v015 saved full precision; v019 as printed to 5 d.p. in the 2026-10-29 console)
anc = struct('v015', containers.Map({'Cook_CTRL','Cook_ASD','Hickman_PLAC','Hickman_HALO','Fraser'}, {0.776131107206, 0.68712559954, 0.225184410568, 0.284024016244, 0.888405946777}), ...
             'v019', containers.Map({'Cook_CTRL','Cook_ASD','Hickman_PLAC','Hickman_HALO','Fraser'}, {0.76654, 0.62753, 0.192, 0.22204, 0.88809}));
tolOf = struct('v015', 1e-9, 'v019', 6e-6);
rows = {};  per = struct();
for d = opts.Datasets
    [cA, vA, idA] = trialCCC_local(dataDir, d, opts.VerA);  [cB, vB, idB] = trialCCC_local(dataDir, d, opts.VerB);
    if ~isequal(idA, idB), error('matched:ids', '%s: trial IDs differ between %s and %s', d, opts.VerA, opts.VerB); end
    for v = [opts.VerA opts.VerB]
        if isfield(anc, v) && isKey(anc.(v), char(d))
            m = median(iff_local(v == opts.VerA, cA(vA), cB(vB)), 'omitnan');  b = anc.(v)(char(d));
            if abs(m - b) > tolOf.(v), error('matched:anchor', '%s %s: cccMed %.9f differs from constellationMetrics_v006 (%.9f)', d, v, m, b); end
        end
    end
    both = vA & vB;  onlyB = vB & ~vA;  onlyA = vA & ~vB;
    dd = cB(both) - cA(both);  dd = dd(isfinite(dd));
    bs = bootstrp(opts.NBoot, @(x) median(x), dd);
    rows(end+1, :) = {d, nnz(vA), nnz(vB), nnz(both), nnz(onlyB), nnz(onlyA), median(cA(vA), 'omitnan'), median(cB(vB), 'omitnan'), median(cA(both), 'omitnan'), median(cB(both), 'omitnan'), ...
        median(cB(both), 'omitnan') - median(cA(both), 'omitnan'), median(dd), prctile(bs, 2.5), prctile(bs, 97.5), median(cB(onlyB), 'omitnan'), median(cA(onlyA), 'omitnan')}; %#ok<AGROW>
    per.(char(d)) = struct('ccA', cA, 'ccB', cB, 'vA', vA, 'vB', vB, 'ids', idA);
end
T = cell2table(rows, 'VariableNames', {'dataset','nA','nB','nBoth','nOnlyB','nOnlyA','allA','allB','matchedA','matchedB','dMedMatched','medPairedD','pdLo','pdHi','onlyB_ccc','onlyA_ccc'});
T.dataset = string(T.dataset);
fprintf('anchors OK: per-trial cccMed over all valid trials reproduces constellationMetrics_v006 for every dataset in both versions\n');
fprintf('\nTable E3: v015 -> v019 split into template effect (same trials) and population effect (newly valid trials)\n');  disp(T);
fprintf('allB - allA = total shift; matchedB - matchedA = template effect on the %s trials valid in both; medPairedD = median per-trial paired difference (95%% bootstrap CI pdLo..pdHi).\n', 'same');
fprintf('onlyB_ccc = cccMed of the trials valid ONLY in %s (what the tail used to exclude).\n', opts.VerB);
out = struct('T', T, 'per', per, 'opts', opts, 'runDate', datetime('now'));
save(fullfile(resDir, 'compareMatchedTrialsEndAnchor_v001.mat'), 'out');
fprintf('Saved results/compareMatchedTrialsEndAnchor_v001.mat\n=== COMPLETED OK ===\n');
end

function [cc, valid, ids] = trialCCC_local(dataDir, ds, ver)
% per-trial Lin's CCC and the validity mask, as constellationMetrics_v006 (shaped_xu)
MDC = 0.03;  DEG = 0.005;  nP = 6;  pr = nchoosek(1:nP, 2);  pP = pr(:, 1)';  pQ = pr(:, 2)';
f = fullfile(dataDir, "loopClosureResults_" + ds + "_all_shaped_xu_" + ver + ".mat");
if ~isfile(f), error('matched:missing', 'No file: %s', f); end
S = load(f, 'results', 'config');  r = S.results;  N = numel(r);
if S.config.N_REPS ~= 200 || S.config.N_BETA ~= 25 || S.config.RNG_SEED ~= 1729, error('matched:config', '%s: N_REPS/N_BETA/seed differ from 200/25/1729', f); end
ids = string({r.trialID})';  cc = NaN(N, 1);  valid = false(N, 1);
for t = 1:N
    bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;  g = r(t).betaGenStarMed;
    if numel(bo) ~= nP || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nP || ~isfinite(g), continue; end
    bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
    dO = bo(pP) - bo(pQ);  dP = bp(pP) - bp(pQ);
    if max(abs(dO)) < DEG || all(abs(dO) < MDC), continue; end
    q = linCCC_v001(dP', dO');  cc(t) = q.ccc;  valid(t) = true;
end
end

function v = iff_local(c, a, b), if c, v = a; else, v = b; end, end
