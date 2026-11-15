function out = loopClosureLOPOBaselines_v002(opts)
% loopClosureLOPOBaselines_v002  v001 on the N=48 Fraser corpus (2026-11-04). ONLY CHANGES: the self-check against Finding #243 (20,514 cells,
% 84.53% within MDC, Fraser 16,104 cells) cannot hold when Fraser is N=48 (16,496 cells), so it is split: the four unchanged
% datasets (Cook CTRL 417, Cook ASD 426, Hickman PLAC 1,903, Hickman HALO 1,664 cells) must reproduce Finding #243's counts EXACTLY,
% and the total / percentage anchor is replaced by printing the N=45 values beside the N=48 ones (Fraser's count and the pooled
% rates are reported, never forced). Output results/lopoBaselines_v002.mat. v001's header follows.
%LOOPCLOSURELOPOBASELINES_V001  D-I1: does the out-of-sample leave-one-pipeline-out figure beat simple baselines?
%
%   The framework's prediction (loopClosureLOPOOutOfSample_v001): for a trial and a held-out pipeline p, take the median
%   of the other five pipelines' beta_gen*, push it through p's OWN per-trial forward map, and ask whether the result lands
%   within the MDC (0.03) of what p observed. This script scores that prediction against four baselines on the SAME cells
%   (every cell where the framework predicts, the Finding #243 domain: Cook CTRL, Cook ASD, Hickman PLAC, Hickman HALO, Fraser):
%     B1  naive mean of the other five pipelines' beta_obs.
%     B2  a linear fit of p's beta_obs on that mean, per dataset and pipeline, cross-validated in folds of whole subjects
%         (so it is fitted on other participants and scored on unseen ones). A generous baseline: it learns each dataset's
%         own offset and slope.
%     B3  mismatched map, same dataset: the framework's own input (median of the other five beta_gen*) pushed through the
%         forward map of ANOTHER participant's trial from the same dataset. Tests whether the per-trial map matters.
%     B4  wrong-colour map: the same, with a donor trial from a different dataset whose noise exponent (alphaIRA) differs by at
%         least ColourGap. Tests whether the map's noise model matters.
%   B3 and B4 draw NDonor donors per cell (seeded); each draw is scored.
%
%   RATES. Primary: within-MDC rate among cells the predictor defines (a donor map that does not cover the input gives NaN,
%   counted and reported, never clamped). Secondary: NaN counted as a miss. Pooled over cells, and with the five datasets
%   weighted equally (Fraser is 78% of cells). MARGIN: framework minus baseline in percentage points with a 95% interval from a
%   paired cluster bootstrap over participants (dataset + subjectID), pooled and stratified by dataset.
%   "Beats" means the interval's lower bound is above zero on both bases.
%
%   SELF-CHECK (fail loud): the framework scored here must reproduce Finding #243's five-dataset figures: 20,514 cells, 84.5%
%   within MDC (84.53), 77.7% equal-weighted, and the per-dataset counts. Any difference is an error, not a result.
%
%   RUNS ON THE HPC (reads the v015 corpora; minutes, no pool). Run from src/. Writes results/lopoBaselines_v001.mat.
%   Fraser, D.S. (2026)  v001

arguments
    opts.Datasets   (1,:) string = ["Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Fraser"]
    opts.NoiseModel (1,1) string = "shaped_xu"
    opts.Version    (1,1) string = "v015"
    opts.MinOthers  (1,1) double = 3
    opts.MDC        (1,1) double = 0.03
    opts.NDonor     (1,1) double = 5
    opts.ColourGap  (1,1) double = 0.5
    opts.NFold      (1,1) double = 10
    opts.NBoot      (1,1) double = 2000
    opts.Seed       (1,1) double = 1729
    opts.Anchor     (1,1) logical = true
    opts.Save       (1,1) logical = true
end
DEFAULT_DS = ["Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Fraser"];
if opts.Anchor && ~isequal(opts.Datasets, DEFAULT_DS)
    error('lopoBase:anchorDomain', '%s', 'The Finding #243 self-check applies to the default five datasets; pass Anchor=false for others.');
end
srcDir = fileparts(mfilename('fullpath'));
nD = numel(opts.Datasets);  R = opts.NDonor;  MDC = opts.MDC;
PRED = ["Framework (own map)" "B1 naive mean of other five obs" "B2 linear fit, CV by participant" ...
        "B3 mismatched map, same dataset" "B4 wrong-colour map"];

%% Load: flat arrays over all trials
bg = [];  allB = [];  allG = [];  allC = [];  allA = [];  allDS = [];  allKey = strings(0, 1);
for k = 1:nD
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_%s_%s.mat', opts.Datasets(k), opts.NoiseModel, opts.Version));
    if ~isfile(f), error('lopoBase:missingCorpus', 'Not found:\n  %s', f); end
    S = load(f, 'results', 'betaGenVec');  r = S.results;
    if ~isfield(r, 'betaRecCurveMed'), error('lopoBase:noCurve', '%s has no betaRecCurveMed (pre-v015 corpus).', f); end
    if isempty(bg), bg = S.betaGenVec(:)'; elseif ~isequal(bg, S.betaGenVec(:)'), error('lopoBase:grid', '%s: betaGenVec differs from the first dataset', f); end
    n = numel(r);  NB = numel(bg);
    B = zeros(n, 6);  G = zeros(n, 6);  C = zeros(n, 6, NB);
    for i = 1:n
        B(i, :) = r(i).betaObs(:)';  G(i, :) = r(i).betaGenStar(:)';
        if ~isequal(size(r(i).betaRecCurveMed), [6 NB]), error('lopoBase:curveShape', '%s trial %d: curve %s', opts.Datasets(k), i, mat2str(size(r(i).betaRecCurveMed))); end
        C(i, :, :) = reshape(r(i).betaRecCurveMed, [1 6 NB]);
    end
    allB = [allB; B]; allG = [allG; G]; allC = [allC; C]; allA = [allA; [r.alphaIRA]']; %#ok<AGROW>
    allDS = [allDS; k * ones(n, 1)]; allKey = [allKey; opts.Datasets(k) + "|" + string({r.subjectID})']; %#ok<AGROW>
    fprintf('%-14s %5d trials\n', opts.Datasets(k), n);
end
[ukey, ~, gid] = unique(allKey);  nG = numel(ukey);  N = numel(gid);
validMap = squeeze(sum(isfinite(allC), 3)) >= 2;                    % N x 6

%% Cells: trial x held-out pipeline where the framework predicts
cds = [];  ci = [];  cp = [];  cx = [];  cy = [];  cm = [];  cF = [];  nAttempt = 0;  nOOR = 0;
for i = 1:N
    for p = 1:6
        if ~isfinite(allB(i, p)), continue, end
        oth = setdiff(1:6, p);  og = allG(i, oth);  og = og(isfinite(og));
        if numel(og) < opts.MinOthers, continue, end
        nAttempt = nAttempt + 1;
        x = median(og);
        pf = interp_local(bg, squeeze(allC(i, p, :))', x);
        if isnan(pf), nOOR = nOOR + 1; continue, end
        ob = allB(i, oth);
        cds(end + 1, 1) = allDS(i); ci(end + 1, 1) = i; cp(end + 1, 1) = p; cx(end + 1, 1) = x; %#ok<AGROW>
        cy(end + 1, 1) = allB(i, p); cm(end + 1, 1) = mean(ob(isfinite(ob))); cF(end + 1, 1) = pf; %#ok<AGROW>
    end
end
nC = numel(cy);  cg = gid(ci);
fprintf('\nCells: %d framework predictions (%d attempts, %d out of range)\n', nC, nAttempt, nOOR);

%% Self-check against Finding #243 (v002: the four unchanged datasets exact; Fraser and pooled reported against N=45)
if opts.Anchor
    want = struct('Cook_CTRL', 417, 'Cook_ASD', 426, 'Hickman_PLAC', 1903, 'Hickman_HALO', 1664, 'Fraser', 16104);   % Finding #243 (N=45)
    for k = 1:nD
        nk = nnz(cds == k);  nm = char(opts.Datasets(k));
        if nm == "Fraser"
            fprintf('Fraser: %d cells (Finding #243, N=45: %d) -- reported, not an anchor in v002.\n', nk, want.Fraser);
        elseif nk ~= want.(nm)
            error('lopoBase:anchorN', '%s: %d cells, Finding #243 has %d', nm, nk, want.(nm));
        end
    end
    pct = 100 * mean(abs(cF - cy) < MDC);
    fprintf('Framework, N=48: %d cells, %.2f%% within MDC (Finding #243, N=45: 20,514 cells, 84.53%%). Unchanged datasets reproduce Finding #243 exactly.\n', nC, pct);
end

%% Baselines
rng(opts.Seed, 'twister');
pB1 = cm;                                                           % B1
pB2 = nan(nC, 1);                                                   % B2: CV by participant, per dataset and pipeline
for k = 1:nD
    for p = 1:6
        m = find(cds == k & cp == p);
        if numel(m) < 5, continue, end
        ug = unique(cg(m));  K = min(opts.NFold, numel(ug));  perm = ug(randperm(numel(ug)));
        fold = zeros(numel(m), 1);
        for q = 1:numel(perm), fold(cg(m) == perm(q)) = mod(q - 1, K) + 1; end
        for f = 1:K
            te = fold == f;  tr = ~te;
            if nnz(tr) < 3 || numel(unique(cm(m(tr)))) < 2, continue, end
            c = polyfit(cm(m(tr)), cy(m(tr)), 1);
            pB2(m(te)) = polyval(c, cm(m(te)));
        end
    end
end
pB3 = nan(nC, R);  pB4 = nan(nC, R);  noDonor3 = 0;  noDonor4 = 0;   % B3, B4
for c = 1:nC
    i = ci(c);  p = cp(c);
    cand3 = find(allDS == cds(c) & validMap(:, p) & gid ~= gid(i));
    cand4 = find(allDS ~= cds(c) & validMap(:, p) & abs(allA - allA(i)) >= opts.ColourGap);
    if isempty(cand3), noDonor3 = noDonor3 + 1; else
        for q = 1:R, j = cand3(randi(numel(cand3))); pB3(c, q) = interp_local(bg, squeeze(allC(j, p, :))', cx(c)); end
    end
    if isempty(cand4), noDonor4 = noDonor4 + 1; else
        for q = 1:R, j = cand4(randi(numel(cand4))); pB4(c, q) = interp_local(bg, squeeze(allC(j, p, :))', cx(c)); end
    end
end
if noDonor3 > 0 || noDonor4 > 0, warning('lopoBase:noDonor', '%d cells had no same-dataset donor, %d no wrong-colour donor (left undefined, counted).', noDonor3, noDonor4); end

%% Scores
P = {cF, pB1, pB2, pB3, pB4};  nP = numel(P);
hit = cell(1, nP);  def = cell(1, nP);  rep = [1 1 1 R R];
for s = 1:nP
    def{s} = isfinite(P{s});  hit{s} = def{s} & abs(P{s} - cy) < MDC;
end
cluH = zeros(nG, nP);  cluD = zeros(nG, nP);  cluA = zeros(nG, nP);
for s = 1:nP
    cluH(:, s) = accumarray(cg, sum(hit{s}, 2), [nG 1]);  cluD(:, s) = accumarray(cg, sum(def{s}, 2), [nG 1]);
    cluA(:, s) = accumarray(cg, rep(s) * ones(nC, 1), [nG 1]);
end
gds = accumarray(gid, allDS, [nG 1], @max);                         % dataset of each participant
tab = struct('name', {}, 'definedPct', {}, 'ratePooled', {}, 'rateNaNmiss', {}, 'rateEqual', {}, 'medAbsErr', {}, 'perDataset', {});
for s = 1:nP
    pdr = arrayfun(@(k) sum(cluH(gds == k, s)) / max(sum(cluD(gds == k, s)), 1), 1:nD);
    e = abs(P{s}(:) - repmat(cy, size(P{s}, 2), 1));
    tab(s) = struct('name', PRED(s), 'definedPct', 100 * sum(cluD(:, s)) / sum(cluA(:, s)), 'ratePooled', 100 * sum(cluH(:, s)) / sum(cluD(:, s)), ...
        'rateNaNmiss', 100 * sum(cluH(:, s)) / sum(cluA(:, s)), 'rateEqual', 100 * mean(pdr), 'medAbsErr', median(e, 'omitnan'), 'perDataset', 100 * pdr);
end

%% Paired cluster bootstrap of the margin (framework minus baseline, percentage points)
rng(opts.Seed + 1, 'twister');
dP = nan(opts.NBoot, nP - 1);  dE = nan(opts.NBoot, nP - 1);
byDS = arrayfun(@(k) find(gds == k), 1:nD, 'UniformOutput', false);
for b = 1:opts.NBoot
    ix = randi(nG, nG, 1);
    rate = @(ii, s) sum(cluH(ii, s)) / sum(cluD(ii, s));
    rk = zeros(nD, nP);
    for k = 1:nD
        gk = byDS{k};  ik = gk(randi(numel(gk), numel(gk), 1));
        for s = 1:nP, rk(k, s) = rate(ik, s); end
    end
    for s = 2:nP
        dP(b, s - 1) = 100 * (rate(ix, 1) - rate(ix, s));
        dE(b, s - 1) = 100 * (mean(rk(:, 1)) - mean(rk(:, s)));
    end
end
ciP = prctile(dP, [2.5 97.5]);  ciE = prctile(dE, [2.5 97.5]);

%% Report
fprintf('\n%-36s %8s %9s %9s %9s %8s | %-26s %-26s\n', 'predictor', 'defined%', 'rate%', 'NaN=miss', 'equal-wt', 'medAbsE', 'margin vs framework, pooled', 'margin, datasets equal');
for s = 1:nP
    t = tab(s);
    if s == 1
        fprintf('%-36s %8.1f %9.2f %9.2f %9.2f %8.4f |\n', t.name, t.definedPct, t.ratePooled, t.rateNaNmiss, t.rateEqual, t.medAbsErr);
    else
        fprintf('%-36s %8.1f %9.2f %9.2f %9.2f %8.4f | %+6.1f pp [%+6.1f, %+6.1f]   %+6.1f pp [%+6.1f, %+6.1f]\n', t.name, t.definedPct, t.ratePooled, ...
            t.rateNaNmiss, t.rateEqual, t.medAbsErr, tab(1).ratePooled - t.ratePooled, ciP(1, s - 1), ciP(2, s - 1), ...
            tab(1).rateEqual - t.rateEqual, ciE(1, s - 1), ciE(2, s - 1));
    end
end
fprintf('\nPer dataset, %% within MDC among defined cells (columns: %s)\n', strjoin(opts.Datasets, ', '));
for s = 1:nP, fprintf('%-36s %s\n', tab(s).name, sprintf('%8.1f', tab(s).perDataset)); end
beats = ciP(1, :) > 0 & ciE(1, :) > 0;
for s = 2:nP, fprintf('Framework beats %-34s : %s\n', PRED(s), string(beats(s - 1))); end

out = struct('predictors', PRED, 'table', tab, 'marginPooled', tab(1).ratePooled - [tab(2:end).ratePooled], 'ciPooled', ciP, ...
    'marginEqual', tab(1).rateEqual - [tab(2:end).rateEqual], 'ciEqual', ciE, 'beats', beats, 'nCells', nC, 'nAttempt', nAttempt, 'nOOR', nOOR, ...
    'datasets', opts.Datasets, 'opts', opts, 'generated', string(datetime('now')));
if opts.Save
    resDir = fullfile(srcDir, '..', 'results');  if ~isfolder(resDir), mkdir(resDir); end
    o = fullfile(resDir, 'lopoBaselines_v002.mat');  save(o, '-struct', 'out');
    fprintf('\nSaved: %s\n', o);
end
end

function y = interp_local(bg, c, x)
% Linear interpolation on a forward map; NaN outside its support (never clamped).
    m = isfinite(c) & isfinite(bg);
    if nnz(m) < 2 || ~isfinite(x), y = NaN; return, end
    bgm = bg(m);  cm = c(m);
    if x < min(bgm) || x > max(bgm), y = NaN; return, end
    [u, iu] = unique(bgm);
    y = interp1(u, cm(iu), x, 'linear');
end
