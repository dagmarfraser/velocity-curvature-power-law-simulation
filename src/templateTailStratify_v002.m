function T = templateTailStratify_v002()
%TEMPLATETAILSTRATIFY_V001  Is the template-window effect an edge artefact? Tail (zero-residual stretch) and crop test of SESSION_LOG Addenda 26/28, no HPC.
%   For each trial and sweep setting est[nH,nCW], the template routine's window rule is replayed exactly (templateSubtract_local, loopClosureCore_v004 lines 478-513): win0 = round(nCW/f0*fs), hop = round(win0/2),
%   win = max(min(win0, N), round(2/f0*fs)); windows s = 1, 1+hop, ... while s+win-1 <= N; crop cl = max(round(win/4), 1); analysed segment cs = cl to ce = min(max(N - cl, cs + 100), N); samples after the last window
%   end keep residual 0, so zero stretch z = max(0, ce - lastEnd) samples (all of the segment if no window fits). Reported in cycles (z*f0/fs) with the analysed length L (cycles).
%   Per dataset and setting: dz = z - z(baseline est[4,4]), dL = L - L(baseline) and the per-trial hit change D = hit(setting) - hit(est[4,4]) (same trial, same run). Output: (1) D by dz group (<= 0.05 cycle, > 0.05) with
%   subject-clustered 95% CIs; (2) pooled over settings, OLS D = b0 + b1*dz + b2*dL with subject-clustered bootstrap CIs and R^2 (the tail and crop mechanisms entered together; b1 < 0 = more zero stretch costs hits);
%   (3) the same fit excluding whole-trace-clamped settings. Descriptive and post hoc; f0 is the geometry-pass f0 (as the stratifiers), trial selection and beta_gen rebuilt and checked as templateCyclesStratify_v002.
%   v002 also saves the baseline's tail share (percent of trials with zero stretch above DZ_TOL), median stretch and analysed length per dataset (v001 printed them only). Run from src/. Saves results/templateTailStratify_v002.mat.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));  addpath(genpath(fullfile(srcDir, 'req')));
resDir = fullfile(srcDir, '..', 'results');  SEED = 3141;  N = 432;  MDC = 0.03;  NBOOT = 1000;  DZ_TOL = 0.05;
rdsRes = '/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results';
runs = {'templateVariantTest_v002.mat', 'templateVariantTest_v003_sweep.mat', 'templateVariantTest_v003_sweep2.mat', 'templateVariantTest_v004_sweep3.mat'};
BAD2 = ["est[8,2]" "est[8,2.5]" "est[8,4]" "est[8,5]"];   % sweep2 columns holding other settings (Addendum 25); sweep3 re-runs them
RR = cell(1, numel(runs));
for r = 1:numel(runs)
    f = fullfile(resDir, runs{r});  if ~isfile(f), f = fullfile(rdsRes, runs{r}); end
    if ~isfile(f), error('tail:missing', '%s not found locally or on the RDS (mounted?)', runs{r}); end
    RR{r} = load(f);
end
T = struct('dataset', {}, 'tab', {}, 'fit', {}, 'fitNoClamp', {}, 'baseTailShare', {}, 'baseMedZ', {}, 'baseMedL', {});
for d = ["Hickman PLAC" "Zarandi" "Fraser"]
    [matFile, sigToMM, FS, importFn] = dataset_local(d);
    bio = load(fullfile(srcDir, matFile), 'bioResults').bioResults;  trials = importFn();  tids = string({trials.trialID})';
    [tf, loc] = ismember(string(bio.trialID), tids);
    rows = find(tf & isfinite(bio.ira_alphaMean) & isfinite(bio.sigmaMean) & bio.sigmaMean > 0 & bio.ira_alphaMean <= 6 & bio.ira_alphaMean >= -2);
    rng(SEED, 'twister');  nU = numel(rows);  sel = zeros(N, 1);
    for q = 1:ceil(N / nU), p = rows(randperm(nU));  k = min(nU, N - (q - 1) * nU);  sel((q - 1) * nU + (1:k)) = p(1:k); end
    bgV = 0.05 + 0.60 * rand(N, 1);  subj = string(bio.subjectID(sel));
    for r = 1:numel(RR), Rr = RR{r}.R(strcmp(string({RR{r}.R.dataset}), d));  if max(abs(Rr.betaGen - bgV)) > 1e-12, error('tail:rebuild', 'Rebuilt beta_gen differs from %s for %s', runs{r}, d); end, end
    uq = unique(sel);  f0u = nan(size(uq));  nu = nan(size(uq));
    for i = 1:numel(uq), tr = trials(loc(uq(i)));  g = loopClosureCore_v004(struct('x', double(tr.x(:)), 'y', double(tr.y(:))), FS, sigToMM, GeometryOnly=true);  f0u(i) = g.f0;  nu(i) = numel(tr.x); end
    [~, ix] = ismember(sel, uq);  f0 = f0u(ix);  Nn = nu(ix);  [zB, LB] = edge_local(Nn, f0, FS, 4);
    rowsOut = {};  stack = [];   % stack: [setting index, subject idx, D, dz, dL, clamped]
    [~, ~, sg] = unique(subj);  seen = strings(0);
    for r = 1:numel(RR)
        Rr = RR{r}.R(strcmp(string({RR{r}.R.dataset}), d));  nA = numel(Rr.labels);  M = Rr.med(:, 1:nA);  ok = all(isfinite(M), 2);  H = abs(M - bgV) <= MDC;
        for a = 2:nA
            lab = string(Rr.labels(a));
            if ismember(lab, ["nor" "orc"]) || (strcmp(runs{r}, 'templateVariantTest_v003_sweep2.mat') && ismember(lab, BAD2)) || ismember(lab, seen), continue, end
            tok = regexp(lab, '^est\[(\d+),([\d.]+)\]$', 'tokens', 'once');  if isempty(tok), continue, end
            seen(end + 1) = lab; %#ok<AGROW>
            nCW = str2double(tok{2});  [z, L, clamped] = edge_local(Nn, f0, FS, nCW);  dz = z - zB;  dL = L - LB;  D = double(H(:, a)) - double(H(:, 1));
            lo = ok & dz <= DZ_TOL;  hi = ok & dz > DZ_TOL;
            rowsOut(end + 1, :) = {char(lab), median(z(ok)), 100 * mean(z(ok) > DZ_TOL), nnz(hi), 100 * mean(D(hi)), ci_local(D, hi, sg, NBOOT, SEED), 100 * mean(D(lo)), ci_local(D, lo, sg, NBOOT, SEED + 1), median(dL(ok))}; %#ok<AGROW>
            stack = [stack; repmat(numel(seen), nnz(ok), 1), sg(ok), D(ok), dz(ok), dL(ok), clamped(ok)]; %#ok<AGROW>
        end
    end
    fprintf('\n%s: N %d trials; baseline est[4,4] median zero stretch %.2f cycles (%.0f%% of trials > %.2f), analysed %.2f cycles\n', d, N, median(zB), 100 * mean(zB > DZ_TOL), DZ_TOL, median(LB));
    fprintf('  %-11s | %7s %6s | %5s %7s %17s | %7s %17s | %6s\n', 'setting', 'med z', '%z>tol', 'n dz+', 'D dz+', '[95% CI]', 'D dz<=0', '[95% CI]', 'dL cyc');
    for i = 1:size(rowsOut, 1), fprintf('  %-11s | %7.2f %6.0f | %5d %+7.1f [%+6.1f,%+6.1f] | %+7.1f [%+6.1f,%+6.1f] | %+6.2f\n', rowsOut{i, 1:4}, rowsOut{i, 5}, rowsOut{i, 6}, rowsOut{i, 7}, rowsOut{i, 8}, rowsOut{i, 9}); end
    fit = fit_local(stack, NBOOT, SEED);  fit2 = fit_local(stack(stack(:, 6) == 0, :), NBOOT, SEED);
    fprintf('  pooled OLS D(pp) = b0 + b1*dz + b2*dL: b1 (per cycle of extra zero stretch) %+.1f [%+.1f,%+.1f], b2 (per cycle of extra analysed length) %+.1f [%+.1f,%+.1f], R2 %.3f, rows %d\n', 100 * fit.b(2), 100 * fit.ci(:, 2)', 100 * fit.b(3), 100 * fit.ci(:, 3)', fit.r2, fit.n);
    fprintf('  excluding whole-trace-clamped settings:                  b1 %+.1f [%+.1f,%+.1f], b2 %+.1f [%+.1f,%+.1f], R2 %.3f, rows %d\n', 100 * fit2.b(2), 100 * fit2.ci(:, 2)', 100 * fit2.b(3), 100 * fit2.ci(:, 3)', fit2.r2, fit2.n);
    T(end + 1) = struct('dataset', d, 'tab', {rowsOut}, 'fit', fit, 'fitNoClamp', fit2, 'baseTailShare', 100 * mean(zB > DZ_TOL), 'baseMedZ', median(zB), 'baseMedL', median(LB)); %#ok<AGROW>
end
save(fullfile(resDir, 'templateTailStratify_v002.mat'), 'T');  fprintf('\nSaved results/templateTailStratify_v002.mat\n=== TAIL DONE ===\n');
end

function [z, L, clamped] = edge_local(Nn, f0, fs, nCW)
% Replays templateSubtract_local's window rule per trial: z = zero-residual samples inside the analysed segment (cycles), L = analysed length (cycles), clamped = window >= trace.
    n = numel(Nn);  z = nan(n, 1);  L = z;  clamped = false(n, 1);
    for i = 1:n
        N = Nn(i);  win0 = round(nCW / f0(i) * fs);  hop = round(win0 / 2);  win = max(min(win0, N), round(2 / f0(i) * fs));  clamped(i) = win0 >= N;
        cl = max(round(win / 4), 1);  cs = cl;  ce = min(max(N - cl, cs + 100), N);
        if win > N, lastEnd = 0; else, lastEnd = (win - 1) + hop * floor((N - win) / hop) + 1; end   % end of the last window: s = 1 + k*hop, e = s + win - 1
        z(i) = max(0, ce - max(lastEnd, cs - 1)) * f0(i) / fs;  L(i) = (ce - cs + 1) * f0(i) / fs;
    end
end

function ci = ci_local(D, m, sg, nb, seed)
% Subject-clustered bootstrap 95% CI [lo hi] (pp) of the mean of D over mask m; NaN if fewer than 5 trials.
    idx = find(m);  if numel(idx) < 5, ci = [NaN NaN];  return, end
    us = unique(sg(idx));  ii = arrayfun(@(u) idx(sg(idx) == u), us, 'UniformOutput', false);  rng(seed, 'twister');  bs = zeros(nb, 1);
    for k = 1:nb, pk = vertcat(ii{randi(numel(us), numel(us), 1)});  bs(k) = 100 * mean(D(pk)); end
    ci = prctile(bs, [2.5 97.5])';  ci = ci(:)';
end

function fit = fit_local(S, nb, seed)
% S columns: setting, subject, D, dz, dL, clamped. Subject-clustered bootstrap of the OLS coefficients.
    X = [ones(size(S, 1), 1) S(:, 4) S(:, 5)];  b = X \ S(:, 3);  e = S(:, 3) - X * b;  r2 = 1 - sum(e .^ 2) / sum((S(:, 3) - mean(S(:, 3))) .^ 2);
    us = unique(S(:, 2));  ii = arrayfun(@(u) find(S(:, 2) == u), us, 'UniformOutput', false);  rng(seed, 'twister');  B = zeros(nb, 3);
    for k = 1:nb, pk = vertcat(ii{randi(numel(us), numel(us), 1)});  B(k, :) = (X(pk, :) \ S(pk, 3))'; end
    fit = struct('b', b, 'ci', prctile(B, [2.5 97.5], 1), 'r2', r2, 'n', size(S, 1));
end

function [matFile, sigToMM, FS, importFn] = dataset_local(name)
    switch name
        case "Fraser",       matFile = 'noiseCharacterisation_fraser.mat';       sigToMM = 1/10.41793; FS = 240; importFn = @() importDB_fraser_v001(Verbose=false);
        case "Hickman PLAC", matFile = 'noiseCharacterisation_hickmanPLAC.mat';  sigToMM = 0.248;      FS = 133; importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Zarandi",      matFile = 'noiseCharacterisation_zarandi.mat';      sigToMM = 10.0;       FS = 100; importFn = @() importDB_zarandi_v001(Verbose=false);
        otherwise, error('tail:dataset', 'Dataset "%s" not supported', name);
    end
end
