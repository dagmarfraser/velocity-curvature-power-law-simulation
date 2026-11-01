function R = templateEndAnchorTest_v001(opts)
%TEMPLATEENDANCHORTEST_V001  Causal test of the template zero-residual tail (SESSION_LOG Addenda 39, 47). EXPLORATORY, post hoc, not a prereg analysis.
%   The template routine (loopClosureCore_v004 templateSubtract_local) advances windows by half a window and requires each to fit, so samples after the last window keep
%   residual exactly 0 inside the analysed segment ("tail"). The END-ANCHORED variant adds one final window ending at the last sample (s = N - win + 1) when the
%   regular windows stop short; windows, harmonics, weights and crop (cl = round(win/4)) are otherwise unchanged. Nothing else differs from the core (text-patched copies).
%   PART 1 (simulation, identical trials/seeds/beta_gen as templateVariantTest_v004: seed 3141, N = 432, runI2_v005 base world, shaped_xu): arms [nH nCW anchored] =
%     [4 4 0] production, [4 4 1], and [8 2], [8 5], [8 8] each without and with anchoring. Score: pooled-median beta_gen* within MDC of beta_gen. Reported: hit %, paired
%     contrasts with subject-clustered bootstrap 95% CIs: (a) anchoring effect per setting, (b) window effect vs [4,4] without anchoring, (c) window effect vs [4,4] with
%     anchoring. If the tail drives the window effects, (c) is near 0 where (b) is not, and (a) is positive where the setting adds tail.
%   PART 2 (real trials; bounds the headline): per core dataset, NREAL random trials (seeded); the core run at the production setting [4,4], unanchored and anchored, same
%     seed per trial, NReps = 20, NR7 = 0 (paired; the v015 headline used N_REPS = 200). Tail trials flagged by the production window rule (zero stretch > 0.05 cycle).
%     Reported by tail / non-tail: median paired change in beta_gen* (anchored - production) with subject-clustered bootstrap CI, change in invertible share, and the
%     dataset median beta_gen*. Non-tail trials should be unchanged (internal control).
%   Chunks (resumable) in results/tea1_<part>_<dataset>_chunks. Output results/templateEndAnchorTest_v001.mat. Anchors: replayed analysed length equals the core's saved M.
%   HPC: open in the MATLAB GUI (sinteractive, 72 cores) and Run, no arguments (see runTemplateEndAnchorHPC_v001). Expect about 6 h; pauses safely at the deadline; re-run resumes.
%   Fraser, D.S. (2026)  v001
arguments
    opts.Datasets (1,:) string = ["Hickman PLAC" "Zarandi" "Fraser"]       % part 1, as templateVariantTest_v004
    opts.RealDatasets (1,:) string = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"]
    opts.NTrials  (1,1) double = 432
    opts.NReal    (1,1) double = 300
    opts.Chunk    (1,1) double = 144
    opts.Seed     (1,1) double = 3141
    opts.Variants (:,3) double = [4 4 0; 4 4 1; 8 2 0; 8 2 1; 8 5 0; 8 5 1; 8 8 0; 8 8 1]
    opts.MDC      (1,1) double = 0.03
    opts.DzTol    (1,1) double = 0.05
    opts.NBoot    (1,1) double = 2000
    opts.NWorkers (1,1) double = 72
    opts.MaxFailFrac (1,1) double = 0.05
    opts.Deadline (1,1) datetime = NaT
    opts.Tag      (1,1) string = ""
end
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));  addpath(genpath(fullfile(srcDir, 'req')));
resDir = fullfile(srcDir, '..', 'results');  if ~isfolder(resDir), mkdir(resDir); end
outFile = fullfile(resDir, 'templateEndAnchorTest_v001' + opts.Tag + '.mat');
if isfile(outFile), error('tea:exists', 'Output exists, not overwriting: %s', outFile); end
if ~isequal(opts.Variants(1, :), [4 4 0]), error('tea:baseline', 'First variant must be the production [4 4 0]'); end
for f = ["synthWorld_v001.m" "loopClosureCore_v004.m"], if isempty(which(f)), error('tea:missing', 'Missing %s on the path', f); end, end
V = opts.Variants;  nV = size(V, 1);  [fns, tplDir] = makeVariants_local(srcDir, V);
if opts.NWorkers > 0 && isempty(gcp('nocreate')), parpool('local', opts.NWorkers, 'IdleTimeout', Inf); end
if opts.NWorkers > 0, pctRunOnAll(sprintf("addpath('%s')", tplDir)); end
sfx = ["" "+EA"];  labs = compose("[%g,%g]", V(:, 1), V(:, 2)) + sfx(V(:, 3) + 1)';
fprintf('Arms: %s\n', strjoin(labs, '  '));
R1 = struct('dataset', {}, 'n', {}, 'hit', {}, 'labels', {}, 'contrast', {}, 'betaGen', {}, 'med', {}, 'subj', {});
R2 = struct('dataset', {}, 'n', {}, 'nTail', {}, 'tab', {}, 'bgProd', {}, 'bgEA', {}, 'isTail', {}, 'subj', {}, 'invProd', {}, 'invEA', {});

%% PART 1 -- simulation
for d = opts.Datasets
    [matFile, sigToMM, FS, importFn] = datasetSim_local(d, srcDir);
    bio = load(matFile, 'bioResults').bioResults;  trials = importFn();  tids = string({trials.trialID})';
    [tf, loc] = ismember(string(bio.trialID), tids);
    rows = find(tf & isfinite(bio.ira_alphaMean) & isfinite(bio.sigmaMean) & bio.sigmaMean > 0 & bio.ira_alphaMean <= 6 & bio.ira_alphaMean >= -2);
    rng(opts.Seed, 'twister');  N = opts.NTrials;  nU = numel(rows);
    sel = zeros(N, 1);  for q = 1:ceil(N / nU), p = rows(randperm(nU));  k = min(nU, N - (q - 1) * nU);  sel((q - 1) * nU + (1:k)) = p(1:k); end
    bgV = 0.05 + 0.60 * rand(N, 1);  alC = double(bio.ira_alphaMean(sel));  sgC = double(bio.sigmaMean(sel));  subjSel = string(bio.subjectID(sel));
    trC = arrayfun(@(i) struct('x', double(trials(loc(i)).x(:)), 'y', double(trials(loc(i)).y(:))), sel);  clear trials
    fprintf('\nPART 1 %s: %d trials from %d unique geometries | FS %d Hz\n', d, N, nU, FS);
    sg = struct('seed', opts.Seed, 'n', N, 'chunk', opts.Chunk, 'var', V, 'bg', bgV, 'ds', d);
    ckDir = fullfile(resDir, "tea1_sim_" + strrep(d, " ", "_") + opts.Tag + "_chunks");  if ~isfolder(ckDir), mkdir(ckDir); end
    nCh = ceil(N / opts.Chunk);  chunkSecs = [];  Med = nan(N, nV);  msg = strings(N, 1);  seed0 = opts.Seed;
    for c = 1:nCh
        idx = (c - 1) * opts.Chunk + 1 : min(c * opts.Chunk, N);  ck = fullfile(ckDir, sprintf('chunk_%03d.mat', c));
        if isfile(ck)
            K = load(ck);  if ~isequal(K.sg, sg) || ~isequal(K.idx, idx), error('tea:chunkMismatch', 'Checkpoint %s was written with different options; move %s away.', ck, ckDir); end
            Med(idx, :) = K.m;  msg(idx) = K.ms;  chunkSecs(end + 1) = K.secs; %#ok<AGROW>
            fprintf('  chunk %d/%d loaded (took %.1f min)\n', c, nCh, K.secs / 60);  continue
        end
        if ~isnat(opts.Deadline) && datetime('now') + seconds(1.15 * max([chunkSecs 0])) > opts.Deadline
            fprintf('\ntea PAUSED before %s chunk %d/%d (not an error). Run again to continue.\n', d, c, nCh);  R = [];  return
        end
        t1 = tic;  m = nan(numel(idx), nV);  ms = strings(numel(idx), 1);  bC = bgV(idx);  aC = alC(idx);  sC = sgC(idx);  tC = trC(idx);
        parfor (j = 1:numel(idx), opts.NWorkers)
            [m(j, :), ms(j)] = oneSim_local(tC(j), bC(j), aC(j), sC(j), FS, sigToMM, seed0 + 1000 * idx(j), fns, nV);
        end
        secs = toc(t1);  tmp = fullfile(ckDir, sprintf('chunk_%03d_tmp.mat', c));
        save(tmp, 'm', 'ms', 'idx', 'sg', 'secs', '-v7.3');  movefile(tmp, ck);
        Med(idx, :) = m;  msg(idx) = ms;  chunkSecs(end + 1) = secs; %#ok<AGROW>
        fprintf('  chunk %d/%d done in %.1f min | ETA %.0f min | failed so far %d\n', c, nCh, secs / 60, mean(chunkSecs) * (nCh - c) / 60, nnz(msg ~= ""));
    end
    nFail = nnz(msg ~= "");
    if nFail > 0, fprintf('WARNING: %d of %d trials FAILED and are excluded (first: %s)\n', nFail, N, msg(find(msg ~= "", 1))); end
    if nFail / N > opts.MaxFailFrac, error('tea:failures', '%d of %d trials failed (limit %.0f%%). First: %s', nFail, N, 100 * opts.MaxFailFrac, msg(find(msg ~= "", 1))); end
    ok = all(isfinite(Med), 2);  hitM = abs(Med - bgV) <= opts.MDC;  idxOk = find(ok);  bsel = boot_local(subjSel, idxOk, opts.NBoot, opts.Seed + 5);
    hb = @(a) cellfun(@(p) mean(hitM(p, a)), bsel);
    C = {};   % contrasts: name, arm a, arm b  (a - b)
    for s = [1 3 5 7], if s + 1 <= nV, C(end+1, :) = {"anchoring effect " + labs(s), s + 1, s}; end, end %#ok<AGROW>   % default 8-arm layout; fewer arms give fewer contrasts
    for s = [3 5 7], if s + 1 <= nV, C(end+1, :) = {"window effect, production: " + labs(s) + " vs [4,4]", s, 1}; end, end %#ok<AGROW>
    for s = [3 5 7], if s + 1 <= nV, C(end+1, :) = {"window effect, anchored: " + labs(s + 1) + " vs [4,4]+EA", s + 1, 2}; end, end %#ok<AGROW>
    CT = table(strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', ["contrast" "diffPP" "ciLo" "ciHi"]);
    fprintf('%s (%d scored of %d)\n  %-14s %6s\n', d, nnz(ok), N, 'arm', 'hit %');
    for a = 1:nV, fprintf('  %-14s %6.1f\n', labs(a), 100 * mean(hitM(idxOk, a))); end
    fprintf('  %-46s %8s  %s\n', 'paired contrast (pp)', 'diff', '[95% CI]');
    for i = 1:size(C, 1)
        dv = 100 * (hb(C{i, 2}) - hb(C{i, 3}));  ci = prctile(dv, [2.5 97.5]);  dp = 100 * (mean(hitM(idxOk, C{i, 2})) - mean(hitM(idxOk, C{i, 3})));
        CT = [CT; table(string(C{i, 1}), dp, ci(1), ci(2), 'VariableNames', CT.Properties.VariableNames)]; %#ok<AGROW>
        fprintf('  %-46s %+8.1f  [%+.1f, %+.1f]\n', C{i, 1}, dp, ci);
    end
    R1(end + 1) = struct('dataset', d, 'n', nnz(ok), 'hit', 100 * mean(hitM(idxOk, :), 1), 'labels', labs, 'contrast', CT, 'betaGen', bgV, 'med', Med, 'subj', subjSel); %#ok<AGROW>
end

%% PART 2 -- real trials, production [4,4] with and without anchoring
for d = opts.RealDatasets
    [FS, sigToMM, importFn] = datasetReal_local(d);
    tr = importFn();  inV015 = load(fullfile(srcDir, "loopClosureResults_" + d + "_all_shaped_xu_v015.mat"), "results").results;
    tr = tr(ismember(string({tr.trialID}), string({inV015.trialID})));   % the headline corpus: the trials in the v015 results
    if numel(tr) ~= numel(inV015), error('tea:corpus', '%s: %d imported trials match %d v015 results', d, numel(tr), numel(inV015)); end
    nT = numel(tr);  rng(opts.Seed + 77, 'twister');  sel = sort(randperm(nT, min(opts.NReal, nT)))';  N = numel(sel);
    trC = arrayfun(@(i) struct('x', double(tr(i).x(:)), 'y', double(tr(i).y(:))), sel);  subj = string(arrayfun(@(i) subjectOf_local(tr(i)), sel));  tid = string({tr(sel).trialID})';  clear tr
    fprintf('\nPART 2 %s: %d random trials of %d | FS %d Hz\n', d, N, nT, FS);
    sg = struct('seed', opts.Seed, 'n', N, 'chunk', opts.Chunk, 'tid', tid, 'ds', d);
    ckDir = fullfile(resDir, "tea1_real_" + d + opts.Tag + "_chunks");  if ~isfolder(ckDir), mkdir(ckDir); end
    nCh = ceil(N / opts.Chunk);  chunkSecs = [];  BG = nan(N, 2);  INV = false(N, 2);  Z = nan(N, 1);  msg = strings(N, 1);  seed0 = opts.Seed;  fP = fns{1};  fE = fns{2};
    for c = 1:nCh
        idx = (c - 1) * opts.Chunk + 1 : min(c * opts.Chunk, N);  ck = fullfile(ckDir, sprintf('chunk_%03d.mat', c));
        if isfile(ck)
            K = load(ck);  if ~isequal(K.sg, sg) || ~isequal(K.idx, idx), error('tea:chunkMismatch', 'Checkpoint %s was written with different options; move %s away.', ck, ckDir); end
            BG(idx, :) = K.bg;  INV(idx, :) = K.iv;  Z(idx) = K.z;  msg(idx) = K.ms;  chunkSecs(end + 1) = K.secs; %#ok<AGROW>
            fprintf('  chunk %d/%d loaded (took %.1f min)\n', c, nCh, K.secs / 60);  continue
        end
        if ~isnat(opts.Deadline) && datetime('now') + seconds(1.15 * max([chunkSecs 0])) > opts.Deadline
            fprintf('\ntea PAUSED before real %s chunk %d/%d (not an error). Run again to continue.\n', d, c, nCh);  R = [];  return
        end
        t1 = tic;  bg = nan(numel(idx), 2);  iv = false(numel(idx), 2);  z = nan(numel(idx), 1);  ms = strings(numel(idx), 1);  tC = trC(idx);
        parfor (j = 1:numel(idx), opts.NWorkers)
            [bg(j, :), iv(j, :), z(j), ms(j)] = oneReal_local(tC(j), FS, sigToMM, seed0 + 1000 * idx(j), fP, fE);
        end
        secs = toc(t1);  tmp = fullfile(ckDir, sprintf('chunk_%03d_tmp.mat', c));
        save(tmp, 'bg', 'iv', 'z', 'ms', 'idx', 'sg', 'secs', '-v7.3');  movefile(tmp, ck);
        BG(idx, :) = bg;  INV(idx, :) = iv;  Z(idx) = z;  msg(idx) = ms;  chunkSecs(end + 1) = secs; %#ok<AGROW>
        fprintf('  chunk %d/%d done in %.1f min | ETA %.0f min | failed so far %d\n', c, nCh, secs / 60, mean(chunkSecs) * (nCh - c) / 60, nnz(msg ~= ""));
    end
    nFail = nnz(msg ~= "");
    if nFail > 0, fprintf('WARNING: %d of %d trials FAILED and are excluded (first: %s)\n', nFail, N, msg(find(msg ~= "", 1))); end
    if nFail / N > opts.MaxFailFrac, error('tea:failures', '%d of %d real trials failed. First: %s', nFail, N, msg(find(msg ~= "", 1))); end
    good = msg == "" & isfinite(Z);  isT = Z > opts.DzTol;
    tab = table();
    for g = ["tail" "non-tail" "all"]
        k = good & ((g == "tail" & isT) | (g == "non-tail" & ~isT) | g == "all");  ik = find(k);
        dB = BG(:, 2) - BG(:, 1);  fin = ik(isfinite(dB(ik)));
        ci = [NaN NaN];   % an empty group is reported as NaN, not dropped
        if ~isempty(fin), bs = boot_local(subj, fin, opts.NBoot, opts.Seed + 9);  bm = cellfun(@(p) median(dB(p)), bs);  ci = prctile(bm, [2.5 97.5]); end
        tab = [tab; table(g, numel(ik), numel(fin), median(dB(fin)), ci(1), ci(2), 100 * mean(INV(ik, 1)), 100 * mean(INV(ik, 2)), median(BG(ik, 1), 'omitnan'), median(BG(ik, 2), 'omitnan'), ...
            'VariableNames', ["group" "n" "nBothFinite" "medDeltaBG" "ciLo" "ciHi" "pctInvProd" "pctInvEA" "medBGprod" "medBGea"])]; %#ok<AGROW>
    end
    fprintf('%s: tail %d of %d (%.0f%%). Anchored - production, beta_gen* (paired):\n', d, nnz(good & isT), nnz(good), 100 * mean(isT(good)));  disp(tab)
    R2(end + 1) = struct('dataset', d, 'n', nnz(good), 'nTail', nnz(good & isT), 'tab', tab, 'bgProd', BG(:, 1), 'bgEA', BG(:, 2), 'isTail', isT, 'subj', subj, 'invProd', INV(:, 1), 'invEA', INV(:, 2)); %#ok<AGROW>
end
R = struct('part1', {R1}, 'part2', {R2});
save(outFile, 'R', 'R1', 'R2', 'opts');  fprintf('\nSaved %s\n=== COMPLETED OK ===\n', outFile);
end

function [m, ms] = oneSim_local(tr, bg, al, sg, FS, sigToMM, seed, fns, nV)
    m = nan(1, nV);  ms = "";
    try
        rng(seed, 'twister');
        g = loopClosureCore_v004(tr, FS, sigToMM, GeometryOnly=true);  M = numel(tr.x);
        [xc, yc, nx, ny] = synthWorld_v001("base", g, M, FS, bg, al, sg);
        syn = struct('x', xc + nx, 'y', yc + ny);
        for v = 1:nV, rng(seed + 1, 'twister');  fv = fns{v};  o = fv(syn, FS, sigToMM, NoiseModel="shaped_xu", NR7=0);  m(v) = o.betaGenStarMed; end
    catch ME
        ms = string(ME.identifier) + ": " + string(ME.message);
    end
end

function [bg, iv, z, ms] = oneReal_local(tr, FS, sigToMM, seed, fP, fE)
    bg = nan(1, 2);  iv = false(1, 2);  z = NaN;  ms = "";
    try
        g = loopClosureCore_v004(tr, FS, sigToMM, GeometryOnly=true);  z = zeroStretch_local(numel(tr.x), g.f0, FS, 4);
        for a = 1:2
            fv = {fP, fE};  rng(seed, 'twister');  o = fv{a}(tr, FS, sigToMM, NoiseModel="shaped_xu", NReps=20, NR7=0);
            bg(a) = o.betaGenStarMed;  iv(a) = isfinite(o.betaGenStarMed);
        end
    catch ME
        ms = string(ME.identifier) + ": " + string(ME.message);
    end
end

function z = zeroStretch_local(N, f0, fs, nCW)
% templateTailStratify_v002 edge_local for one trial: zero-residual stretch inside the analysed segment, in cycles.
    z = NaN;  if ~isfinite(f0) || f0 <= 0, return, end
    win0 = round(nCW / f0 * fs);  hop = round(win0 / 2);  win = max(min(win0, N), round(2 / f0 * fs));
    cl = max(round(win / 4), 1);  cs = cl;  ce = min(max(N - cl, cs + 100), N);
    if win > N, lastEnd = 0; else, lastEnd = (win - 1) + hop * floor((N - win) / hop) + 1; end
    z = max(0, ce - max(lastEnd, cs - 1)) * f0 / fs;
end

function bsel = boot_local(subj, idxOk, nBoot, seed)
% Subject-clustered bootstrap index sets over the scored rows.
    us = unique(subj(idxOk));  ii = cell(numel(us), 1);  for s = 1:numel(us), ii{s} = idxOk(subj(idxOk) == us(s)); end
    rng(seed);  bsel = cell(nBoot, 1);  for b = 1:nBoot, bsel{b} = vertcat(ii{randi(numel(us), numel(us), 1)}); end
end

function [fns, d] = makeVariants_local(srcDir, V)
% Text-patched copies of the CURRENT core. The single template call gets [nH nCW]; anchored variants call templateSubtractEA_local (appended below), which
% is templateSubtract_local plus one final window s = N - win + 1 when the regular windows stop short of the last sample. Unique names + path check (Addendum 25).
    txt = fileread(fullfile(srcDir, 'loopClosureCore_v004.m'));
    anchor = 'double(tr.x), double(tr.y), FS, res.f0, 4, 4);';  callHead = 'ceIdx] = templateSubtract_local( ...';
    if numel(strfind(txt, anchor)) ~= 1 || numel(strfind(txt, callHead)) ~= 1, error('tea:anchor', 'Template call not found exactly once in loopClosureCore_v004.m'); end
    d = tempname;  mkdir(d);  addpath(d);  fns = cell(1, size(V, 1));  stamp = char(datetime('now', 'Format', 'yyyyMMddHHmmssSSS'));
    for v = 1:size(V, 1)
        nm = sprintf('loopClosureCore_tea%s_%d', stamp, v);
        t = strrep(txt, anchor, sprintf('double(tr.x), double(tr.y), FS, res.f0, %g, %g);', V(v, 1), V(v, 2)));
        if V(v, 3) ~= 0, t = strrep(t, callHead, 'ceIdx] = templateSubtractEA_local( ...');  t = [t newline eaFunction_local()]; end %#ok<AGROW>
        t = strrep(t, 'function out = loopClosureCore_v004(', ['function out = ' nm '(']);
        fid = fopen(fullfile(d, [nm '.m']), 'w');  fwrite(fid, t);  fclose(fid);  fns{v} = str2func(nm);
        if ~startsWith(string(which(nm)), string(d)), error('tea:shadow', 'Patched core %s resolves to %s, not the run folder %s', nm, which(nm), d); end
    end
end

function s = eaFunction_local()
    L = {
    'function [resX, resY, fitX, fitY, cs, ce] = templateSubtractEA_local(x, y, fs, f0, nH, nCW)'
    '    N   = numel(x);'
    '    win = round(nCW/f0*fs);'
    '    hop = round(win/2);'
    '    win = min(win, N);'
    '    win = max(win, round(2/f0*fs));'
    '    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);'
    '    starts = 1:hop:(N - win + 1);'
    '    if win > N, starts = []; end   % no window fits (trace shorter than the minimum window): as production, residual stays 0'
    '    if ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % end-anchored final window'
    '    for s = starts'
    '        e   = s + win - 1;'
    '        tw  = (0:win-1)''/fs;'
    '        han = 0.5*(1 - cos(2*pi*(0:win-1)''/(win-1)));'
    '        D   = [ones(win,1), tw];'
    '        for h = 1:nH'
    '            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>'
    '            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>'
    '        end'
    '        Dw        = D .* han;'
    '        bX        = Dw \ (x(s:e) .* han);'
    '        bY        = Dw \ (y(s:e) .* han);'
    '        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;'
    '        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;'
    '        ww(s:e)   = ww(s:e) + han;'
    '    end'
    '    ok = ww > 0;'
    '    resX(ok) = resX(ok) ./ ww(ok);'
    '    resY(ok) = resY(ok) ./ ww(ok);'
    '    cl   = max(round(win/4), 1);'
    '    cs   = cl;'
    '    ce   = min(max(N - cl, cs + 100), N);'
    '    fitX = x(cs:ce) - resX(cs:ce);'
    '    fitY = y(cs:ce) - resY(cs:ce);'
    '    resX = resX(cs:ce);'
    '    resY = resY(cs:ce);'
    'end'};
    s = strjoin(L, newline);
end

function [matFile, sigToMM, FS, importFn] = datasetSim_local(name, srcDir)
    switch name
        case "Fraser",       matFile = 'noiseCharacterisation_fraser.mat';       sigToMM = 1/10.41793; FS = 240; importFn = @() importDB_fraser_v001(Verbose=false);
        case "Hickman PLAC", matFile = 'noiseCharacterisation_hickmanPLAC.mat';  sigToMM = 0.248;      FS = 133; importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Zarandi",      matFile = 'noiseCharacterisation_zarandi.mat';      sigToMM = 10.0;       FS = 100; importFn = @() importDB_zarandi_v001(Verbose=false);
        otherwise, error('tea:dataset', 'Dataset "%s" not supported', name);
    end
    matFile = fullfile(srcDir, matFile);
end

function [FS, sigToMM, importFn] = datasetReal_local(name)
    switch name
        case "Fraser",       FS = 240; sigToMM = 1/10.41793; importFn = @() importDB_fraser_v001(Verbose=false);
        case "Cook_CTRL",    FS = 133; sigToMM = 0.248;      importFn = @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook_ASD",     FS = 133; sigToMM = 0.248;      importFn = @() importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false);
        case "Hickman_PLAC", FS = 133; sigToMM = 0.248;      importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman_HALO", FS = 133; sigToMM = 0.248;      importFn = @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        otherwise, error('tea:dataset', 'Dataset "%s" not supported', name);
    end
end

function s = subjectOf_local(t)
    if isfield(t, 'subjectID'), s = string(t.subjectID); else, s = extractBefore(string(t.trialID), "_"); end
end
