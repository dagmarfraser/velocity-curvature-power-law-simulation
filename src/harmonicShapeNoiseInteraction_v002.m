function out = harmonicShapeNoiseInteraction_v002(opts)
%HARMONICSHAPENOISEINTERACTION_V002 How the end-matched noise step interacts with the number of harmonics stripped, K = 0 included (2026-11-20).
%   v002 = v001 (2026-11-18) with ONE purpose: K = 0 is read from the v033 nH0 files, and every input is checked before anything is computed.
%   The analysis code (classify_local, residualAlpha_local, oneTrial_local, estimateF0_local, templateSubtract_local) is v001's, unchanged.
%   What changed versus v001:
%     1. Preflight: every dataset x K file (K = opts.Ks, the primary arm K = 4 and fc) must exist, or the run stops with one error listing ALL
%        missing files (hsni:missingInputs). v001 marked a missing rung "missing" and carried on; v002 does not.
%     2. A rung file whose trial list differs from the primary arm's stops the run (hsni:partial). v001 marked it "partial" and excluded it.
%     3. K = 0 files are checked against the md5 recorded in the 2026-11-20 manifest (section D5); a mismatch stops the run (hsni:md5).
%     4. The table has an md5 column and a full-path file column for EVERY row, so each K level of each dataset is traceable to one file.
%     5. opts.OutDir and opts.SrcDir (defaults: ../results and the script's own folder, as v001), opts.K0Dir, opts.K0Override.
%   Per dataset and per K (0, 1, 2, 3, 4, 6, 10, fc) the table reports, from the saved v028/v033 results:
%     ran, failure types (noSegment, short, emShort, segAlpha, notReached, other), drawReject, segAlpha (median, quartiles), Lfrac, and
%     resAlpha (alphaPWLW over [1, min(20, fs/2 - 1)] Hz of the end-anchored template residual with K harmonics stripped, all trials).
%   Files: K = 1, 2, 3, 6, 10 as loopClosureResults_<dataset>_all_endMatched_segdetect_nH<K>_v028.mat in opts.LadderDir; K = 4 is the primary arm
%   (..._segdetect_v028.mat) and fc is ..._nHfc_v028.mat, both in opts.SrcDir; K = 0 is ..._nH0_v033.mat in opts.K0Dir, or the path given for that
%   dataset in opts.K0Override ("Dataset=/full/path.mat"), used where a dataset's proper-named file is not in K0Dir.
%   out = harmonicShapeNoiseInteraction_v002(LadderDir=".../scratch_K", K0Dir=".../outputs_K0", OutDir=".../outputs_harmonic", Tag="");
arguments
    opts.LadderDir (1,1) string = ""
    opts.K0Dir (1,1) string = ""
    opts.K0Override (1,:) string = string.empty(1, 0)
    opts.SrcDir (1,1) string = ""
    opts.OutDir (1,1) string = ""
    opts.Datasets (1,:) string = ["FraserN48" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"]
    opts.Ks (1,:) double = [0 1 2 3 4 6 10]
    opts.Tag (1,1) string = ""
    opts.ResidualAlpha (1,1) logical = true
    opts.UseParfor (1,1) logical = true
end
here = fileparts(mfilename('fullpath'));
srcDir = opts.SrcDir;  if srcDir == "", srcDir = string(here); end
addpath(char(srcDir), genpath(fullfile(char(srcDir), 'functions')));
resDir = opts.OutDir;  if resDir == "", resDir = string(fullfile(fileparts(here), 'results')); end
if ~isfolder(resDir), mkdir(char(resDir)); end
ws = warning('off', 'all');  wCleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
rungs = [string(opts.Ks), "fc"];
K0md5 = k0Md5_local();
ov = parseOverride_local(opts.K0Override);
% ---- preflight: every input must exist (and K = 0 must match its recorded md5) before anything is computed ----
problems = string.empty(0, 1);
for d = opts.Datasets
    problems = [problems; checkInput_local("primary", rungPath_local(srcDir, opts, ov, d, "4"), "")]; %#ok<AGROW>
    for r = 1:numel(rungs)
        K = rungs(r);  if K == "4", continue, end
        exp0 = "";  if K == "0", exp0 = lookupMd5_local(K0md5, d); end
        problems = [problems; checkInput_local("K=" + K, rungPath_local(srcDir, opts, ov, d, K), exp0)]; %#ok<AGROW>
    end
end
if ~isempty(problems)
    error('hsni:missingInputs', 'harmonicShapeNoiseInteraction_v002 stopped before computing anything; %d input problem(s):\n%s', numel(problems), strjoin(problems, newline));
end
rows = table();
for d = opts.Datasets
    tag = d;  nice = strrep(d, "_", " ");
    P = loadPrimary_local(rungPath_local(srcDir, opts, ov, d, "4"));   % the K = 4 primary arm: also the reference trial list
    ids = string({P.res.trialID});  n = numel(ids);
    resAlpha = NaN(n, numel(rungs));
    if opts.ResidualAlpha, resAlpha = residualAlpha_local(d, ids, opts.Ks, opts.UseParfor, rungs); end
    for r = 1:numel(rungs)
        K = rungs(r);
        f = rungPath_local(srcDir, opts, ov, d, K);
        [R, cfg, state] = loadRung_local(f, K, P.res, P.cfg, tag);
        row = struct('dataset', nice, 'K', K, 'state', state, 'n', n, 'ran', NaN, 'noSegment', NaN, 'short', NaN, 'emShort', NaN, 'segAlpha', NaN, 'notReached', NaN, 'other', NaN, ...
            'drawReject', NaN, 'segAlphaMed', NaN, 'segAlphaQ1', NaN, 'segAlphaQ3', NaN, 'LfracMed', NaN, 'resAlphaMed', median(resAlpha(:, r), 'omitnan'), ...
            'resAlphaQ1', prctile(resAlpha(:, r), 25), 'resAlphaQ3', prctile(resAlpha(:, r), 75), 'file', f, 'md5', md5_local(f));
        st = string([R.emStatus]);  ok = st == "ok";  cat = classify_local(st);
        row.ran = mean(ok);  row.noSegment = mean(cat == "noSegment");  row.short = mean(cat == "short");  row.emShort = mean(cat == "emShort");
        row.segAlpha = mean(cat == "segAlpha");  row.notReached = mean(cat == "notReached");  row.other = mean(cat == "other");
        a = cell2mat(arrayfun(@(x) mean(x.emAlphaSeg, 'omitnan'), R(ok), 'UniformOutput', false));
        if any(ok)
            row.segAlphaMed = median(a, 'omitnan');  row.segAlphaQ1 = prctile(a, 25);  row.segAlphaQ3 = prctile(a, 75);
            row.LfracMed = median([R(ok).emLfrac], 'omitnan');
            nd = cfg.N_BETA * cfg.N_REPS;  row.drawReject = mean([R(ok).nNoiseRejected] / nd, 'omitnan');
        end
        rows = [rows; struct2table(row)]; %#ok<AGROW>
    end
    fprintf('%-13s done (%d trials)\n', nice, n);
end
out = struct('table', rows, 'opts', opts, 'generated', string(datetime('now')));
base = fullfile(resDir, "harmonicShapeNoiseInteraction_v002" + opts.Tag);
writetable(rows, base + ".csv");  save(base + ".mat", '-struct', 'out', '-v7.3');
disp(rows(:, {'dataset', 'K', 'state', 'ran', 'noSegment', 'emShort', 'segAlpha', 'notReached', 'drawReject', 'segAlphaMed', 'resAlphaMed'}));
plot_local(rows, base + ".png");
fprintf('Saved %s.{csv,mat,png}\n', base);
end

function m = k0Md5_local()   % md5 of the seven K = 0 (v033, nH0) files, Manifest_BRM_DRAFT_20261120_v002.md section D5, re-measured 2026-11-20
m = struct('FraserN48', "dee84a0aa1548ab8a9d6833c7ad2f3de", 'Zarandi', "7f042dffaf1548659621c9a9d6fe4f8f", 'Dhieb', "2e13d3aba36ac9d360c276c60fe6a200", ...
    'Hickman_PLAC', "938f024545dd6ce380efe4ad10744a22", 'Hickman_HALO', "54d8295b8ebcc675f3749ab9de1d7d5c", ...
    'Cook_CTRL', "79926c183d3f0611b0dafaed8d441255", 'Cook_ASD', "6fcd9582456780f75fa1e531a09cd4ec");
end

function h = lookupMd5_local(m, d)
if ~isfield(m, char(d)), error('hsni:md5', 'no recorded K = 0 md5 for dataset %s', d); end
h = m.(char(d));
end

function ov = parseOverride_local(s)
ov = struct();
for i = 1:numel(s)
    p = split(s(i), "=");  if numel(p) < 2, error('hsni:override', 'K0Override entries must be "Dataset=/path": %s', s(i)); end
    ov.(char(p(1))) = strjoin(p(2:end), "=");
end
end

function f = rungPath_local(srcDir, opts, ov, d, K)
base = "loopClosureResults_" + d + "_all_endMatched_segdetect_";
if K == "4", f = fullfile(srcDir, base + "v028.mat");
elseif K == "fc", f = fullfile(srcDir, base + "nHfc_v028.mat");
elseif K == "0"
    if isfield(ov, char(d)), f = ov.(char(d)); else, f = fullfile(opts.K0Dir, base + "nH0_v033.mat"); end
else, f = fullfile(opts.LadderDir, base + "nH" + K + "_v028.mat");
end
f = string(f);
end

function p = checkInput_local(label, f, expMd5)
p = string.empty(0, 1);
if ~isfile(f), p = label + " MISSING: " + f;  return, end
if expMd5 ~= ""
    h = md5_local(f);
    if h ~= expMd5, p = label + " md5 MISMATCH: " + f + " has " + h + ", expected " + expMd5; end
end
end

function h = md5_local(f)
md = java.security.MessageDigest.getInstance('MD5');  fid = fopen(f, 'r');  c = onCleanup(@() fclose(fid));
while ~feof(fid), b = fread(fid, 2^24, '*uint8');  if ~isempty(b), md.update(typecast(b, 'int8')); end, end
h = string(lower(reshape(dec2hex(typecast(md.digest, 'uint8'))', 1, [])));
end

function P = loadPrimary_local(f)
S = load(f, 'results', 'config');  P = struct('res', S.results, 'cfg', S.config);
end

function [R, cfg, state] = loadRung_local(f, K, ref, refCfg, tag)
if K == "4", R = ref;  cfg = refCfg;  state = "complete";  return, end
S = load(f, 'results', 'config');
if numel(S.results) ~= numel(ref) || ~isequal(string({S.results.trialID}), string({ref.trialID}))
    error('hsni:partial', '%s K=%s: trial list differs from the primary arm (%d against %d trials): %s', tag, K, numel(S.results), numel(ref), f);
end
R = S.results;  cfg = S.config;  state = "complete";
end

function cat = classify_local(st)
cat = repmat("other", size(st));
cat(st == "ok") = "ok";
cat(contains(st, "no length of")) = "noSegment";
cat(contains(st, "shorter than MinLength")) = "short";
cat(contains(st, "below the runner minimum")) = "emShort";
cat(contains(st, "segment alpha not finite")) = "segAlpha";
cat(startsWith(st, "not reached")) = "notReached";
end

function A = residualAlpha_local(d, ids, Ks, usePar, rungs)
switch d
    case "FraserN48",      tr = importDB_fraser_v002(Verbose=false);
    case "Cook_CTRL",      tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
    case "Cook_ASD",       tr = importDB_cook_v002(Group="ASD",  Tasks=7, Verbose=false);
    case "Hickman_PLAC",   tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
    case "Hickman_HALO",   tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
    case "Zarandi",        tr = importDB_zarandi_v001(Verbose=false);
    case "Dhieb",          tr = importDB_dhieb_v001(Verbose=false);
end
[tf, loc] = ismember(ids, string({tr.trialID}));  if ~all(tf), error('hsni:ids', '%s: %d result trials not in the importer', d, nnz(~tf)); end
tr = tr(loc);  n = numel(tr);  nK = numel(Ks);  A = NaN(n, numel(rungs));
X = {tr.x};  Y = {tr.y};  FS = [tr.fs];
cols = NaN(n, nK);
if usePar
    try, ensurePool_v001(); catch, end   %#ok<TRYNC>
    parfor i = 1:n, cols(i, :) = oneTrial_local(double(X{i}), double(Y{i}), FS(i), Ks); end
else
    for i = 1:n, cols(i, :) = oneTrial_local(double(X{i}), double(Y{i}), FS(i), Ks); end
end
A(:, 1:nK) = cols;
end

function a = oneTrial_local(x, y, fs, Ks)
a = NaN(1, numel(Ks));  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
band = [1.0, min(20.0, fs/2 - 1)];
for k = 1:numel(Ks)
    try
        [rx, ry] = templateSubtract_local(x, y, fs, f0, Ks(k), 4);
        a(k) = mean([alphaPWLW_v001(rx(:) - mean(rx), fs, band), alphaPWLW_v001(ry(:) - mean(ry), fs, band)]);
    catch
    end
end
end

function f0 = estimateF0_local(x, y, fs)   % as runLoopClosureFftnoise_v028
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

function [resX, resY] = templateSubtract_local(x, y, fs, f0, nH, nCW)   % as runLoopClosureFftnoise_v028, END-ANCHORED
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);
    if win > N, starts = []; end
    if ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end
    for s = starts
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e) .* han);
        bY        = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
        ww(s:e)   = ww(s:e) + han;
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    ce   = min(max(N - cl, cl + 100), N);
    resX = resX(cl:ce);
    resY = resY(cl:ce);
end

function plot_local(T, f)
fig = figure('Visible', 'off', 'Position', [50 50 1400 520], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Harmonics stripped (K) against the end-matched noise step, per dataset', 'FontWeight', 'bold');
ds = unique(T.dataset, 'stable');  Ks = unique(T.K, 'stable');  x = 1:numel(Ks);
panels = {'ran', 'share of trials with an end-matched segment', [0 1]; 'resAlphaMed', 'median exponent of the template residual (all trials)', []; 'drawReject', 'share of noise draws rejected (trials that ran)', []};
for p = 1:3
    nexttile; hold on;
    for i = 1:numel(ds)
        v = NaN(1, numel(Ks));
        for k = 1:numel(Ks), m = T.dataset == ds(i) & T.K == Ks(k) & T.state == "complete";  if any(m), v(k) = T.(panels{p, 1})(m); end, end
        plot(x, v, '-o', 'LineWidth', 1.4, 'DisplayName', char(ds(i)));
    end
    set(gca, 'XTick', x, 'XTickLabel', cellstr(Ks)); xlabel('K stripped'); ylabel(panels{p, 2}); grid on;
    if ~isempty(panels{p, 3}), ylim(panels{p, 3}); end
    if p == 1, legend('Location', 'southwest', 'FontSize', 7); end
end
exportgraphics(fig, f, 'Resolution', 150);  close(fig);
end
