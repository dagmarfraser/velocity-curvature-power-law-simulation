%% analyseAlphaRateLadder_v001.m
% Read-only analysis of results/obsAlphaRateLadder_v001_<corpus>_stage<k>.mat (docs/SPEC_AlphaRateLadder_v001.md section 7).
%   Alpha: per-rung medians with person-cluster bootstrap CIs for arms A4, A12, B; the model comparison by leave-one-person-out
%          RMSE (log2 fs, fs, 1/fs, quadratic in log2 fs); local slopes between adjacent rungs; slopes over fixed windows;
%          per-shape slopes when the file holds several shapes.
%   Beta_obs: paired change against 240 Hz for the BW, SGp and SGw pipelines, the share of trials past the SEM, and the split by
%          the in-domain rule f0 <= fs / 11 (cutoff at least 5.5 f0 and below Nyquist; not tested for shapes with angular frequency > 2).
% Reads one results file; writes results/analyseAlphaRateLadder_v001_<corpus>_stage<k>.mat and two figures in figures/. Run from src/.

%% CONFIG
CORPUS = "fraser";  STAGE = 0;
FS_MIN = 8;  FS_MIN_SENS = 30;  N_BOOT = 2000;  SEED = 20261102;  SEM = 0.0108;  RULE_FS = 11;  FS0 = 240;
ROOT = fileparts(fileparts(mfilename("fullpath")));  TAG = sprintf("%s_stage%d", CORPUS, STAGE);
IN_MAT  = fullfile(ROOT, "results", "obsAlphaRateLadder_v001_" + TAG + ".mat");
OUT_MAT = fullfile(ROOT, "results", "analyseAlphaRateLadder_v001_" + TAG + ".mat");
if ~isfile(IN_MAT), error("ladderA:input", "%s", "Missing: " + IN_MAT); end
R = load(IN_MAT);  S = R.S;  nT = numel(R.ids);
[uD, ~, ~] = unique(S.D);  fr = FS0 ./ uD;  [fr, o] = sort(fr);  uD = uD(o);  nR = numel(fr);   % rungs, ascending rate
[~, ~, pIdx] = unique(R.person);  nP = max(pIdx);
fprintf("%s: %d trials, %d persons, %d rungs (%s Hz)\n", TAG, nT, nP, nR, num2str(fr'));
if ~any(fr == FS0), error("ladderA:no240", "%s", "The 240 Hz rung is needed as the reference"); end

%% Per-rung arrays (mean over phases within a trial)
AL  = cat(3, byRung_local(R.aA, S, uD), byRung_local(R.aB, S, uD));            % trials x rungs x [A4 A12 B]
ARM = ["A4" "A12" "B"];
BR  = byRung_local(R.bet, S, uD);                                              % trials x rungs x 9
dB  = BR - BR(:, fr == FS0, :);                                                % paired change against 240 Hz
PIPE = R.PIPE;

%% Alpha: medians, CIs, windows, models
bA   = bootMed_local(AL, pIdx, N_BOOT, SEED);
medA = squeeze(median(AL, 1, "omitnan"));  loA = squeeze(prctile(bA, 2.5, 1));  hiA = squeeze(prctile(bA, 97.5, 1));
T1 = table(fr, medA(:, 1), medA(:, 2), medA(:, 3), loA(:, 1), hiA(:, 1), 'VariableNames', ["fs" "A4" "A12" "B" "A4lo" "A4hi"]);
fprintf("\nMedian alpha per rung (person-cluster bootstrap CI for A4):\n"); disp(T1);

W = struct("name", {"fs >= 8" "fs >= 30" "60 to 240" "60, 120, 240 only"}, ...
    "mask", {fr >= FS_MIN, fr >= FS_MIN_SENS, fr >= 60, ismember(fr, [60 120 240])});
T2 = table();
for w = 1:numel(W)
    m = W(w).mask;  for a = 1:3
        sl = slope_local(fr(m), medA(m, a));  sb = NaN(N_BOOT, 1);
        for b = 1:N_BOOT, sb(b) = slope_local(fr(m), squeeze(bA(b, m, a))'); end
        T2 = [T2; table(string(W(w).name), ARM(a), sum(m), sl, prctile(sb, 2.5), prctile(sb, 97.5), 'VariableNames', ["window" "arm" "nRungs" "slopePerOctave" "lo" "hi"])]; %#ok<AGROW>
    end
end
fprintf("\nAlpha slope per octave of rate (dalpha / dlog2 fs):\n"); disp(T2);

prim = find(fr >= FS_MIN);  frp = fr(prim);
BAS = {@(f) [ones(size(f)) log2(f)], @(f) [ones(size(f)) f], @(f) [ones(size(f)) 1 ./ f], @(f) [ones(size(f)) log2(f) log2(f) .^ 2]};
MODEL = ["M1 log2 fs" "M2 fs" "M3 1/fs" "M4 quadratic in log2 fs"];
rm = lopoRMSE_local(AL(:, prim, 1), pIdx, frp(:), BAS);  T3 = table(MODEL', rm, rm ./ min(rm), 'VariableNames', ["model" "lopoRMSE" "ratioToBest"]);
T3.log2Adequate = repmat(rm(1) <= 1.05 * min(rm), numel(rm), 1);
fprintf("\nLeave-one-person-out RMSE of the per-rung median, arm A4, rungs >= %g Hz (rule: log2 adequate if M1 within 5%% of the best):\n", FS_MIN); disp(T3);

L = (1:numel(prim) - 1)';  lo2 = log2(frp);
dl = (squeeze(bA(:, prim(L + 1), 1)) - squeeze(bA(:, prim(L), 1))) ./ (lo2(L + 1) - lo2(L))';
df = (squeeze(bA(:, prim(L + 1), 1)) - squeeze(bA(:, prim(L), 1))) ./ (frp(L + 1) - frp(L))';
T4 = table(frp(L), frp(L + 1), (medA(prim(L + 1), 1) - medA(prim(L), 1)) ./ (lo2(L + 1) - lo2(L)), prctile(dl, 2.5, 1)', prctile(dl, 97.5, 1)', ...
    (medA(prim(L + 1), 1) - medA(prim(L), 1)) ./ (frp(L + 1) - frp(L)), prctile(df, 2.5, 1)', prctile(df, 97.5, 1)', ...
    'VariableNames', ["fsLow" "fsHigh" "dAlphaPerOctave" "octLo" "octHi" "dAlphaPerHz" "hzLo" "hzHi"]);
fprintf("\nLocal slopes between adjacent rungs, arm A4:\n"); disp(T4);

T5 = table();
if numel(unique(R.shapeCode(isfinite(R.shapeCode)))) > 1
    for sh = unique(R.shapeCode(isfinite(R.shapeCode)))'
        j = R.shapeCode == sh;  [~, ~, pj] = unique(R.person(j));  bj = bootMed_local(AL(j, :, 1), pj, N_BOOT, SEED);  mj = median(AL(j, :, 1), 1, "omitnan")';
        for w = 1:2
            m = W(w).mask;  sb = NaN(N_BOOT, 1);  for b = 1:N_BOOT, sb(b) = slope_local(fr(m), squeeze(bj(b, m))'); end
            T5 = [T5; table(sh, string(W(w).name), sum(j), slope_local(fr(m), mj(m)), prctile(sb, 2.5), prctile(sb, 97.5), 'VariableNames', ["shape" "window" "nTrials" "slope" "lo" "hi"])]; %#ok<AGROW>
        end
    end
    fprintf("\nAlpha slope per shape, arm A4:\n"); disp(T5);
end

%% Beta_obs
bB = bootMed_local(dB, pIdx, N_BOOT, SEED);  medB = squeeze(median(dB, 1, "omitnan"));
share = squeeze(mean(abs(dB) > SEM, 1, "omitnan"));
isD = false(nT, nR);  for r = 1:nR, isD(:, r) = R.f0 <= fr(r) / RULE_FS; end
sel = [3 6 9];  T6 = table();
for r = 1:nR
    for q = sel
        dq = dB(:, r, q);  T6 = [T6; table(fr(r), PIPE(q), medB(r, q), prctile(bB(:, r, q), 2.5), prctile(bB(:, r, q), 97.5), share(r, q), ...
            median(dq(isD(:, r)), "omitnan"), median(dq(~isD(:, r)), "omitnan"), mean(isD(:, r)), ...
            'VariableNames', ["fs" "pipeline" "medianChange" "lo" "hi" "shareBeyondSEM" "medianInDomain" "medianOutDomain" "inDomainShare"])]; %#ok<AGROW>
    end
end
fprintf("\nPaired change in beta_obs against 240 Hz (IRLS; SEM %.4f); in-domain = f0 <= fs/%d:\n", SEM, RULE_FS); disp(T6);
fprintf("Rungs where the SG frame is floor-bound: %s Hz; rungs where the 10 Hz Butterworth is beyond Nyquist: %s Hz.\n", ...
    num2str(unique(S.fs(S.floorBound))'), num2str(unique(S.fs(~S.bwValid))'));

%% Figures
fig = figure("Visible", "off", "Position", [100 100 1100 420]);
for pnl = 1:2
    ax = subplot(1, 2, pnl);  hold(ax, "on");  col = lines(3);
    for a = 1:3
        if a == 1, fill([fr; flipud(fr)], [loA(:, 1); flipud(hiA(:, 1))], col(1, :), "FaceAlpha", 0.15, "EdgeColor", "none", "HandleVisibility", "off"); end
        plot(ax, fr, medA(:, a), "-o", "Color", col(a, :), "DisplayName", ARM(a), "MarkerSize", 4);
    end
    xline(ax, FS_MIN, ":k", "HandleVisibility", "off");  xlabel(ax, "sampling rate (Hz)");  ylabel(ax, "median residual-IRASA \alpha");
    if pnl == 2, set(ax, "XScale", "log"); xticks(ax, fr); xticklabels(ax, string(fr)); title(ax, "log axis"); else, title(ax, "linear axis"); end
    legend(ax, "Location", "best");
end
exportgraphics(fig, fullfile(ROOT, "figures", "alphaRateLadder_v001_" + TAG + "_alpha.png"), "Resolution", 200);  close(fig);
fig = figure("Visible", "off", "Position", [100 100 1100 420]);
ax = subplot(1, 2, 1);  hold(ax, "on");  col = lines(3);
for q = 1:3
    plot(ax, fr, medB(:, sel(q)), "-o", "Color", col(q, :), "DisplayName", PIPE(sel(q)), "MarkerSize", 4);
end
yline(ax, [-SEM SEM], ":k", "HandleVisibility", "off");  set(ax, "XScale", "log");  xticks(ax, fr);  xticklabels(ax, string(fr));
xlabel(ax, "sampling rate (Hz)");  ylabel(ax, "median change in \beta_{obs} against 240 Hz");  legend(ax, "Location", "best");
ax = subplot(1, 2, 2);  p1 = S(S.phase == 1, :);  [~, so] = sort(p1.fs);  p1 = p1(so, :);
plot(ax, p1.fs, p1.durProd, "-o", p1.fs, p1.durWM, "-s", "MarkerSize", 4);  set(ax, "XScale", "log");  xticks(ax, p1.fs);  xticklabels(ax, string(p1.fs));
xlabel(ax, "sampling rate (Hz)");  ylabel(ax, "SG window (s)");  legend(ax, ["SGp (case 6)" "SGw (case 4)"], "Location", "best");
exportgraphics(fig, fullfile(ROOT, "figures", "alphaRateLadder_v001_" + TAG + "_beta.png"), "Resolution", 200);  close(fig);

save(OUT_MAT, "fr", "ARM", "medA", "loA", "hiA", "T1", "T2", "T3", "T4", "T5", "T6", "medB", "share", "SEM", "RULE_FS", "FS_MIN", "N_BOOT", "SEED", "-v7.3");
fprintf("\nSaved: %s\n", OUT_MAT);

%% =========================================================================
function M = byRung_local(X, S, uD)
% Trials x schemes (x k) to trials x rungs (x k): the mean over a rung's phases within each trial.
    sz = size(X);  X = reshape(X, sz(1), sz(2), []);  M = NaN(sz(1), numel(uD), size(X, 3));
    for r = 1:numel(uD), M(:, r, :) = mean(X(:, S.D == uD(r), :), 2, "omitnan"); end
end

function B = bootMed_local(A, pIdx, nB, seed)
% Person-cluster bootstrap of the column medians of A (trials x rungs x arms); pIdx is the person index of each trial.
    rng(seed);  nP = max(pIdx);  L = cell(nP, 1);  for u = 1:nP, L{u} = find(pIdx == u); end
    B = NaN([nB, size(A, 2), size(A, 3)]);
    for b = 1:nB
        j = vertcat(L{randi(nP, nP, 1)});  B(b, :, :) = median(A(j, :, :), 1, "omitnan");
    end
end

function s = slope_local(f, a)
    g = isfinite(a(:)) & isfinite(f(:));  s = NaN;  if sum(g) < 2, return, end
    p = polyfit(log2(f(g)), a(g), 1);  s = p(1);
end

function rm = lopoRMSE_local(A, pIdx, fr, BAS)
% Leave-one-person-out: fit each model to the other persons' per-rung medians, score it on the held-out person's per-rung medians.
    nP = max(pIdx);  se = zeros(numel(BAS), 1);  n = 0;
    for u = 1:nP
        tr = pIdx ~= u;  mt = median(A(tr, :), 1, "omitnan");  mh = median(A(~tr, :), 1, "omitnan");  g = isfinite(mt) & isfinite(mh);
        for m = 1:numel(BAS)
            X = BAS{m}(fr(g));  beta = X \ mt(g)';  e = mh(g)' - X * beta;  se(m) = se(m) + sum(e .^ 2);
        end
        n = n + sum(g);
    end
    rm = sqrt(se / n);
end
