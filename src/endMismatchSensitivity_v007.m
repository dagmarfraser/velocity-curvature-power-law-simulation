function [R, S] = endMismatchSensitivity_v007(datasetName, opts)
%ENDMISMATCHSENSITIVITY_V007 Segment-length rule and full-length extension of the end-matched random-phase noise.
%   v007 (2026-11-15, coherence check): v006 with the production template's uncovered tail removed from the residual
%   (the zero stretch after the last overlap-add window; checkTemplateTailSegment_v001 found it inside 95-98% of the
%   segments of tailed trials). Everything else is v006, same seeds. v006 stays the record.
%   [R, S] = endMismatchSensitivity_v007("Hickman PLAC")   (also "FraserN48", "Hickman HALO", "Cook CTRL", "Cook ASD")
%   Follows Dagmar's decisions of 2026-11-14: end-matched random phases is the primary noise; the segment length must come
%   from the data, not a fixed floor; and extending the noise to the full trial length must be tested.
%   Residuals are rebuilt exactly as in endMismatchSensitivity_v002 (estimateF0_local and templateSubtract_local copied
%   verbatim from runLoopClosureFftnoise_v020; nH = 4, nCW = 4; major/minor rotation). All kinds use uniformly random
%   phases (the v002 "endToEndWhite" construction) and alpha selection (|dAlpha| <= Tol, at most MaxDraw draws):
%     e2e75      segment by endToEndSegment_v002 Rule "minFraction", MinFraction 0.75 (= v002 endToEndWhite)
%     e2e90      the same with MinFraction 0.90
%     e2eFair    segment by Rule "fairShare" (Tol 1): the longest segment whose wrap step is no larger than its mean step
%     extWrapFair  e2eFair surrogate (length L) continued circularly to the trial length M (u is periodic, so the joins
%                are ordinary steps; it repeats M - L samples)
%     extSpecFair  the fairShare segment's periodogram, interpolated linearly in frequency onto the M-point grid, with
%                random phases: a length-M surrogate with no repeated samples
%     extSpec75  the same from the 75% minFraction segment
%   Scoring, per surrogate (medians over NSurr): the e2e kinds against their segment (as v002), the ext kinds against the
%   whole residual r; every kind is ALSO scored against r (columns *Full). e2e kinds are rescaled to std(segment), ext
%   kinds to std(r). Alpha is selected against the segment's alpha (the spectrum's source) in every kind.
%     velRatio*/accRatio*: SD of the production-filtered first / second derivative, surrogate over reference
%       (differentiateKinematicsEBR, SG [4 17] and BWFD [2 10 1], edge clip 20; as v002)
%   Seeds (docs/README_SeedPractice_v001.md): threefry4x64_20, one master (opts.Master), one substream per
%   (dataset, trial, axis, kind, surrogate) from a fixed dataset list; the substream formula is written to the summary
%   file and the manifest CSV holds one row per cell. NSurr = 5 makes this a screening run, not a three-master result.
%   Writes results/endMismatchSensitivity/endMismatch_<dataset>_v007.csv, _summary_v007.csv, _manifest_v007.csv.
%   Every failure is in the status column and counted; nothing is dropped.
arguments
    datasetName (1,1) string
    opts.NSurr (1,1) double {mustBeInteger, mustBeInRange(opts.NSurr, 1, 99)} = 5
    opts.MaxDraw (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.Tol (1,1) double {mustBePositive} = 0.3
    opts.Master (1,1) double {mustBeInteger, mustBeNonnegative} = 20261114
    opts.OutDir (1,1) string = ""
    opts.MaxTrials (1,1) double {mustBePositive} = Inf     % smoke test: first MaxTrials trials only (output tagged _smoke)
end
DATASETS = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD"];     % fixed: indexes substreams
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
switch datasetName
    case "FraserN48",    FS = 240; importFn = @() importDB_fraser_v002(Verbose=false);
    case "Hickman PLAC", FS = 133; importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
    case "Hickman HALO", FS = 133; importFn = @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
    case "Cook CTRL",    FS = 133; importFn = @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
    case "Cook ASD",     FS = 133; importFn = @() importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false);
    otherwise, error('endMismatch:dataset', 'Unknown dataset "%s".', datasetName);
end
dIdx = find(DATASETS == datasetName);
for fn = ["alphaPWLW_v001", "differentiateKinematicsEBR", "endToEndSegment_v002", "ensurePool_v001"]
    if isempty(which(fn)), error('endMismatch:path', '%s not on the path.', fn); end
end
if opts.OutDir == "", opts.OutDir = fullfile(here, '..', 'results', 'endMismatchSensitivity'); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end
trials = importFn();  trials = trials(1:min(numel(trials), opts.MaxTrials));  nT = numel(trials);
band = [1 min(20, FS/2 - 1)];  nS = opts.NSurr;  maxDraw = opts.MaxDraw;  tol = opts.Tol;  master = opts.Master;
pool = ensurePool_v001();
wait(parfevalOnAll(pool, @addpath, 0, here, fullfile(here, 'functions')));
out = cell(nT, 1);  man = cell(nT, 1);
parfor t = 1:nT
    [out{t}, man{t}] = oneTrial(trials(t), t, dIdx, master, FS, band, nS, maxDraw, tol);
end
R = vertcat(out{:});  R.dataset = repmat(datasetName, height(R), 1);  MF = vertcat(man{:});  MF.dataset = repmat(datasetName, height(MF), 1);
tag = strrep(datasetName, " ", "");  if isfinite(opts.MaxTrials), tag = tag + "_smoke"; end
writetable(R, fullfile(opts.OutDir, "endMismatch_" + tag + "_v007.csv"));
writetable(MF, fullfile(opts.OutDir, "endMismatch_" + tag + "_manifest_v007.csv"));
ok = R.status == "ok";  g = ["axis", "kind"];
v = ["Lfrac", "cJump", "cSlip", "dAlphaMed", "velRatioSG", "accRatioSG", "velRatioBW", "accRatioBW", "velRatioBWFull", "accRatioBWFull", "nDrawMed"];
S = groupsummary(R(ok, :), g, @(q) median(q, 'omitnan'), v);
p = groupsummary(R(ok, :), g, {@(q) pct(q, 0.1), @(q) pct(q, 0.9)}, ["Lfrac", "velRatioBW", "accRatioBW", "velRatioBWFull", "accRatioBWFull"]);
S.Properties.VariableNames = replace(S.Properties.VariableNames, "fun1_", "med_");
p.Properties.VariableNames = replace(replace(p.Properties.VariableNames, "fun1_", "p10_"), "fun2_", "p90_");
S = join(S, removevars(p, "GroupCount"), 'Keys', g);
m = groupsummary(R(ok, :), g, "mean", ["shareNoAccept"]);
S = join(S, removevars(m, "GroupCount"), 'Keys', g);
S.master = repmat(master, height(S), 1);
S.substreamFormula = repmat("dIdx*1e9 + trial*1e4 + axis*1e3 + kind*100 + surrogate; DATASETS = " + strjoin(DATASETS, ","), height(S), 1);
writetable(S, fullfile(opts.OutDir, "endMismatch_" + tag + "_summary_v007.csv"));
fprintf('%s: fs %d, %d trials, median M %g, %d rows, %d with a FAIL status (%s)\n', datasetName, FS, nT, median(R.M), height(R), nnz(~ok), ...
    strjoin(unique(extractBefore(R.status(~ok) + ":", ":")), "; "));
disp(S);
if any(~ok), warning('endMismatch:fails', '%d of %d rows have a FAIL status; see the status column.', nnz(~ok), height(R)); end
end

function q = pct(v, p)
v = sort(v(~isnan(v)));  if isempty(v), q = NaN; else, q = v(min(numel(v), max(1, round(p * numel(v))))); end
end

function [R, MF] = oneTrial(tr, tIdx, dIdx, master, FS, band, nS, maxDraw, tol)
x = double(tr.x(:));  y = double(tr.y(:));
f0 = estimateF0_local(x, y, FS);
[xR, yR] = templateSubtract_local(x, y, FS, f0, 4, 4, false);
xc = x - mean(x);  yc = y - mean(y);
[V, D] = eig(cov(xc, yc));  [~, ord] = sort(diag(D), 'descend');  V = V(:, ord);  th = atan2(V(2,1), V(1,1));
res = {xR*cos(th) + yR*sin(th), -xR*sin(th) + yR*cos(th)};  ax = ["maj", "min"];
names = ["e2e75", "e2e90", "e2eFair", "extWrapFair", "extSpecFair", "extSpec75"];
segOf = ["m75", "m90", "fair", "fair", "fair", "m75"];  isExt = [false false false true true true];
nM = numel(names);  rows = cell(2, nM, 2);  MF = table();
for j = 1:2
    r = res{j}(:);  M = numel(r);  aR = NaN;  st = repmat("ok", 1, nM);
    X = NaN(nS, nM, 2, 7);                % [dAlpha velSG accSG velBW accBW velBWFull accBWFull]; 3rd dim: first draw, selected draw
    nD = NaN(nS, nM);  nA = NaN(nS, nM);  aRef = NaN(1, nM);  Lf = NaN(1, nM);  cJ = NaN(1, nM);  cS = NaN(1, nM);
    try
        if std(r) == 0, error('endMismatch:zeroResidual', 'residual is identically zero (trial shorter than 2 movement cycles)'); end
        aR = alphaPWLW_v001(r, FS, band);  frFull = filt_local(r, FS);
    catch ME
        st(:) = "FAIL residual: " + ME.message;
    end
    seg = struct();
    if all(st == "ok")
        for key = ["m75", "m90", "fair"]
            try
                switch key
                    case "m75",  [i0, i1, info] = endToEndSegment_v002(r, Rule="minFraction", MinFraction=0.75);
                    case "m90",  [i0, i1, info] = endToEndSegment_v002(r, Rule="minFraction", MinFraction=0.90);
                    case "fair", [i0, i1, info] = endToEndSegment_v002(r, Rule="fairShare", Tol=1);
                end
                seg.(key) = struct('ref', r(i0:i1), 'info', info, 'err', "");
            catch ME
                seg.(key) = struct('ref', [], 'info', [], 'err', ME.message);
            end
        end
    end
    for m = find(st == "ok")
        try
            sg = seg.(segOf(m));  if sg.err ~= "", error('endMismatch:segment', '%s', sg.err); end
            ref = sg.ref;  L = numel(ref);  Lf(m) = L / M;  cJ(m) = sg.info.cJump;  cS(m) = sg.info.cSlip;
            aRef(m) = alphaPWLW_v001(ref, FS, band);
            if isExt(m), tgt = r;  sR = std(r);  fr = frFull; else, tgt = ref;  sR = std(ref);  fr = filt_local(ref, FS); end
            for s = 1:nS
                ss = dIdx * 1e9 + tIdx * 1e4 + j * 1e3 + m * 100 + s;
                rs = RandStream("threefry4x64_20", "Seed", master);  rs.Substream = ss;  RandStream.setGlobalStream(rs);
                MF = [MF; table(string(tr.trialID), tIdx, ax(j), names(m), s, master, ss, 'VariableNames', {'trialID','trial','axis','kind','surrogate','master','substream'})]; %#ok<AGROW>
                nA(s, m) = 1;
                for dI = 1:maxDraw
                    switch names(m)
                        case {"e2e75", "e2e90", "e2eFair"}, u = phaseRand_local(ref);
                        case "extWrapFair", u0 = phaseRand_local(ref);  u = u0(mod(0:M - 1, L) + 1);
                        case {"extSpecFair", "extSpec75"}, u = specExtend_local(ref, M);
                    end
                    u = u - mean(u);  u = u / std(u) * sR;  dA = NaN;
                    try, dA = alphaPWLW_v001(u, FS, band) - aRef(m); catch, end      % a failed estimate is a rejected draw
                    fu = filt_local(u, FS);  q = [dA, fu ./ fr, fu(3:4) ./ frFull(3:4)];
                    if dI == 1, X(s, m, 1, :) = q; end
                    if abs(dA) <= tol, X(s, m, 2, :) = q;  nA(s, m) = 0;  break, end
                end
                nD(s, m) = dI;
            end
        catch ME
            st(m) = "FAIL " + names(m) + ": " + ME.message + " (in " + string(ME.stack(1).name) + ")";
        end
    end
    for m = 1:nM
        for w = 1:2
            md = @(k) median(X(:, m, w, k), 'omitnan');  kind = names(m) + extractBefore("Sel", 1 + 3 * (w - 1));
            nDw = 1;  nAw = 0;  if w == 2, nDw = median(nD(:, m));  nAw = mean(nA(:, m)); end
            rows{j, m, w} = table(string(tr.trialID), ax(j), kind, M, Lf(m), cJ(m), cS(m), f0, aR, aRef(m), md(1), ...
                md(2), md(3), md(4), md(5), md(6), md(7), nDw, nAw, nnz(isnan(X(:, m, w, 1))), st(m), ...
                'VariableNames', {'trialID','axis','kind','M','Lfrac','cJump','cSlip','f0','alphaResid','alphaRef','dAlphaMed', ...
                'velRatioSG','accRatioSG','velRatioBW','accRatioBW','velRatioBWFull','accRatioBWFull','nDrawMed','shareNoAccept','nSurrFailed','status'});
        end
    end
end
R = vertcat(rows{:});
end

function u = phaseRand_local(ref)
% |FFT(ref)| on uniformly random phases (Hermitian, so u is real); the mean is removed by the caller
L = numel(ref);  A = abs(fft(ref - mean(ref)));  h = floor((L - 1) / 2);
Z = zeros(L, 1);  Z(2:h + 1) = A(2:h + 1) .* exp(2i * pi * rand(h, 1));  Z(L:-1:L - h + 1) = conj(Z(2:h + 1));
if mod(L, 2) == 0, Z(L/2 + 1) = A(L/2 + 1) * sign(randn); end
u = real(ifft(Z));
end

function u = specExtend_local(ref, M)
% the segment's periodogram, interpolated linearly in normalised frequency onto the M-point grid, random phases
L = numel(ref);  P = abs(fft(ref - mean(ref))).^2 / L;  kL = (0:floor(L/2))';  fL = kL / L;
kM = (0:floor(M/2))';  fM = kM / M;  PM = interp1(fL, P(kL + 1), fM, 'linear', P(floor(L/2) + 1));  A = sqrt(PM * M);
h = floor((M - 1) / 2);  Z = zeros(M, 1);
Z(2:h + 1) = A(2:h + 1) .* exp(2i * pi * rand(h, 1));  Z(M:-1:M - h + 1) = conj(Z(2:h + 1));
if mod(M, 2) == 0, Z(M/2 + 1) = A(M/2 + 1) * sign(randn); end
u = real(ifft(Z));
end

function q = filt_local(u, FS)
% SD of the filtered velocity and acceleration, SG then BWFD, with the production call pattern and edge clip (v002)
E = 20;  q = NaN(1, 4);
d6 = differentiateKinematicsEBR(u, u, 6, [4 17], FS);     q(1) = std(d6(E:end-E, 2));  q(2) = std(d6(E:end-E, 3));
d2 = differentiateKinematicsEBR(u, u, 2, [2 10 1], FS);   q(3) = std(d2(E:end-E, 2));  q(4) = std(d2(E:end-E, 3));
end

% ---- copied verbatim from runLoopClosureFftnoise_v020.m (do not edit here; edit both or neither) ----
function f0 = estimateF0_local(x, y, fs)
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

function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);   % the production windows (identical to: s = 1; while s + win - 1 <= N, s = s + hop)
    if win > N, starts = []; end
    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window; v020: only when Template = 'endAnchored'
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
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    fitX = x(cs:ce) - resX(cs:ce);
    fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
    lastEnd = 0;  if ~isempty(starts), lastEnd = starts(end) + win - 1; end   % v007: drop the uncovered (zero) tail
    tailN = max(0, ce - lastEnd);  resX = resX(1:end - tailN);  resY = resY(1:end - tailN);
end
