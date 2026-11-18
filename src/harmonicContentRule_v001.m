function out = harmonicContentRule_v001(opts)
%HARMONICCONTENTRULE_V001 The content rule of Amendment 01 (fallback of section 3 item 4): the template's harmonic count per dataset from position data only (2026-11-17).
%   Frozen definition (docs/V028_PreCommitment_AMENDMENT01_20261117_FROZEN_v001.md, section 3 item 4): per dataset, the highest harmonic k whose pooled amplitude-to-inter-harmonic-noise MEDIAN ratio is at
%   least 2, with all lower harmonics included (the counting stops at the first harmonic below 2); method of harmonicCountDiagnostic_v002 (harmAmp_local: one least-squares fit per record and axis,
%   constant and linear trend plus a cosine and sine at each half-multiple of f0; amp(k) at k f0; ampNull(k) the geometric mean of the amplitudes at (k - 0.5) f0 and (k + 0.5) f0; ratio = amp / ampNull,
%   mean of the two axes; f0 from the runner's FFT-peak estimate). The count is at least 1.
%   Trials: the dataset's trials in the v028 primary result file (trial IDs only are read from it; no beta_gen* or other result field is used). Position data only.
%   out = harmonicContentRule_v001();   Writes results/harmonicContentRule_v001<Tag>.mat
arguments
    opts.Datasets (1,:) string = ["FraserN48" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"]
    opts.Cut (1,1) double = 2
    opts.Kmax (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.Tag (1,1) string = ""
    opts.UseParfor (1,1) logical = true
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
ws = warning('off', 'all');  wCleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
K = opts.Kmax;  S = table();  detail = struct();
for d = opts.Datasets
    switch d
        case "FraserN48",      tr = importDB_fraser_v002(Verbose=false);
        case "Cook_CTRL",      tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook_ASD",       tr = importDB_cook_v002(Group="ASD",  Tasks=7, Verbose=false);
        case "Hickman_PLAC",   tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman_HALO",   tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        case "Zarandi",        tr = importDB_zarandi_v001(Verbose=false);
        case "Dhieb",          tr = importDB_dhieb_v001(Verbose=false);
    end
    R = loadV028_v001(strrep(d, "FraserN48", "Fraser"));  keep = string({R.trialID});   % trial IDs only
    tr = tr(ismember(string({tr.trialID}), keep));  n = numel(tr);
    X = {tr.x};  Y = {tr.y};  FS = [tr.fs];  f0 = NaN(n, 1);  amp = NaN(n, K);  ampNull = NaN(n, K);
    if opts.UseParfor
        try, ensurePool_v001(); catch, end   %#ok<TRYNC>
        parfor i = 1:n, [f0(i), amp(i, :), ampNull(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), K); end
    else
        for i = 1:n, [f0(i), amp(i, :), ampNull(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), K); end
    end
    snr = amp ./ ampNull;  med = median(snr, 1, 'omitnan');  share = mean(snr > opts.Cut, 1, 'omitnan');
    below = find(~(med >= opts.Cut), 1);  if isempty(below), Kc = K; else, Kc = max(1, below - 1); end
    nUsable = nnz(all(isfinite(snr(:, 1:min(3, K))), 2));
    fprintf('%-13s n %5d (usable %5d)  median f0 %.2f Hz  medians k=1..6: %s  -> content count K = %d\n', d, n, nUsable, median(f0, 'omitnan'), mat2str(round(med(1:6), 2)), Kc);
    S = [S; table(d, n, nUsable, median(f0, 'omitnan'), Kc, med(1), med(2), med(3), med(4), share(1), share(2), share(3), ...
        'VariableNames', {'dataset', 'nTrials', 'nUsable', 'f0Median', 'Kcontent', 'snrMed1', 'snrMed2', 'snrMed3', 'snrMed4', 'shareAbove1', 'shareAbove2', 'shareAbove3'})]; %#ok<AGROW>
    detail.(matlab.lang.makeValidName(d)) = struct('snrMedian', med, 'shareAboveCut', share, 'f0', f0);
end
disp(S);
out = struct('summary', S, 'detail', detail, 'opts', opts, 'generated', string(datetime('now')));
f = fullfile(resDir, "harmonicContentRule_v001" + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
end

function [f0, amp, ampNull] = oneTrial_local(x, y, fs, K)
amp = NaN(1, K);  ampNull = NaN(1, K);  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
x = detrend(x(:), 'linear');  y = detrend(y(:), 'linear');  t = (0:numel(x)-1)' / fs;
fr = (0.5:0.5:K + 0.5) * f0;  fr = fr(fr < fs / 2 - 0.5);  nf = numel(fr);  A = NaN(2, nf);
D0 = [ones(numel(t), 1), t];
for j = 1:nf
    D = [D0, cos(2*pi*fr(j)*t), sin(2*pi*fr(j)*t)];
    bx = D \ x;  by = D \ y;  A(1, j) = hypot(bx(3), bx(4));  A(2, j) = hypot(by(3), by(4));
end
Am = mean(A, 1);
for k = 1:K
    ik = 2 * k;  if ik + 1 > nf, break, end
    amp(k) = Am(ik);  ampNull(k) = sqrt(Am(ik - 1) * Am(ik + 1));
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
