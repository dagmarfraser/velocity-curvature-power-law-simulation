function T = compareFraserN45N48Descriptives_v001()
%COMPAREFRASERN45N48DESCRIPTIVES_V001  Descriptive Fraser N=45 (v015 corpus) against N=48 (importDB_fraser_v002) table (2026-11-03).
%   Only quantities that need no forward map: participants, trials, sessions, trials per person, tempo f0,
%   beta_obs per pipeline (dataset median of trial values and median of participant medians), noise alpha and sigma.
%   sigma is given in mm at BOTH scales for N=45: 1/10.41793 (as v015 and v008 used) and 1/5.209 (UIKit points per mm,
%   the Step 1 audit's correction; it applies to the old corpus too, which is in the same points). N=48 uses 1/5.209.
%   Sources: N=45 src/loopClosureResults_Fraser_all_shaped_xu_v015.mat and src/noiseCharacterisation_fraser.mat;
%   N=48 results/auditFraserN48ObservedFits_v001.mat (beta_obs, same code path as the runner),
%   src/noiseCharacterisation_fraserN48.mat, importDB_fraser_v002 (f0 computed with the runner's estimateF0_local, copied verbatim).
%   Self-check: N=45 f0 recomputed here must equal v015's stored f0 on trials present in both imports and unchanged by v002.
%   Writes results/compareFraserN45N48Descriptives_v001.mat.   Fraser, D.S. (2026)  v001

srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));
resDir = fullfile(fileparts(srcDir), 'results');
PT_PER_MM = 5.209;  OLD_PX_PER_MM = 10.41793;  FS = 240;
labels = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];

%% N=45 (v015)
S = load(fullfile(srcDir, 'loopClosureResults_Fraser_all_shaped_xu_v015.mat'), 'results');  R = S.results;
b45   = cell2mat(arrayfun(@(r) r.betaObs(:)', R, 'UniformOutput', false)');
f045  = [R.f0]';
id45  = string({R.trialID})';  sub45 = string({R.subjectID})';
n45   = load(fullfile(srcDir, 'noiseCharacterisation_fraser.mat'), 'bioResults').bioResults;

%% N=48 (importDB_fraser_v002 + Step 2b audit)
A = load(fullfile(resDir, 'auditFraserN48ObservedFits_v001.mat'), 'A').A;
trials = importDB_fraser_v002(Verbose=false);
id48 = string({trials.trialID})';
if ~isequal(id48, A.trialID), error('compareN45N48:order', '%s', 'Audit trial order differs from importDB_fraser_v002: rerun the audit.'); end
b48 = A.betaObs;  sub48 = string({trials.subjectID})';
f048 = arrayfun(@(t) estimateF0_local(double(t.x), double(t.y), FS), trials);
n48 = load(fullfile(srcDir, 'noiseCharacterisation_fraserN48.mat'), 'bioResults').bioResults;
if ~isequal(string(n48.trialID), id48), error('compareN45N48:noiseOrder', '%s', 'FraserN48 noise file rows do not match the import.'); end

%% Self-check: f0 on unchanged common trials
[inBoth, loc] = ismember(id48, id45);
same = inBoth & ~A.changed;
d = abs(f048(same) - f045(loc(same)));
if max(d) > 0, error('compareN45N48:f0', 'f0 recomputed here differs from v015 on unchanged trials (max %.3g Hz).', max(d)); end
fprintf('Self-check: f0 identical to v015 on %d unchanged common trials.\n', sum(same));

%% Table
q = @(v) sprintf('%.3f [%.3f, %.3f]', median(v, 'omitnan'), prctile(v, 25), prctile(v, 75));
pm = @(v, s) median(splitapply(@(x) median(x, 'omitnan'), v, findgroups(s)), 'omitnan');   % median of participant medians
sess = @(ids) numel(unique(extractBefore(ids, "_T")));
tpp  = @(s) median(groupcounts(s));
rows = {
    'Participants',                    numel(unique(sub45)),            numel(unique(sub48));
    'Trials',                          numel(id45),                     numel(id48);
    'Sessions (participant x visit x session)', sess(id45),             sess(id48);
    'Trials per participant (median)', tpp(sub45),                      tpp(sub48);
    'Tempo f0, Hz: median [IQR]',      q(f045),                         q(f048);
    'Noise alpha (ira_alphaMean): median [IQR]', q(n45.ira_alphaMean),  q(n48.ira_alphaMean);
    'sigma, mm at 10.418/mm (as v015): median [IQR]', q(n45.sigmaMean / OLD_PX_PER_MM), '-';
    'sigma, mm at 5.209/mm (corrected): median [IQR]', q(n45.sigmaMean / PT_PER_MM),    q(n48.sigmaMean / PT_PER_MM)};
for p = 1:6
    rows(end + 1, :) = {sprintf('beta_obs %s: trial median [IQR]', labels(p)), q(b45(:, p)), q(b48(:, p))}; %#ok<AGROW>
    rows(end + 1, :) = {sprintf('beta_obs %s: median of participant medians', labels(p)), pm(b45(:, p), sub45), pm(b48(:, p), sub48)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', ["quantity" "N45_v015" "N48_v002"]);
disp(T)
out = fullfile(resDir, 'compareFraserN45N48Descriptives_v001.mat');
save(out, 'T', 'f048', 'id48', 'PT_PER_MM', 'OLD_PX_PER_MM');
fprintf('Saved %s\n', out);
end

% --------------------------------------------------------------------------
function f0 = estimateF0_local(x, y, fs)   % verbatim from runLoopClosureFftnoise_v015/v019/v020
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
