function T = harmonicShiftIntervals_v002(opts)
%HARMONICSHIFTINTERVALS_V002 Harmonic-count shift in the headline beta_gen*: K = nH arm minus K = 4 primary, per dataset, with paired subject-cluster bootstrap CI (2026-11-21).
%   Statistic S(arm) = median over subjects of each subject's median betaGenStarMed (NaN omitted), exactly as shiftIntervals_v001 (Table 6b primary statistic, no inversion draw).
%   Shift = S(K) - S(4). Reference K = 4 is the v028 primary file in the old tree src/ (loopClosureResults_<ds>_all_endMatched_segdetect_v028.mat).
%   Arms: nH1, nH2, nH3, nH6, nH10 (BRM/scratch_K, v028) and nH0 (BRM/outputs_K0, v033). nH1 comes first so that, with one shared stream, the Fraser and Hickman_PLAC nH1 rows
%   are the first two draw sets and reproduce the 2026-11-17 csv (validated below at the csv's 15 significant digits, raw difference also reported; a mismatch is an error, tables then not written).
%   Bootstrap: ONE threefry4x64_20 stream, seed 20261116, B = 2000 resamples of subjects with replacement (randi(rs,nS,nS,1)), both arms recomputed on the same resample.
%   Row order = arm (outer: 1,2,3,6,10,0) then dataset (inner: Fraser, PLAC, HALO, CTRL, ASD, Zarandi, Dhieb); every row consumes B draws even if its statistic is undefined.
%   Fail loud: missing file, md5 mismatch with the recorded lists, differing trialIDs/subjectIDs between arms, or a failed validation all raise errors; nothing is caught or dropped.
%   A statistic is "undefined" (status column, shift NaN) only when either arm has no finite subject median; nSubFiniteArm/nSubFiniteRef are reported for every row (K = 0 is mostly NaN).
%   Reads only the named local paths (no RDS). Writes BRM/outputs_harmonic/harmonicShiftIntervals_v002.{csv,mat,log}.
arguments
    opts.Root (1,1) string = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee"
    opts.B (1,1) double {mustBeInteger, mustBePositive} = 2000
    opts.Seed (1,1) double = 20261116
    opts.MDC (1,1) double = 0.03
end
oldTree = fullfile(opts.Root, "PowerLawSimulationPreReg");  brm = fullfile(opts.Root, "PowerLaw_Repro_2026", "BRM");
refDir = fullfile(oldTree, "src");  kDir = fullfile(brm, "scratch_K");  k0Dir = fullfile(brm, "outputs_K0");
od = fullfile(brm, "outputs_harmonic");  if ~isfolder(od), error('hsi:outdir', 'missing %s', od); end
oldCsv = fullfile(oldTree, "results", "v028", "harmonicShiftIntervals_20261117.csv");
if ~isfile(oldCsv), error('hsi:oldcsv', 'missing %s', oldCsv); end
logf = fullfile(od, "harmonicShiftIntervals_v002.log");  fid = fopen(logf, 'w');  cleanup = onCleanup(@() fclose(fid));
lg = @(varargin) fprintf(fid, varargin{:});
ds = ["FraserN48", "Hickman_PLAC", "Hickman_HALO", "Cook_CTRL", "Cook_ASD", "Zarandi", "Dhieb"];
dsLab = ["Fraser", "Hickman_PLAC", "Hickman_HALO", "Cook_CTRL", "Cook_ASD", "Zarandi", "Dhieb"];
Ks = [1 2 3 6 10 0];
% recorded md5 list for the K >= 1 scratch copies (written when they were copied from RDS)
rec = readlines(fullfile(kDir, "md5_local_K1plus.txt"));  recMap = containers.Map();
for r = rec(strlength(strtrim(rec)) > 0)'
    t = split(strtrim(r));  recMap(char(t(end))) = char(t(1));
end
rs = RandStream("threefry4x64_20", "Seed", opts.Seed);
T = table();  md5cache = containers.Map();  refCache = struct();
for K = Ks
    for i = 1:numel(ds)
        d = ds(i);
        % ---- reference K = 4 (loaded once per dataset)
        if ~isfield(refCache, d)
            fR = fullfile(refDir, "loopClosureResults_" + d + "_all_endMatched_segdetect_v028.mat");
            R = ldres(fR);  [~, ~, g] = unique(string({R.subjectID})');
            refCache.(d) = struct('f', fR, 'md5', md5of(fR), 'trial', string({R.trialID})', 'subj', string({R.subjectID})', 'g', g, 'b', [R.betaGenStarMed]');
        end
        ref = refCache.(d);  nS = max(ref.g);
        % ---- arm
        if K == 0
            fA = fullfile(k0Dir, "loopClosureResults_" + d + "_all_endMatched_segdetect_nH0_v033.mat");  ver = "v033";
        else
            fA = fullfile(kDir, "loopClosureResults_" + d + "_all_endMatched_segdetect_nH" + K + "_v028.mat");  ver = "v028";
            [~, nm, ex] = fileparts(fA);  key = char(nm + ex);
            if ~isKey(recMap, key), error('hsi:md5list', '%s not in md5_local_K1plus.txt', key); end
        end
        mA5 = md5of(fA);
        if K ~= 0 && ~strcmp(mA5, recMap(key)), error('hsi:md5', '%s md5 %s differs from recorded %s', fA, mA5, recMap(key)); end
        A = ldres(fA);
        if ~isequal(string({A.trialID})', ref.trial), error('hsi:ids', '%s K=%d: trialIDs differ from the K=4 primary', d, K); end
        if ~isequal(string({A.subjectID})', ref.subj), error('hsi:subj', '%s K=%d: subjectIDs differ from the K=4 primary', d, K); end
        ba = [A.betaGenStarMed]';  bp = ref.b;
        mA = subjMediansOf(ba, ref.g, nS);  mP = subjMediansOf(bp, ref.g, nS);
        SA = median(mA, 'omitnan');  SP = median(mP, 'omitnan');  shift = SA - SP;
        bs = NaN(opts.B, 1);
        for b = 1:opts.B
            pick = randi(rs, nS, nS, 1);  bs(b) = median(mA(pick), 'omitnan') - median(mP(pick), 'omitnan');
        end
        if isnan(shift), status = "undefined: no finite subject median in an arm"; lo = NaN; hi = NaN;
        else, status = "ok"; lo = prctile(bs, 2.5);  hi = prctile(bs, 97.5); end
        row = table(dsLab(i), "nH" + K, shift, lo, hi, lo > 0 | hi < 0, abs(shift) >= opts.MDC, nS, nnz(isfinite(mA)), nnz(isfinite(mP)), nnz(isfinite(bp) & isfinite(ba)), SA, SP, ...
            ver, mA5, ref.md5, string(fA), string(ref.f), status, ...
            'VariableNames', {'dataset', 'arm', 'shiftVs4', 'ciLo', 'ciHi', 'ciExcludesZero', 'atOrAboveMDC', 'nSubjects', 'nSubFiniteArm', 'nSubFiniteRef', 'nPairedTrials', ...
            'S_arm', 'S_K4', 'armVersion', 'md5Arm', 'md5K4', 'fileArm', 'fileK4', 'status'});
        T = [T; row]; %#ok<AGROW>
        lg('%s %s shift=%.15g ci=[%.15g, %.15g] excl0=%d nSubFin=%d/%d %s\n', dsLab(i), "nH" + K, shift, lo, hi, lo > 0 | hi < 0, nnz(isfinite(mA)), nnz(isfinite(mP)), status);
    end
end
% ---- validation against the 2026-11-17 csv: Fraser and Hickman_PLAC nH1 must match exactly (==)
O = readtable(oldCsv, 'TextType', 'string');  V = table();  bad = false;
for dd = ["Fraser", "Hickman_PLAC"]
    o = O(O.dataset == dd & O.arm == "nH1", :);  n = T(T.dataset == dd & T.arm == "nH1", :);
    if height(o) ~= 1 || height(n) ~= 1, error('hsi:val', 'validation row missing for %s nH1', dd); end
    dif = [n.shiftVs4 - o.shiftVs4, n.ciLo - o.ciLo, n.ciHi - o.ciHi];  ex = double(n.ciExcludesZero) == o.ciExcludesZero;
    % the old csv holds 15 significant digits (writetable default), so exact = identical when the new values are formatted %.15g; raw max abs diff is reported alongside
    txtMatch = isequal(compose("%.15g", [n.shiftVs4, n.ciLo, n.ciHi]), compose("%.15g", [o.shiftVs4, o.ciLo, o.ciHi]));
    V = [V; table(dd, o.shiftVs4, n.shiftVs4, o.ciLo, n.ciLo, o.ciHi, n.ciHi, max(abs(dif)), txtMatch && ex, 'VariableNames', {'dataset', 'oldShift', 'newShift', 'oldLo', 'newLo', 'oldHi', 'newHi', 'maxAbsDiff', 'exactMatch'})]; %#ok<AGROW>
    bad = bad || ~(txtMatch && ex);
end
disp(V);  lg('VALIDATION\n');  lg('%s\n', evalc('disp(V)'));
if bad, error('hsi:validation', 'nH1 rows do not reproduce harmonicShiftIntervals_20261117.csv exactly (max raw abs diff %g); outputs NOT written. See %s', max(V.maxAbsDiff), logf); end
% ---- outputs
csvf = fullfile(od, "harmonicShiftIntervals_v002.csv");  matf = fullfile(od, "harmonicShiftIntervals_v002.mat");
writetable(T, csvf);
meta = struct('script', "harmonicShiftIntervals_v002.m", 'seed', opts.Seed, 'B', opts.B, 'MDC', opts.MDC, 'created', string(datetime('now')), 'matlab', version);
save(matf, 'T', 'V', 'meta');
m5 = ["csv  " + md5of(csvf); "mat  " + md5of(matf)];  lg('%s\n', strjoin(m5, newline));
disp(T(:, {'dataset', 'arm', 'shiftVs4', 'ciLo', 'ciHi', 'ciExcludesZero', 'atOrAboveMDC', 'nSubFiniteArm', 'status'}));
end

function m = subjMediansOf(b, g, nS)
% each subject's median of b (NaN omitted); NaN for a subject with no finite value
m = NaN(nS, 1);
for i = 1:nS, v = b(g == i);  v = v(isfinite(v));  if ~isempty(v), m(i) = median(v); end, end
end

function R = ldres(f)
if ~isfile(f), error('hsi:missing', 'missing file %s', f); end
S = load(f, 'results');  R = S.results;
end

function h = md5of(f)
[st, out] = system("/sbin/md5 -q '" + f + "'");
if st ~= 0, error('hsi:md5cmd', 'md5 failed for %s: %s', f, out); end
h = string(strtrim(out));
end
