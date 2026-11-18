function T = shiftIntervals_v001(opts)
%SHIFTINTERVALS_V001 Paired bootstrap intervals for the shift in the dataset beta_gen* between the v028 primary arm and each sensitivity arm (2026-11-16).
%   Pre-commitment draft docs/V028_PreCommitment_DRAFT_20261116_v001.md section 3.3: two criteria are reported side by side, per dataset.
%     (a) the data-derived one: the shift's 95% interval (paired cluster bootstrap over subjects) includes zero or not;
%     (b) the fixed convention: |shift| below 0.005 (a stated convention, not derived), 0.005 to 0.03, and at least the MDC 0.03.
%   Statistic S(arm) = median over subjects of the subject's median betaGenStarMed (window basis, NaN omitted), the Table 6b primary statistic
%   without the inversion draw; the flat median over trials is reported beside it. Shift = S(arm) - S(primary), the same trials in both arms (checked).
%   Bootstrap: B resamples of subjects with replacement from one threefry stream (seed 20261116), both arms recomputed on the same resample.
%   Cross-check, not a rule: the seed-to-seed SD of S across the three masters (primary r200, 271828 r50, 314159 r50), times 2*sqrt(2) = 2.83 (the
%   simulation errors of two arms combine), against the bootstrap interval half-width. Arms: budget1 (endMatched, 1% budget), shapeXu (shapedXuSeg).
%   Opens nothing about one-third, cccMed, the 42-cell table, the tempo attenuation or the other outcomes. Writes results/v028/shift_intervals_v001.csv.
arguments
    opts.Dir (1,1) string = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main_N48/src"
    opts.B (1,1) double {mustBeInteger, mustBePositive} = 2000
end
here = fileparts(mfilename('fullpath'));  od = fullfile(here, '..', 'results', 'v028');  if ~isfolder(od), mkdir(od); end
ds = ["FraserN48", "Hickman_PLAC", "Hickman_HALO", "Cook_CTRL", "Cook_ASD", "Zarandi", "Dhieb"];
rs = RandStream("threefry4x64_20", "Seed", 20261116);  T = table();
for d = ds
    P = ld(opts.Dir, d, "endMatched_segdetect", "");  if isempty(P), continue, end
    subj = string({P.subjectID})';  bp = [P.betaGenStarMed]';  [~, ~, g] = unique(subj);  nS = max(g);
    Sp = subjMedian(bp, g, nS);  seeds = [Sp, NaN, NaN];  sd = ["271828", "314159"];
    for k = 1:2
        Q = ld(opts.Dir, d, "endMatched_segdetect", "_seed" + sd(k) + "_r50");
        if ~isempty(Q) && isequal(string({Q.trialID}), string({P.trialID})), seeds(k + 1) = subjMedian([Q.betaGenStarMed]', g, nS); end
    end
    for arm = ["budget1", "shapeXu"]
        if arm == "budget1", A = ld(opts.Dir, d, "endMatched_bud1", ""); else, A = ld(opts.Dir, d, "shapedXuSeg_segdetect", ""); end
        if isempty(A), continue, end
        if ~isequal(string({A.trialID}), string({P.trialID})), error('shift:ids', '%s %s: trial lists differ from the primary', d, arm); end
        ba = [A.betaGenStarMed]';  shift = subjMedian(ba, g, nS) - Sp;  flat = median(ba, 'omitnan') - median(bp, 'omitnan');
        mA = subjMediansOf(ba, g, nS);  mP = subjMediansOf(bp, g, nS);  bs = NaN(opts.B, 1);
        for b = 1:opts.B
            pick = randi(rs, nS, nS, 1);  bs(b) = median(mA(pick), 'omitnan') - median(mP(pick), 'omitnan');
        end
        lo = prctile(bs, 2.5);  hi = prctile(bs, 97.5);  mcSD = std(seeds, 'omitnan');
        T = [T; table(d, arm, nS, nnz(isfinite(bp) & isfinite(ba)), shift, lo, hi, (hi - lo) / 2, lo > 0 | hi < 0, abs(shift) < 0.005, abs(shift) >= 0.03, flat, mcSD, 2 * sqrt(2) * mcSD, ...
            'VariableNames', {'dataset', 'arm', 'nSubjects', 'nPairedTrials', 'shiftSubjectMedian', 'ci95Lo', 'ci95Hi', 'halfWidth', 'ciExcludesZero', 'below0005', 'atOrAboveMDC', 'shiftFlatMedian', 'seedSD', 'seedSDx283'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile(od, 'shift_intervals_v001.csv'));  disp(T);
end

function m = subjMediansOf(b, g, nS)
% each subject's median of b (NaN omitted); NaN for a subject with no finite value
m = NaN(nS, 1);
for i = 1:nS, v = b(g == i);  v = v(isfinite(v));  if ~isempty(v), m(i) = median(v); end, end
end

function s = subjMedian(b, g, nS)
s = median(subjMediansOf(b, g, nS), 'omitnan');
end

function R = ld(dir, d, tag, extra)
f = fullfile(dir, "loopClosureResults_" + d + "_all_" + tag + extra + "_v028.mat");
if isfile(f), R = load(f, 'results').results; else, R = []; fprintf('missing %s\n', f); end
end
