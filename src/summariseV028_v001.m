function [A, Sd] = summariseV028_v001(opts)
%SUMMARISEV028_V001 First summary of the v028 end-matched runs (2026-11-16): arms, failure reasons, window vs whole trial, seeds.
%   Reads loopClosureResults_<dataset>_all_<model>_<rule>[_seed<S>_r50]_v028.mat from the RDS N48 src folder (mounted) or
%   opts.Dir, and writes results/v028/summary_arms_v001.csv and summary_seeds_v001.csv.
%   Arms: A "primary" (endMatched, detectability rule), B "budget1" (endMatched, 1% budget), C "shapeXu" (shapedXuSeg, detectability).
%   Per dataset and arm: trials, trials that ran (emStatus "ok"), trials with a finite beta_gen* on the window basis (primary
%   basis) and the whole-trial basis, median and IQR of beta_gen* (betaGenStarMed, constellation median over the six
%   pipelines), median SG-IRLS beta_gen*, the paired per-trial window-minus-whole difference, and why trials did not run
%   (segment shorter than 280 / no segment / residual below the minimum / other).
%   Seeds: the primary arm at masters 1729 (N_REPS 200), 271828 and 314159 (N_REPS 50): per dataset the median beta_gen* per
%   master, their range, and the median |per-trial difference| of each seed rerun from the primary.
%   Descriptive only; no template or noise comparison is drawn here.
arguments
    opts.Dir (1,1) string = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main_N48/src"
end
here = fileparts(mfilename('fullpath'));  od = fullfile(here, '..', 'results', 'v028');  if ~isfolder(od), mkdir(od); end
ds = ["FraserN48", "Hickman_PLAC", "Hickman_HALO", "Cook_CTRL", "Cook_ASD", "Zarandi", "Dhieb"];
arms = struct('name', {"primary", "budget1", "shapeXu"}, 'tag', {"endMatched_segdetect", "endMatched_bud1", "shapedXuSeg_segdetect"});
A = table();  P = struct();
for d = ds
    for a = arms
        f = fullfile(opts.Dir, "loopClosureResults_" + d + "_all_" + a.tag + "_v028.mat");
        if ~isfile(f), fprintf('missing %s\n', f);  continue, end
        r = load(f, 'results').results;  P.(d + "_" + a.name) = r;
        st = [r.emStatus];  ran = st == "ok";  bw = [r.betaGenStarMed]';  bwh = [r.betaGenStarMedWhole]';
        sg = arrayfun(@(x) x.betaGenStar(6), r)';  sgh = arrayfun(@(x) x.betaGenStarWhole(6), r)';
        dd = bw - bwh;  dd = dd(isfinite(dd));
        rs = struct('short', nnz(startsWith(st, "FAIL: end-matched length")), 'none', nnz(contains(st, "no length of")), ...
            'resid', nnz(startsWith(st, "not reached")), 'other', 0);
        rs.other = nnz(~ran) - rs.short - rs.none - rs.resid;
        A = [A; table(d, a.name, numel(r), nnz(ran), nnz(isfinite(bw)), nnz(isfinite(bwh)), median(bw, 'omitnan'), prctile(bw, 25), prctile(bw, 75), ...
            median(bwh, 'omitnan'), median(sg, 'omitnan'), median(sgh, 'omitnan'), median(dd), prctile(abs(dd), 90), rs.short, rs.none, rs.resid, rs.other, ...
            median([r(ran).emLfrac], 'omitnan'), ...
            'VariableNames', {'dataset', 'arm', 'nTrials', 'nRan', 'nInvWindow', 'nInvWhole', 'medBetaGenWindow', 'q25', 'q75', 'medBetaGenWhole', 'medSGIRLSWindow', ...
            'medSGIRLSWhole', 'medWindowMinusWhole', 'p90AbsWindowMinusWhole', 'failShort', 'failNoSegment', 'failResidual', 'failOther', 'medSegmentShare'})]; %#ok<AGROW>
    end
end
writetable(A, fullfile(od, 'summary_arms_v001.csv'));  disp(A);
Sd = table();
for d = ds
    key = d + "_primary";  if ~isfield(P, key), continue, end
    r0 = P.(key);  ids0 = string({r0.trialID});  b0 = [r0.betaGenStarMed];  meds = median(b0, 'omitnan');  dif = NaN(1, 2);  n = NaN(1, 2);
    for k = 1:2
        s = ["271828", "314159"];
        f = fullfile(opts.Dir, "loopClosureResults_" + d + "_all_endMatched_segdetect_seed" + s(k) + "_r50_v028.mat");
        if ~isfile(f), continue, end
        r = load(f, 'results').results;  [~, i0, i1] = intersect(ids0, string({r.trialID}));
        x = [r(i1).betaGenStarMed] - b0(i0);  x = x(isfinite(x));  dif(k) = median(abs(x));  n(k) = numel(x);  meds(end + 1) = median([r.betaGenStarMed], 'omitnan'); %#ok<AGROW>
    end
    Sd = [Sd; table(d, meds(1), meds(min(2, end)), meds(min(3, end)), max(meds) - min(meds), dif(1), dif(2), n(1), n(2), ...
        'VariableNames', {'dataset', 'med1729r200', 'med271828r50', 'med314159r50', 'rangeOfMedians', 'medAbsDiff271828', 'medAbsDiff314159', 'nPaired271828', 'nPaired314159'})]; %#ok<AGROW>
end
writetable(Sd, fullfile(od, 'summary_seeds_v001.csv'));  disp(Sd);
end
