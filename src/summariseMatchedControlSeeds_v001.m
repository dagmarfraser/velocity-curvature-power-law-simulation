function S = summariseMatchedControlSeeds_v001(opts)
% SUMMARISEMATCHEDCONTROLSEEDS_V001  Combine positiveControlMatchedByDataset_v006 chunk files and summarise the per-dataset matched controls over
%   replicates and master seeds (2026-11-07; seed practice of the fractal-noise validation: summaries over many replicates, never one draw,
%   and agreement across master seeds within Monte Carlo error as the pass condition).
%   Per dataset, primary master (default 1): median and 5-95% range of medA, medR and lambdaStar_own; share of replicates meeting the pre-set
%   windows (|medA| <= 0.30, 0.70 <= medR <= 1.30); share with lambdaStar_own > 0.5. Monte Carlo SE of each median by bootstrap over replicates
%   (B = 2000, rng 7). STABILITY: for every other master and every dataset, |median(other) - median(primary)| <= 2 * sqrt(SE_p^2 + SE_o^2) for
%   medA, medR and lambdaStar_own; reported per cell and overall (PASS only if every cell passes).
%   Duplicate (master, dataset, rep) rows across files are an error (Fail Loud).
%   S = summariseMatchedControlSeeds_v001(Dir=<results folder>, Primary=1). Saves results/summariseMatchedControlSeeds_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Primary (1,1) double = 1
    opts.B (1,1) double = 2000
end
srcDir = fileparts(mfilename("fullpath"));  if opts.Dir == "", opts.Dir = string(fullfile(srcDir, "..", "results")); end
f = dir(fullfile(opts.Dir, "positiveControlMatchedByDataset_v006_m*.mat"));  f = f(~contains({f.name}, "_pilot"));
if isempty(f), error("mcs:none", "%s", "No v006 chunk files in " + opts.Dir); end
rows = table();
for k = 1:numel(f), L = load(fullfile(f(k).folder, f(k).name), "rows");  rows = [rows; L.rows]; end %#ok<AGROW>
key = rows.master * 1e6 + double(findgroups(rows.dataset)) * 1e4 + rows.rep;
if numel(unique(key)) ~= height(rows), error("mcs:dup", "%s", "Duplicate (master, dataset, rep) rows across chunk files"); end
rng(7, "twister");  V = ["medA" "medR" "lambdaStar_own"];  ds = unique(rows.dataset, "stable");  ms = unique(rows.master)';
sumT = table();  stab = table();
for m = ms
    for d = ds'
        x = rows(rows.master == m & rows.dataset == d, :);  if isempty(x), continue, end
        r = table(m, d, height(x), mean(x.calOK), mean(x.outside), 'VariableNames', ["master","dataset","nReps","shareCalOK","shareOutside"]);
        for v = V
            y = x.(v);  bs = bootstrp(opts.B, @(z) median(z, "omitnan"), y);
            r.(v + "_med") = median(y, "omitnan");  r.(v + "_p5") = prctile(y, 5);  r.(v + "_p95") = prctile(y, 95);  r.(v + "_mcse") = std(bs);
        end
        sumT = [sumT; r]; %#ok<AGROW>
    end
end
P = sumT(sumT.master == opts.Primary, :);
for m = ms(ms ~= opts.Primary)
    for d = ds'
        a = P(P.dataset == d, :);  b = sumT(sumT.master == m & sumT.dataset == d, :);  if isempty(a) || isempty(b), continue, end
        for v = V
            dlt = b.(v + "_med") - a.(v + "_med");  tol = 2 * hypot(a.(v + "_mcse"), b.(v + "_mcse"));
            stab = [stab; table(m, d, v, dlt, tol, abs(dlt) <= tol, 'VariableNames', ["master","dataset","stat","delta","tol2se","pass"])]; %#ok<AGROW>
        end
    end
end
S = struct("rows", rows, "summary", sumT, "stability", stab, "stablePASS", ~isempty(stab) && all(stab.pass), "files", string({f.name})');
disp(sumT(:, ["master","dataset","nReps","shareCalOK","shareOutside","medA_med","medR_med","lambdaStar_own_med","lambdaStar_own_p5","lambdaStar_own_p95"]));
if ~isempty(stab), disp(stab); fprintf("STABILITY across masters: %s (%d of %d cells within 2 MC SE)\n", string(iff(S.stablePASS, "PASS", "FAIL")), nnz(stab.pass), height(stab)); end
save(fullfile(opts.Dir, "summariseMatchedControlSeeds_v001.mat"), "S");  fprintf("Saved %s\n", fullfile(opts.Dir, "summariseMatchedControlSeeds_v001.mat"));
end

function y = iff(c, a, b)
if c, y = a; else, y = b; end
end
