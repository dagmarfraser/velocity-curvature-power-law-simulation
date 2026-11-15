function T = summariseEndMismatchSeeds_v006()
%SUMMARISEENDMISMATCHSEEDS_V006 Three-master-seed check of endMismatchSensitivity_v006 (README_SeedPractice rule 7).
%   Masters: 20261114 (results/endMismatchSensitivity/), 1729 and 314159 (subfolders master1729, master314159). Same machine
%   (M5). For each arm x axis x kind (e2e75Sel, e2e90Sel, e2eFairSel, extSpecFairSel) and statistic (BWFD velocity and
%   acceleration ratio, median over trial-axes), the median per master, its bootstrap SE (1000 resamples of trials, fixed
%   stream), and PASS if every pair of masters agrees within 2 combined SEs. Writes endMismatch_seedCheck_v006.csv.
here = fileparts(mfilename('fullpath'));  d = fullfile(here, '..', 'results', 'endMismatchSensitivity');
masters = [20261114 1729 314159];  dirs = [string(d), string(fullfile(d, "master1729")), string(fullfile(d, "master314159"))];
arms = ["FraserN48", "HickmanPLAC", "HickmanHALO", "CookCTRL", "CookASD"];  kinds = ["e2e75Sel", "e2e90Sel", "e2eFairSel", "extSpecFairSel"];
stats = ["velRatioBW", "accRatioBW"];  T = table();  bs = RandStream("threefry4x64_20", "Seed", 1);
for a = arms
    Rm = cell(1, 3);
    for k = 1:3, Rm{k} = readtable(fullfile(dirs(k), "endMismatch_" + a + "_v006.csv"), 'TextType', 'string'); end
    for ax = ["maj", "min"]
        for kd = kinds
            for s = stats
                med = NaN(1, 3);  se = NaN(1, 3);
                for k = 1:3
                    v = Rm{k}.(s)(Rm{k}.kind == kd & Rm{k}.axis == ax & Rm{k}.status == "ok");  v = v(isfinite(v));
                    med(k) = median(v);  bm = NaN(1000, 1);
                    for b = 1:1000, bm(b) = median(v(randi(bs, numel(v), numel(v), 1))); end
                    se(k) = std(bm);
                end
                ok = true;
                for p = 1:2, for q = p + 1:3, ok = ok && abs(med(p) - med(q)) <= 2 * hypot(se(p), se(q)); end, end
                T = [T; table(a, ax, kd, s, med(1), med(2), med(3), max(med) - min(med), max(se), ok, ...
                    'VariableNames', {'arm', 'axis', 'kind', 'stat', 'med20261114', 'med1729', 'med314159', 'range', 'maxSE', 'pass'})]; %#ok<AGROW>
            end
        end
    end
end
writetable(T, fullfile(d, "endMismatch_seedCheck_v006.csv"));
fprintf('%d of %d arm x axis x kind x statistic rows agree across three masters within 2 combined SEs\n', nnz(T.pass), height(T));
fprintf('largest range across masters: %.4f\n', max(T.range));
disp(T(~T.pass, :));
G = groupsummary(T, "kind", {"max", "mean"}, "range");  disp(G);
end
