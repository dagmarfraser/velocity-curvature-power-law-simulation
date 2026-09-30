%% checkAbstractBetaGenRange_v001.m
% Coherence item A4.14: the Abstract's "generator estimates of 0.28-0.33 in the five core
% datasets" is Table 6b's N_REPS = 20 basis, while the Abstract's held-out figure (84.5% of
% 20,514, #243) is on v015 (N_REPS = 200, the paper's basis since Session 114, D19). This
% script prints the per-dataset constellation median, median over trials of betaGenStarMed
% (the row-wise NaN-omitting median of the six pipelines; checkFraserConstellationMedian_v001),
% on both bases, so one basis can be stated.
%   V0 N_REPS = 20   Table 6b's corpora (v007; Fraser v008, as checkTable6bSourcePipeline_v001)
%                    regression anchor: must reproduce Table 6b's flat column
%                    (0.3292, 0.294, 0.284, 0.299, 0.298) to its printed precision.
%   V1 N_REPS = 200  the v015 corpora.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_{v007,v008,v015}.mat
% Writes: results/checkAbstractBetaGenRange_v001.mat
% USAGE:  from the project root: checkAbstractBetaGenRange_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
DS   = ["Fraser" "Cook_ASD" "Cook_CTRL" "Hickman_HALO" "Hickman_PLAC"];
V0F  = ["v008" "v007" "v007" "v007" "v007"];
PUB  = [0.3292 0.294 0.284 0.299 0.298];            % Table 6b, flat median
DP   = [4 3 3 3 3];                                 % printed precision
OUT_MAT = fullfile(ROOT, "results", "checkAbstractBetaGenRange_v001.mat");

%% Both bases
T = table();
for k = 1:numel(DS)
    m = NaN(1, 2);  n = NaN(1, 2);
    for b = 1:2
        ver = ifelse_local(b == 1, V0F(k), "v015");
        f = fullfile(ROOT, "src", "loopClosureResults_" + DS(k) + "_all_shaped_xu_" + ver + ".mat");
        if ~isfile(f), error("absBeta:input", "%s", "FAILED PATH: " + f); end
        R = load(f, "results").results;
        B = cell2mat(arrayfun(@(r) double(r.betaGenStar(:))', R(:), "UniformOutput", false));
        if size(B, 2) ~= 6, error("absBeta:width", "%s", "betaGenStar is not 6 wide in " + f); end
        bm = arrayfun(@(r) double(r.betaGenStarMed), R(:));
        dd = max(abs(bm - median(B, 2, "omitnan")), [], "omitnan");
        if dd > 1e-12, error("absBeta:med", "%s", sprintf("%s: betaGenStarMed is not the row-wise median (%.3g)", f, dd)); end
        m(b) = median(bm, "omitnan");  n(b) = nnz(isfinite(bm));
    end
    T = [T; table(DS(k), V0F(k), m(1), n(1), m(2), n(2), m(2) - m(1), ...
        'VariableNames', ["dataset" "basisN20" "medN20" "nN20" "medN200" "nN200" "delta"])]; %#ok<AGROW>
end

%% Regression anchor
if any(abs(arrayfun(@(x, d) round(x, d), T.medN20', DP) - PUB) > 1e-9)
    disp(T);  error("absBeta:anchor", "%s", "V0 does not reproduce Table 6b's flat column");
end
fprintf("REGRESSION ANCHOR passed: N_REPS = 20 reproduces Table 6b (%s).\n\n", strjoin(string(PUB), ", "));
disp(T)
fprintf("Range, N_REPS = 20: %.4f-%.4f (%.2f-%.2f)\n", min(T.medN20), max(T.medN20), min(T.medN20), max(T.medN20));
fprintf("Range, N_REPS = 200 (v015): %.4f-%.4f (%.2f-%.2f)\n", min(T.medN200), max(T.medN200), min(T.medN200), max(T.medN200));

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "T", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
