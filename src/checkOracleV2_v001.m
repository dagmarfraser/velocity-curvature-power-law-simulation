%% checkOracleV2_v001.m
% Gate V2 of SPEC_TempoFoldSweep_v001: was Finding #209's oracle right for its inputs, and
% does the tempo account (Findings #224/#225) explain why it found no zero-noise fold?
%   1. Reproduce #209's sigma = 0 residuals from the oracle's own saved output
%      (measureRegressionSaturation_v005.mat, written by saturationSweepEngine_v001:
%      SG case 4 [4 17], edge clip 20, OLS-first seeding; worst case reported 0.0197).
%   2. Re-run the oracle geometries through tempoFoldEngine_v001 at their own f0
%      (SG case 6, edge clip 20, runner seeding). At fs 100 (Maoz, Zarandi) SG cases 4 and 6
%      coincide, so these must match; at 133/240 Hz the difference is the SG window alone.
%   3. Same geometries at fixed f0 of 1, 2.5 and 4 Hz: residual vs tempo.
% sigma = 0 throughout (deterministic; one replicate). Runs locally (serial engine).
% Reads:  src/measureRegressionSaturation_v005.mat
% Writes: results/oracleV2_v001.mat

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
ORACLE   = fullfile(ROOT, "src", "measureRegressionSaturation_v005.mat");
OUT_MAT  = fullfile(ROOT, "results", "oracleV2_v001.mat");
% Geometries exactly as measureRegressionSaturation_v005.m (name, a, b, f0, FS, nCycles)
G = table(["Maoz"; "Cook CTRL"; "Hickman"; "Pilot"; "Zarandi"], [50.0; 57.5; 57.0; 30.8; 64.0], ...
          [25.0; 27.1; 26.5; 13.8; 15.6], [1.000; 0.369; 0.828; 0.732; 0.647], [100; 133; 133; 240; 100], ...
          10 * ones(5, 1), 'VariableNames', ["name" "a" "b" "f0" "FS" "nCycles"]);
F0_EXTRA = [1 2.5 4];
BETA     = linspace(0, 0.75, 10);           % measureRegressionSaturation_v005.m betaGenSweep
LABELS   = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];   % engine and oracle order
CLIP     = 20;

%% Setup
addpath(fullfile(ROOT, "src"));  addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
if ~isfile(ORACLE), error("oracleV2:input", "%s", "Missing oracle output: " + ORACLE); end
if ~isfolder(fileparts(OUT_MAT)), error("oracleV2:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
O = load(ORACLE, "meanBetaObs", "betaAnalytic", "sigmaSweep", "betaGenSweep", "pipelineNames", "geometries");
if ~isequal(string(O.pipelineNames), LABELS), error("oracleV2:labels", "%s", "Oracle pipeline order differs: " + strjoin(string(O.pipelineNames), ", ")); end
if max(abs(O.betaGenSweep(:)' - BETA)) > 1e-12, error("oracleV2:beta", "Oracle beta_gen sweep differs from CONFIG"); end
if ~isequal(string({O.geometries.name})', G.name), error("oracleV2:geom", "Oracle geometry order differs"); end
s0 = find(O.sigmaSweep == 0, 1);
if isempty(s0), error("oracleV2:sigma0", "Oracle output has no sigma = 0 column"); end

%% 1. The oracle's own sigma = 0 residuals (#209)
orc = squeeze(O.meanBetaObs(:, :, :, s0));                     % nGeom x 6 x nBeta
resO = orc - reshape(BETA, 1, 1, []);
fprintf("1. Oracle's own sigma = 0 residual |beta_rec - beta_gen| (max over beta_gen, finite values):\n");
t1 = array2table(round(max(abs(resO), [], 3, "omitnan"), 4), 'VariableNames', LABELS, 'RowNames', cellstr(G.name));
disp(t1)
[mx, ix] = max(abs(resO(:)), [], "omitnan");  [gi, pi_, bi] = ind2sub(size(resO), ix);
fprintf("   worst: %.4f (%s, %s, beta_gen %.3f)  [#209 reports 0.0197]\n", mx, G.name(gi), LABELS(pi_), BETA(bi));

%% 2-3. Engine at own f0 and at fixed higher tempos
J = table();
for g = 1:height(G)
    f0s = unique([G.f0(g) F0_EXTRA]);
    [f0, bg] = ndgrid(f0s, BETA);  n = numel(f0);
    J = [J; table(repmat("V2", n, 1), repmat("ellipse", n, 1), repmat(G.a(g)/G.b(g), n, 1), ...
        repmat(G.a(g), n, 1), repmat(G.b(g), n, 1), f0(:), NaN(n, 1), repmat(G.FS(g), n, 1), bg(:), ...
        zeros(n, 1), zeros(n, 1), ones(n, 1), repmat(CLIP, n, 1), repmat(G.nCycles(g), n, 1), ...
        repmat("runner", n, 1), zeros(n, 1), repmat(G.name(g), n, 1), 'VariableNames', ["block" "generator" ...
        "ab" "a" "b" "f0" "vgf" "FS" "betaGen" "sigma" "alpha" "nReps" "edgeClip" "nCycles" "seedMode" "jobID" "geom"])]; %#ok<AGROW>
end
J.jobID = (1:height(J))';
fprintf("\nEngine jobs: %d (sigma = 0)\n", height(J));
R = tempoFoldEngine_v001(J, struct("nWorkers", 0, "rngSeed", 1729));
fprintf("Pipeline failures: %s\n", strjoin(compose("%s %d", LABELS', sum(R.nFail, 1)'), ", "));

%% 2. Engine at the oracle's own f0 vs the oracle
t2 = table();
for g = 1:height(G)
    r = sortrows(R(R.geom == G.name(g) & R.f0 == G.f0(g), :), "betaGen");
    d = abs(r.betaMean - squeeze(orc(g, :, :))');                % nBeta x 6
    t2 = [t2; array2table(round(max(d, [], 1, "omitnan"), 5), 'VariableNames', LABELS, 'RowNames', cellstr(G.name(g) + " (fs " + G.FS(g) + ")"))]; %#ok<AGROW>
end
fprintf("\n2. Engine (SG case 6) vs oracle (SG case 4), own f0, max |diff| over beta_gen\n");
fprintf("   (fs 100 rows must be ~0: the SG cases coincide there)\n");  disp(t2)

%% 3. Residual vs tempo
t3 = table();
for g = 1:height(G)
    for f0 = unique([G.f0(g) F0_EXTRA])
        r = sortrows(R(R.geom == G.name(g) & R.f0 == f0, :), "betaGen");
        res = abs(r.betaMean - r.betaGen);
        [~, i23] = min(abs(r.betaGen - 2/3));
        t3 = [t3; table(G.name(g), f0, f0 == G.f0(g), max(res(:, 1)), max(res(:, 6)), max(res, [], "all", "omitnan"), ...
            r.betaMean(i23, 6), 'VariableNames', ["geometry" "f0" "ownTempo" "maxRes_BWFD_OLS" "maxRes_SG_IRLS" ...
            "maxRes_any" "SG_IRLS_at_2over3"])]; %#ok<AGROW>
    end
end
t3{:, 4:end} = round(t3{:, 4:end}, 4);
fprintf("\n3. sigma = 0 residual |beta_rec - beta_gen| by tempo (beta_gen 0-0.75; 2/3 read at the nearest node, %.4f)\n", ...
    BETA(find(abs(BETA - 2/3) == min(abs(BETA - 2/3)), 1)));
disp(t3)

save(OUT_MAT, "R", "J", "G", "t1", "t2", "t3", "BETA", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
