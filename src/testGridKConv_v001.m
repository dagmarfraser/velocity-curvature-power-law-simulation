%% testGridKConv_v001.m
% Verifies gridKConv_v001 (Fail Loud: any failure errors).
%   1. beta is required: no-argument, negative, NaN and complex calls all error.
%   2. K(beta) is strictly decreasing on the grid's beta nodes and K(0) equals the path length.
%   3. Agreement with the DB's realised tempo, f0 = VGF / K(beta), over every perNode row
%      (results/gridTempoAtFold_v004.mat): asserted for beta <= 0.40, reported above it.
%   4. The retired K_conv (perimeter / mean(kappa^-1/3), 175.2636) is 5.28% below K(1/3).
% TOLERANCES: MAX_ERR_LOW_BETA and MED_ERR_LOW_BETA are set from the observed values on
%   the 504 nodes at beta <= 0.40 (median 0.041%, max 0.553%; 2026-09-29), not derived
%   bounds. They guard against regression, not against a theory of the discretisation.
% Writes: results/testGridKConv_v001.mat
% USAGE:  from the project root: testGridKConv_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
V0_CANDS = [fullfile(ROOT, "results", "gridTempoAtFold_v004.mat")
            "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results/gridTempoAtFold_v004.mat"];
OUT_MAT  = fullfile(ROOT, "results", "testGridKConv_v001.mat");
BETA_MAX_ASSERT   = 0.40 + 1e-3;
MAX_ERR_LOW_BETA  = 0.006;        % observed 0.00553
MED_ERR_LOW_BETA  = 0.001;        % observed 0.00041
KCONV_RETIRED     = 175.2636;     % perimeter / mean(kappa^-1/3), retired
RETIRED_GAP_EXPECT = [-0.0535 -0.0520];   % KCONV_RETIRED / K(1/3) - 1, observed -0.0528
addpath(ROOT); addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));

%% 1. beta is required and validated
bad = {{}, {-0.1}, {NaN}, {Inf}, {1 + 2i}, {"a"}};
labels = ["no argument" "negative" "NaN" "Inf" "complex" "string"];
for i = 1:numel(bad)
    threw = false;
    try
        gridKConv_v001(bad{i}{:});
    catch
        threw = true;
    end
    if ~threw, error("testKConv:noError", "%s", "gridKConv_v001 accepted: " + labels(i)); end
end
fprintf("1. required/validated beta: %d bad calls all error\n", numel(bad));

%% 2. Shape of K(beta)
b = (0:21) / 30;                                    % the grid's 22 nodes, 0 to 0.70
K = gridKConv_v001(b);
if ~isequal(size(K), size(b)), error("testKConv:size", "%s", "K does not match beta's size"); end
if any(diff(K) >= 0), error("testKConv:monotone", "%s", "K(beta) is not strictly decreasing"); end
addpath(genpath(fullfile(ROOT, "src")));
G = gridShapeGeometry_v002(ROOT);
if abs(K(1) / G.perimPx - 1) > 1e-6
    error("testKConv:K0", "%s", sprintf("K(0) %.6f differs from the shape's perimeter %.6f", K(1), G.perimPx));
end
fprintf("2. K strictly decreasing over %d nodes; K(0) = perimeter = %.3f px\n", numel(b), K(1));

%% 3. Realised tempo from the DB
iV = find(isfile(V0_CANDS), 1);
if isempty(iV), error("testKConv:V0", "%s", "V0 output not found: " + strjoin(V0_CANDS, " | ")); end
V = load(V0_CANDS(iV), "out");  F = V.out.perNode;
F = F(F.generated_beta > 0 & isfinite(F.f0), :);
fa  = F.vgf_value ./ gridKConv_v001(F.generated_beta);
err = abs(F.f0 ./ fa - 1);
lo  = F.generated_beta <= BETA_MAX_ASSERT;
fprintf("3. realised vs VGF/K(beta): beta <= 0.40  n = %d  median %.3f%%  max %.3f%%\n", nnz(lo), 100 * median(err(lo)), 100 * max(err(lo)));
fprintf("                            beta  > 0.40  n = %d  median %.3f%%  max %.3f%% (reported, not asserted)\n", nnz(~lo), 100 * median(err(~lo)), 100 * max(err(~lo)));
if max(err(lo)) > MAX_ERR_LOW_BETA
    error("testKConv:maxErr", "%s", sprintf("max relative error %.4f exceeds %.4f at beta <= 0.40", max(err(lo)), MAX_ERR_LOW_BETA));
end
if median(err(lo)) > MED_ERR_LOW_BETA
    error("testKConv:medErr", "%s", sprintf("median relative error %.4f exceeds %.4f at beta <= 0.40", median(err(lo)), MED_ERR_LOW_BETA));
end

%% 4. The retired constant
K13 = gridKConv_v001(1/3);
BQ = [0.28 0.30 1/3 0.40 2/3];                      % the values the manuscript quotes (Part 1 (c), #244)
fprintf("   K at quoted beta: %s\n", strjoin(compose("K(%.3f) = %.1f", [BQ; gridKConv_v001(BQ)]'), ", "));
gap = KCONV_RETIRED / K13 - 1;
fprintf("4. K(1/3) = %.4f; retired K_conv %.4f is %+.2f%% from it\n", K13, KCONV_RETIRED, 100 * gap);
if gap < RETIRED_GAP_EXPECT(1) || gap > RETIRED_GAP_EXPECT(2)
    error("testKConv:retired", "%s", sprintf("retired-constant gap %.4f outside expected [%.4f, %.4f]", gap, RETIRED_GAP_EXPECT));
end

%% Save
tbl = table(b(:), K(:), 'VariableNames', ["beta" "K"]);
runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "tbl", "K13", "gap", "runDate");
fprintf("PASS. Saved %s\n", OUT_MAT);
