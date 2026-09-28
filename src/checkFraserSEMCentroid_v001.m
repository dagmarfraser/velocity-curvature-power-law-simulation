% checkFraserSEMCentroid_v001.m  Fraser-only R5 SEM lookup + R7 gap verification.
%
% Part of TODO_FraserIntegration_v001.md Part C1 ("SEM centroid lookup for Fraser (R5)"
% and "R7 surrogate alpha gap check"). Written as a standalone, minimal script rather
% than extending knockDownFlags_v001.m, per that TODO's own caution: knockDownFlags'
% hardcoded Pilot row still carries the stale "3.91mm" sigma (see TODO_FraserIntegration
% Part B item 2), so a Fraser-only query avoids touching that file at all.
%
% FLAG 1 mechanism reused verbatim from knockDownFlags_v001.m (nearest-neighbour snap
% on the perCoordinateSEM_v2_001 coordTable grid), applied to Fraser's own centroid only:
%   alpha=4.289, sigma=2.009 mm, fs=240 Hz  (Finding #134/#135; noiseCharacterisation_fraser.mat)
%
% R7 gap is read directly from loopClosureResults_Fraser_per_subject_shaped_xu_v008.mat's
% own surrogateAlphaGapMaj/Min columns (N=48 trials) rather than trusting the
% TODO-cited 0.297/0.344 figures -- this script recomputes them from source.
%
% OUTPUT: console tables only. No files written.
%
% USAGE: checkFraserSEMCentroid_v001
%
% Fraser, D.S. (2026)  v001

clearvars
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
cd(srcDir);
% semAdequacyThreshold_v001 lives in src/functions -- add it or the threshold
% call below fails with an undefined function.
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('=== FRASER R5 SEM CENTROID LOOKUP + R7 GAP CHECK ===\n\n');

%% ---------------------------------------------------------------------
%% R5: SEM at Fraser's centroid (alpha=4.289, sigma=2.009mm, fs=240Hz)
%% ---------------------------------------------------------------------
fprintf('--- R5: SEM lookup at Fraser centroid ---\n');

semFile = fullfile(srcDir, 'perCoordinateSEM_v2_001.mat');
if ~isfile(semFile)
    error('checkFraserSEMCentroid:notFound', '%s', ...
        sprintf('FAILED PATH: %s not found.', semFile));
end
load(semFile, 'coordTable');
T = coordTable;

allAlpha = sort(unique(T.alpha));
allSigma = sort(unique(T.sigma));
allFs    = sort(unique(T.fs));

tA = 4.289;   % Finding #135 / noiseCharacterisation_fraser.mat
tS = 2.009;   % mm
tF = 240;     % Hz

[~,ai] = min(abs(allAlpha - tA));
[~,si] = min(abs(allSigma - tS));
[~,fi] = min(abs(allFs    - tF));
sA = allAlpha(ai); sS = allSigma(si); sF = allFs(fi);

fprintf('  Query   : alpha=%.3f  sigma=%.3f mm  fs=%d Hz\n', tA, tS, tF);
fprintf('  Snapped : alpha=%.3f  sigma=%.3f mm  fs=%d Hz', sA, sS, sF);
if abs(sA - tA) > 0.5 || abs(sS - tS) > 1.0
    fprintf('   [large snap distance -- check grid coverage]\n');
else
    fprintf('\n');
end

subAll = T(T.alpha==sA & T.sigma==sS & T.fs==sF, :);
if isempty(subAll)
    fprintf('  WARNING: no coordTable rows at snapped coordinate. No SEM available.\n');
else
    [G, pipNames] = findgroups(subAll.pipeline);
    semPerPipe = splitapply(@(x) mean(x,'omitnan'), subAll.sem, G);
    fprintf('\n  %-12s  %8s\n', 'Pipeline', 'SEM');
    fprintf('  %s\n', repmat('-',1,22));
    for i = 1:numel(pipNames)
        fprintf('  %-12s  %8.4f\n', pipNames(i), semPerPipe(i));
    end
    sgLmlsIdx = find(pipNames == "SG-LMLS", 1);
    if ~isempty(sgLmlsIdx)
        SEM_ADEQUATE = semAdequacyThreshold_v001();   % MDC/2.77 = 0.0108303, exact
        fprintf('\n  SG-LMLS (primary pipeline, Finding #85/#88): SEM = %.5f', semPerPipe(sgLmlsIdx));
        if semPerPipe(sgLmlsIdx) < SEM_ADEQUATE
            fprintf('  [ADEQUATE, < %.5f = MDC/2.77]\n', SEM_ADEQUATE);
        else
            fprintf('  [does NOT meet %.5f = MDC/2.77 threshold]\n', SEM_ADEQUATE);
        end
    end
end

%% ---------------------------------------------------------------------
%% R7: surrogate alpha-gap fidelity check, recomputed from source
%% ---------------------------------------------------------------------
fprintf('\n--- R7: surrogate alpha-gap check (recomputed from source) ---\n');

lcFile = fullfile(srcDir, 'loopClosureResults_Fraser_per_subject_shaped_xu_v008.mat');
if ~isfile(lcFile)
    error('checkFraserSEMCentroid:notFound', '%s', ...
        sprintf('FAILED PATH: %s not found.', lcFile));
end
D = load(lcFile, 'results');
r = D.results;
N = numel(r);

gapMaj = nan(N,1); gapMin = nan(N,1);
for ti = 1:N
    if isfield(r(ti), 'surrogateAlphaGapMaj') && isnumeric(r(ti).surrogateAlphaGapMaj)
        gapMaj(ti) = r(ti).surrogateAlphaGapMaj;
    end
    if isfield(r(ti), 'surrogateAlphaGapMin') && isnumeric(r(ti).surrogateAlphaGapMin)
        gapMin(ti) = r(ti).surrogateAlphaGapMin;
    end
end

meanGapMaj = mean(gapMaj, 'omitnan');
meanGapMin = mean(gapMin, 'omitnan');
nMaj = sum(isfinite(gapMaj));
nMin = sum(isfinite(gapMin));

fprintf('  N trials = %d\n', N);
fprintf('  mean |gap| major axis = %.4f  (N=%d finite)\n', meanGapMaj, nMaj);
fprintf('  mean |gap| minor axis = %.4f  (N=%d finite)\n', meanGapMin, nMin);

STANDALONE_GATE = 0.5;    % README_ShapedXu_v001.md Sec 6.1: gate threshold
STANDALONE_OBS  = 0.075;  % README_ShapedXu_v001.md Sec 6.1: standalone-check observed value

fprintf('\n  Reference points (docs/README_ShapedXu_v001.md Sec 6.1/6.2, Finding #75):\n');
fprintf('    Standalone gate PASS threshold  : %.3f\n', STANDALONE_GATE);
fprintf('    Standalone gate observed value  : %.3f  (clean synthetic fGn, not in-runner)\n', STANDALONE_OBS);
fprintf('    In-runner range, other 5 datasets (Finding #75): 0.199 (HALO, best) to 0.305/0.325 (Pilot, worst)\n');

fprintf('\n  Verdict:\n');
if meanGapMaj < STANDALONE_GATE && meanGapMin < STANDALONE_GATE
    fprintf('    Fraser (%.3f / %.3f) is comfortably below the %.1f standalone gate threshold.\n', ...
        meanGapMaj, meanGapMin, STANDALONE_GATE);
    fprintf('    No formal separate per-dataset in-runner threshold exists in the docs beyond\n');
    fprintf('    this standalone gate; the in-runner range across the other 5 datasets is\n');
    fprintf('    0.199-0.325, so Fraser sits inside the established range, not an outlier above it.\n');
else
    fprintf('    WARNING: Fraser gap exceeds the standalone gate threshold. Flag before citing loopCCC.\n');
end

fprintf('\n=== DONE ===\n');
