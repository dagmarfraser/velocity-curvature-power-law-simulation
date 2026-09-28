function T = computeInversionNoiseRepetitionFloor_v001(opts)
% COMPUTEINVERSIONNOISEREPETITIONFLOOR_V001  Theoretical lower bound on
% trial count from noise/sampling rate alone, via closed-form algebra on
% the ALREADY-COMPUTED per-coordinate SEM grid -- no new simulation.
%
% MOTIVATION (2026-09-20). Dagmar's question: can noise colour/magnitude
% and sampling rate be mapped to a minimal repetition count? Answer:
% only for the INVERSION-noise component (the smallest of the three
% terms Finding #216 decomposes, 10-20% typically) -- within-subject and
% between-subject variance are population properties no noise model can
% predict, and Finding #216 shows the latter dominates (62-83%). This
% function computes the inversion-noise floor only, precisely so its
% smallness (or, for harsher noise, its non-smallness) can be shown
% rather than asserted.
%
% METHOD, closed-form, no new analysis. perCoordinateSEM_v2_001.mat's
% own `sem` column is the standard error of the mean betaGenStar/betaRec
% recovery at a fixed (filterType, regressType, betaGen, VGF, fs, alpha,
% sigma), computed from nReps=5 simulated replicates at that exact grid
% coordinate. Standard SEM scaling, SEM_N = SD/sqrt(N), gives the
% per-single-replicate SD as SD = sem_5 * sqrt(5); solving for the N at
% which SEM_N first reaches a target then gives
%   N_required = (SD / target)^2 = 5 * (sem_5 / target)^2.
% This is arithmetic on an existing table (already cited throughout this
% paper for SEM-adequacy classification), not a new simulation or a new
% statistical method.
%
% GRID MISMATCH, stated not hidden: the empirical datasets' own native
% sampling rates (Fraser 240Hz -- an exact grid match; Cook 133Hz;
% Hickman/Zarandi/Dhieb ~100Hz) mostly do not coincide with the grid's
% own three levels (60/120/240Hz) -- this is the SAME mismatch already
% disclosed in Part 2 Method's "grid is not a cheaper alternative"
% passage, reused here rather than re-argued. Nearest available grid fs
% is used per dataset; Fraser alone gets an exact match.
%
% PIPELINE CHOICE, stated not defaulted silently: SG-IRLS is used as the
% pipeline for the headline table below because it was, empirically,
% the tightest or near-tightest pipeline for every core dataset in
% exploreTier3PipelineChoice_v001.m's own six-pipeline comparison -- an
% OPTIMISTIC bound, not a claim that SG-IRLS is the right choice in
% general (see that script's own header for why an "auto-pick the
% tightest pipeline" feature needs an identifiability safeguard this
% inversion-noise-only calculation does not by itself provide). All six
% pipelines are computed and returned for transparency, not just the
% headline one.
%
% VGF, marginalised: the grid crosses 14 VGF levels per (fs,alpha,sigma,
% pipeline,betaGen) coordinate; this function takes the median sem
% across VGF rather than picking one arbitrarily, since no single real
% dataset's own VGF was matched to a specific grid level for this
% illustration.
%
% SYNTAX:
%   T = computeInversionNoiseRepetitionFloor_v001()
%   T = computeInversionNoiseRepetitionFloor_v001(Pipeline="BWFD-OLS")
%
% INPUTS (name-value):
%   SrcDir     (1,1) string = "."      Folder containing perCoordinateSEM_v2_001.mat
%   TargetSEM  (1,1) double = 0.03/2.77  MDC/2.77, the paper's own SEM-adequacy threshold
%   TargetBeta (1,1) double = 1/3      Generating beta_gen to query (nearest grid point used)
%   Datasets   table = the seven canonical empirical (name, alpha, sigma, native fs) coordinates below
%
% OUTPUT: table T, one row per dataset x pipeline, columns:
%   dataset, alpha_real, sigma_real, fs_native, fs_grid, alpha_grid, sigma_grid,
%   pipeline, sem_N5, K_required
%
% See also: exploreTier3PipelineChoice_v001, tier3DesignPrecision_v001
%
% Fraser, D.S. (2026)

arguments
    opts.SrcDir (1,1) string = "."
    opts.TargetSEM (1,1) double {mustBePositive} = 0.03/2.77
    opts.TargetBeta (1,1) double = 1/3
    opts.Pipeline (1,1) string = "SG-IRLS"
    opts.AllPipelines (1,1) logical = true
end

% Seven canonical empirical coordinates, alpha/sigma as reported in this
% paper's own text (Part 2 Method's real-vs-shaped_xu comparison and the
% dataset-registration table), native fs as reported in the
% dataset-registration table (Fraser, Cook, Zarandi) or corrected
% 2026-09-20 by Dagmar's own direct knowledge and a Zotero fulltext check:
% Hickman shares Cook's own WACOM paradigm (133Hz, not the ~100Hz
% originally assumed by analogy to Zarandi/Dhieb's own devices -- a
% wrong assumption caught and corrected before this went into the
% manuscript); Dhieb confirmed at exactly 100Hz by reading Dhieb et al.
% (2022) Methods 2.3.2 directly ("sampled at 100 Hz", GENIUS MousePen
% i608X, not a WACOM device) via Zotero, not assumed by analogy.
names   = ["Fraser";"Zarandi";"Cook CTRL";"Cook ASD";"Hickman PLAC";"Hickman HALO";"Dhieb"];
alphas  = [4.289;    3.18;     4.77;       5.06;      5.34;          5.42;          2.50];
sigmas  = [2.009;    4.77;     8.15;       7.84;      7.17;          7.42;          7.50];
fsNative = [240;      100;      133;        133;       133;           133;           100];
fsNativeNote = ["exact grid match";"stated (Zarandi WACOM, 100Hz)";"stated (Cook WACOM, 133Hz)"; ...
    "stated (Cook WACOM, 133Hz)";"same paradigm as Cook, Dagmar 2026-09-20 (WACOM, 133Hz)";...
    "same paradigm as Cook, Dagmar 2026-09-20 (WACOM, 133Hz)";"confirmed via Zotero, Dhieb et al. 2022 Methods 2.3.2: sampled at 100Hz"];

D = load(fullfile(opts.SrcDir, 'perCoordinateSEM_v2_001.mat'));
Tgrid = D.coordTable;

if opts.AllPipelines
    pipelines = categories(Tgrid.pipeline);
else
    pipelines = {char(opts.Pipeline)};
end

rows = {};
for d = 1:numel(names)
    fsG = nearestValue_local(Tgrid.fs, fsNative(d));
    maskFs = Tgrid.fs == fsG;
    bgGridAtFs = unique(Tgrid.betaGen(maskFs));
    bgG = nearestValue_local(bgGridAtFs, opts.TargetBeta);
    alphaGridAtFs = unique(Tgrid.alpha(maskFs));
    aG = nearestValue_local(alphaGridAtFs, alphas(d));
    sigmaGridAtFs = unique(Tgrid.sigma(maskFs));
    sG = nearestValue_local(sigmaGridAtFs, sigmas(d));

    for p = 1:numel(pipelines)
        pp = string(pipelines{p});
        mask = maskFs & abs(Tgrid.betaGen-bgG)<1e-6 & abs(Tgrid.alpha-aG)<1e-9 ...
            & abs(Tgrid.sigma-sG)<1e-9 & Tgrid.pipeline==pp;
        rowsHere = Tgrid(mask,:);
        if isempty(rowsHere)
            error('computeInversionNoiseRepetitionFloor_v001:NoGridMatch', '%s', sprintf( ...
                'No grid rows found for %s at fs=%.0f, betaGen=%.4f, alpha=%.2f, sigma=%.3f, pipeline=%s.', ...
                names(d), fsG, bgG, aG, sG, pp));
        end
        semMed = median(rowsHere.sem);
        Kreq = 5 * (semMed / opts.TargetSEM)^2;
        rows(end+1,:) = {names(d), alphas(d), sigmas(d), fsNative(d), fsNativeNote(d), ...
            fsG, aG, sG, pp, semMed, Kreq}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'dataset','alpha_real','sigma_real','fs_native', ...
    'fs_native_note','fs_grid','alpha_grid','sigma_grid','pipeline','sem_N5','K_required'});

fprintf('Target SEM = %.4f (MDC/2.77)\n', opts.TargetSEM);
fprintf('%-14s %-10s %-12s %-10s\n', 'dataset', 'pipeline', 'sem_N5', 'K_required');
for i = 1:height(T)
    fprintf('%-14s %-10s %-12.6f %-10.4f\n', T.dataset(i), T.pipeline(i), T.sem_N5(i), T.K_required(i));
end
end

function v = nearestValue_local(candidates, target)
    [~, idx] = min(abs(candidates - target));
    v = candidates(idx);
end
