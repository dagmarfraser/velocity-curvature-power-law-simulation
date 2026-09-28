function results = runLoopClosureAll_v015(opts)
%RUNLOOPCLOSUREALL_V015  Run the v015 loop closure over the datasets, interactively.
%
%   For an interactive MATLAB session with a 72-worker pool already up.
%   No batch plumbing: no Slurm, no log files, no scheduler assumptions.
%
%   v015 is v013 plus two persisted fields -- betaRecCurveMed / betaRecCurveSD,
%   the per-trial forward map the inversion is performed against, which v013
%   computes and discards. That map is what the out-of-sample form of
%   leave-one-pipeline-out needs and cannot get from any existing corpus.
%
%   THE ONE THING THIS ADDS BEYOND A FOR LOOP: v015 is supposed to change no
%   numerics, so betaObs / betaGenStar / ciLo / ciHi / invertStatus must match
%   v013 exactly. Verify() checks that per dataset and errors on any real
%   difference, because a difference means a generator bug, not a new result.
%
%   USAGE, on the HPC in an interactive MATLAB session:
%
%       runLoopClosureAll_v015
%
%   That is the whole thing. It generates v015 from v013 if absent, opens a
%   72-worker pool if none is up, runs the seven datasets, and verifies each
%   against v013.
%
%   Datasets run cheapest-first, Fraser last (~72% of the total). The runner
%   skips any dataset whose output already exists, so this is restartable.
%   Benchmarked ~586 core-hours total, regression-bound.
%
%   Provenance: written 2026-09-16, Session 110.

arguments
    % NOTE: these must match runLoopClosureFftnoise's own switch, which uses
    % SPACES, not underscores (v013 lines 178-216). The corpus FILENAMES use
    % underscores -- the runner converts via strrep(DATASET,' ','_') at line
    % 227. verify_local below applies the same conversion.
    opts.Datasets   (1,:) string = ["Zarandi","Cook CTRL","Cook ASD", ...
                                    "Hickman PLAC","Hickman HALO","Dhieb","Fraser"]
    opts.NoiseModel (1,1) string = "shaped_xu"
    opts.N_BETA     (1,1) double = 25
    opts.N_REPS     (1,1) double = 200
    opts.N_R7       (1,1) double = 10
    opts.RngSeed    (1,1) double = 1729
    opts.NWorkers   (1,1) double = 72
    opts.Verify     (1,1) logical = true
    % Tolerance on the STOCHASTIC fields. betaGenStar is not run-reproducible
    % even at a fixed seed: the synthetic side's RNG streams depend on the
    % parallel configuration, so an August 72-worker run and a September one
    % differ. Measured on Zarandi v013 vs v015: median 0.0016, p90 0.0098,
    % max 0.0437, 0.7% of cells above MDC, NaN flips 12/798 all in LMLS
    % (fitnlm convergence at the margin). See Finding #204.
    opts.TolMedian  (1,1) double = 0.005    % median |diff|, betaGenStar/ciLo/ciHi
    opts.TolNaNFrac (1,1) double = 0.05     % fraction of cells whose NaN state may flip
end

srcDir = fileparts(mfilename("fullpath"));
addpath(srcDir);

if ~isfile(fullfile(srcDir,"runLoopClosureFftnoise_v015.m"))
    fprintf("v015 not present -- generating from v013.\n");
    makeRunLoopClosureV015_v001("SrcDir", srcDir);
else
    fprintf("v015 present.\n");
end

% IdleTimeout must be disabled. verify_local runs single-threaded between
% datasets (17 min for Cook CTRL), the pool idles through it, and on
% 2026-09-16 it timed out mid-run -- after which v015 line 374's
% gcp('nocreate').NumWorkers hit a shutting-down pool and threw
% "Unrecognized method, property, or field 'NumWorkers'". That line is
% inherited from v013 and assumes a live pool, so the pool must not die.
p = ensurePool_local(opts.NWorkers);

results = struct("dataset",{},"nTrials",{},"minutes",{},"verify",{});
tAll = tic;

for ds = opts.Datasets
    fprintf("\n=== %s ===\n", ds);
    p = ensurePool_local(opts.NWorkers);   % re-open if a previous dataset's
    t0 = tic;                              % verify step let it lapse
    r = runLoopClosureFftnoise_v015(char(ds), ...
        'TrialSelection','all', 'NoiseModel',char(opts.NoiseModel), ...
        'N_BETA',opts.N_BETA, 'N_REPS',opts.N_REPS, 'N_R7',opts.N_R7, ...
        'RngSeed',opts.RngSeed, 'UseParfor',true);
    mins = toc(t0)/60;

    v = "skipped";
    if opts.Verify
        v = verify_local(srcDir, ds, opts.NoiseModel, r, opts);
    end
    fprintf("%s: %d trials, %.1f min, verify=%s\n", ds, numel(r), mins, v);
    results(end+1) = struct("dataset",ds,"nTrials",numel(r),"minutes",mins,"verify",v); %#ok<AGROW>
end

fprintf("\n=== done in %.1f h ===\n", toc(tAll)/3600);
fprintf("%-16s %8s %9s  %s\n","dataset","nTrials","minutes","verify");
for i = 1:numel(results)
    fprintf("%-16s %8d %9.1f  %s\n", results(i).dataset, results(i).nTrials, ...
        results(i).minutes, results(i).verify);
end
end

% =========================================================================
function p = ensurePool_local(nWorkers)
% Return a live pool with IdleTimeout disabled, creating one if needed.
% Fail Loud: if the pool cannot be established, stop -- do not silently fall
% back to serial execution, which would turn an 8-hour run into a 24-day one.
p = gcp("nocreate");
if ~isempty(p)
    try
        alive = p.Connected && p.NumWorkers > 0;
    catch
        alive = false;          % shutting down: property access throws
    end
    if ~alive
        fprintf("pool is dead or shutting down -- discarding it.\n");
        try, delete(p); catch, end
        p = [];
    end
end

if isempty(p)
    fprintf("opening pool (%d workers, IdleTimeout disabled)...\n", nWorkers);
    p = parpool("local", nWorkers, "IdleTimeout", Inf);
else
    try
        if ~isinf(p.IdleTimeout)
            p.IdleTimeout = Inf;
            fprintf("existing pool: IdleTimeout disabled.\n");
        end
    catch ME
        warning("runAllV015:idleTimeout", ...
            "Could not disable IdleTimeout on the existing pool (%s). The pool " + ...
            "may lapse during the between-dataset verify step.", ME.message);
    end
end

if isempty(p) || ~p.Connected
    error("runAllV015:noPool", "Could not establish a parallel pool.");
end
fprintf("pool: %d workers\n", p.NumWorkers);
end

function status = verify_local(srcDir, ds, model, r15, opts)
% What v015 must reproduce EXACTLY: betaObs. It is computed from the real
% trajectory with no stochastic component, so any difference there is a
% genuine defect and errors.
%
% What it CANNOT reproduce exactly: betaGenStar, ciLo, ciHi, invertStatus.
% These come from the synthetic forward map, whose RNG streams depend on the
% parallel configuration. Same seed, different run, different realisation.
% These are checked against a tolerance and REPORTED, not asserted equal.
% Asserting equality here halted an 8-hour run over a 0.0016 median on
% 2026-09-16; see Finding #204.
tag = strrep(ds, " ", "_");
f = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_v013.mat", tag, model));
if ~isfile(f)
    f = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_v014.mat", tag, model));
end
if ~isfile(f), status = "no-comparator"; return; end

S = load(f,"results"); r13 = S.results;
if numel(r13) ~= numel(r15)
    error("runAllV015:count","%s: %d trials vs comparator's %d.", ds, numel(r15), numel(r13));
end
if ~isfield(r15,"betaRecCurveMed")
    error("runAllV015:noCurve","%s: no betaRecCurveMed -- v015's entire purpose.", ds);
end

% ---- betaObs: must be exact -------------------------------------------
worstObs = 0;
for ti = 1:numel(r13)
    a = r13(ti).betaObs(:); b = r15(ti).betaObs(:);
    if ~isequal(isfinite(a), isfinite(b))
        error("runAllV015:obsNaN","%s trial %d: betaObs NaN pattern differs. " + ...
            "betaObs is deterministic -- this is a real defect.", ds, ti);
    end
    m = isfinite(a);
    if any(m), worstObs = max(worstObs, max(abs(a(m)-b(m)))); end
end
if worstObs > 1e-12
    error("runAllV015:obsDiff","%s: betaObs differs by %.3g. betaObs is " + ...
        "deterministic -- v015 must reproduce it exactly.", ds, worstObs);
end

% ---- stochastic fields: tolerance, and report ---------------------------
d = []; nanFlip = 0; nCell = 0;
for ti = 1:numel(r13)
    for fld = ["betaGenStar","ciLo","ciHi"]
        a = r13(ti).(fld)(:); b = r15(ti).(fld)(:);
        nCell = nCell + numel(a);
        nanFlip = nanFlip + nnz(isfinite(a) ~= isfinite(b));
        m = isfinite(a) & isfinite(b);
        if any(m), d = [d; abs(a(m)-b(m))]; end %#ok<AGROW>
    end
end
statFlip = 0; statN = 0;
for ti = 1:numel(r13)
    x = string(r13(ti).invertStatus); y = string(r15(ti).invertStatus);
    statN = statN + numel(x);
    statFlip = statFlip + nnz(x ~= y);
end

medD = median(d); p90 = quantile(d,0.90); maxD = max(d);
fracNaN = nanFlip / max(nCell,1);

if medD > opts.TolMedian
    error("runAllV015:stochDrift", ...
        "%s: median |diff| on the stochastic fields is %.4f, above TolMedian=%.4f. " + ...
        "Expected run-to-run noise is ~0.0016. This is larger than resampling " + ...
        "explains -- check the generator before trusting v015.", ds, medD, opts.TolMedian);
end
if fracNaN > opts.TolNaNFrac
    error("runAllV015:nanDrift", ...
        "%s: %.1f%% of cells flipped NaN state, above TolNaNFrac=%.1f%%.", ...
        ds, 100*fracNaN, 100*opts.TolNaNFrac);
end

status = sprintf("obs exact; stoch med %.4f p90 %.4f max %.4f, NaN-flip %d/%d, status-flip %d/%d", ...
    medD, p90, maxD, nanFlip, nCell, statFlip, statN);
end
