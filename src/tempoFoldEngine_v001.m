function R = tempoFoldEngine_v001(J, cfg)
% tempoFoldEngine_v001  Six-pipeline beta recovery for a table of trajectory jobs.
% SPEC_TempoFoldSweep_v001 section 4.1. Replaces saturationSweepEngine_v001 for
% this sweep: SG is case 6 ([4 17], fs-scaled, as runner and v058), edge clip and
% regression seeding are per-job, pipeline labels are stored explicitly, and every
% failure is counted with its first message (v001 swallowed them in bare catch).
%
% J   table, one row per job. Required variables:
%       jobID, block (string), generator ("ellipse" | "v058"), a, b, f0, vgf,
%       FS, betaGen, sigma, alpha, nReps, edgeClip, nCycles, seedMode ("runner" | "truth")
%     generator "ellipse": generatePowerLawEllipse_v001(a, b, f0, FS, betaGen), fixed f0.
%     generator "v058": generateSyntheticData_v011(6, [1920;1080], FS, betaGen, vgf,
%       nCycles, 20, 0, 0) with yt = -yt, exactly as Toolchain_func_v032 L135-141;
%       fixed VGF, tempo realised. a, b, f0 are ignored.
%     seedMode "runner": OLS seeded [1 -1/3], LMLS/IRLS seeded from OLS
%       (runLoopClosureFftnoise_v012 L856-866). "truth": all seeded
%       [VGF_gen beta_gen] (Toolchain_func_v032 L107).
% cfg struct: nWorkers (0 = serial), rngSeed.
% R   J plus: M, nKept, lengthExcluded, f0Real, vgfGen, betaRec {nReps x 6},
%     betaMean, betaSD, betaAnalytic (1 x 6; ellipse only), nFail, firstErr.
%
% Length rule as the v058 harness: kept points = M - 2*edgeClip + 1 must be >= 10.
% Fraser, D.S. (2026)  v001

    req = ["jobID" "block" "generator" "a" "b" "f0" "vgf" "FS" "betaGen" "sigma" ...
           "alpha" "nReps" "edgeClip" "nCycles" "seedMode"];
    miss = setdiff(req, string(J.Properties.VariableNames));
    if ~isempty(miss), error("tempoFold:jobs", "%s", "J missing: " + strjoin(miss, ", ")); end

    n = height(J);
    out = cell(n, 1);
    Jc = table2struct(J);
    seed = cfg.rngSeed;
    parfor (k = 1:n, cfg.nWorkers)
        out{k} = runJob_local(Jc(k), seed);
    end
    R = [J, struct2table(vertcat(out{:}))];
end

%% =========================================================================
function o = runJob_local(j, rngSeed)
    LABELS = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];
    DERIV  = {2, [2 10 1]; 6, [4 17]};     % {filterType, filterParams}: BWFD, SG
    REG    = [3 4 5];                      % OLS, LMLS, IRLS
    MIN_PTS = 10;

    o = struct("M", NaN, "nKept", NaN, "lengthExcluded", false, "f0Real", NaN, ...
        "vgfGen", NaN, "betaRec", {NaN(j.nReps, 6)}, "betaMean", NaN(1, 6), ...
        "betaSD", NaN(1, 6), "betaAnalytic", NaN(1, 6), "nFail", zeros(1, 6), ...
        "firstErr", {strings(1, 6)}, "pipelines", {LABELS});

    %% Generate the noiseless trajectory
    ka = []; va = [];
    switch string(j.generator)
        case "ellipse"
            [x, y, ka, va] = generatePowerLawEllipse_v001(j.a, j.b, j.f0, j.FS, ...
                j.betaGen, 'nCycles', j.nCycles);
            o.vgfGen = va(1) * ka(1)^j.betaGen;   % v = VGF*kappa^-beta, exact by construction
            o.f0Real = j.f0;
        case "v058"
            [x, y] = generateSyntheticData_v011(6, [1920; 1080], j.FS, j.betaGen, ...
                j.vgf, j.nCycles, 20, 0, 0);
            x = x(:); y = -y(:);                  % Toolchain_func_v032 L141
            o.vgfGen = j.vgf;
            o.f0Real = j.nCycles / (numel(x) / j.FS);
        otherwise
            error("tempoFold:generator", "%s", "Unknown generator: " + string(j.generator));
    end
    x = x(:); y = y(:);
    o.M = numel(x);
    o.nKept = o.M - 2*j.edgeClip + 1;
    if o.nKept < MIN_PTS
        o.lengthExcluded = true;
        return
    end
    seed0 = [1, -1/3];
    if string(j.seedMode) == "truth", seed0 = [o.vgfGen, j.betaGen]; end

    %% Analytic level (exact speed and curvature; ellipse only)
    if ~isempty(ka)
        o.betaAnalytic = regressSix_local(va, ka, seed0, j.seedMode, REG, DERIV);
    end

    %% Noisy (or sigma = 0) replicates through the six pipelines
    if j.sigma > 0
        rs = RandStream("threefry4x64_20", "Seed", rngSeed);
        rs.Substream = j.jobID;
        RandStream.setGlobalStream(rs);
    end
    for rep = 1:j.nReps
        xs = x; ys = y;
        if j.sigma > 0
            xs = xs + reshape(generateCustomNoise_v003(o.M, j.alpha, j.sigma, j.FS), [], 1);
            ys = ys + reshape(generateCustomNoise_v003(o.M, j.alpha, j.sigma, j.FS), [], 1);
        end
        for d = 1:2
            cols = (d-1)*3 + (1:3);
            try
                [dx, dy] = differentiateKinematicsEBR(xs, ys, DERIV{d,1}, DERIV{d,2}, j.FS);
            catch ME
                o = logFail_local(o, cols, "differentiate: " + ME.message);
                continue
            end
            c  = j.edgeClip;
            vx = dx(c:end-c, 2); vy = dy(c:end-c, 2);
            ax = dx(c:end-c, 3); ay = dy(c:end-c, 3);
            sp = hypot(vx, vy);
            kp = curvatureKinematicEBR(vx, vy, ax, ay);
            ok = isfinite(sp) & isfinite(kp) & kp > 0;     % Toolchain_func_v032 L324
            if sum(ok) < MIN_PTS
                o = logFail_local(o, cols, sprintf("only %d valid points", sum(ok)));
                continue
            end
            [b, errs] = regressThree_local(sp(ok), kp(ok), seed0, j.seedMode, REG);
            o.betaRec(rep, cols) = b;
            for r = 1:3
                if errs(r) ~= "", o = logFail_local(o, cols(r), errs(r)); end
            end
        end
    end
    o.betaMean = mean(o.betaRec, 1, "omitnan");
    o.betaSD   = std(o.betaRec, 0, 1, "omitnan");
end

%% =========================================================================
function [b, errs] = regressThree_local(sp, kp, seed0, seedMode, REG)
    b = NaN(1, 3); errs = strings(1, 3); lm = seed0;
    for r = 1:3
        try
            [bb, vv] = regressDataEBR(sp, kp, REG(r), lm, 0, 0);
            b(r) = bb;
            if ~isfinite(bb), errs(r) = "non-finite beta"; end
            if r == 1 && string(seedMode) == "runner" && isfinite(bb) && isfinite(vv)
                lm = [vv, bb];                          % runner L866: seed LMLS/IRLS from OLS
            end
        catch ME
            errs(r) = "regress: " + ME.message;
        end
    end
end

function b6 = regressSix_local(va, ka, seed0, seedMode, REG, DERIV)
    % Analytic regression is derivation-free, so both derivation slots get the same values.
    b3 = regressThree_local(va(:), ka(:), seed0, seedMode, REG);
    b6 = repmat(b3, 1, size(DERIV, 1));
end

function o = logFail_local(o, cols, msg)
    o.nFail(cols) = o.nFail(cols) + 1;
    for c = cols
        if o.firstErr(c) == "", o.firstErr(c) = string(msg); end
    end
end
