function G = gridShapeGeometry_v002(root)
% gridShapeGeometry_v002  Single source for the v058 grid's trajectory geometry.
%   G = gridShapeGeometry_v002(root) returns the nu = 2 Huh & Sejnowski shape exactly as
%   generateSyntheticData_v011 builds it (GENERATE_ALL_Curves_SIMPLE L17: 100*[x+2; -y]),
%   in canvas px, plus the harness's px-per-mm constant.
%   Fields: pathPx (N x 2, open loop), kappaPx (N x 1, 1/px, Menger as the generator),
%   aPx, bPx (semi-extents), ab, aFitPx, bFitPx, abFit (least-squares ellipse fit, the
%   ellipse stand-in used by the tempo sweeps; Finding #236), pixelScale, aMM, bMM,
%   aFitMM, bFitMM, perimPx, source, cacheCheck.
%   The mm values exist only through pixelScale, which the harness applies to sigma
%   alone (Toolchain_caller_v058.m L213, L545); the generator itself has no mm concept.
% v002: adds the least-squares ellipse fit (235.06 x 108.85 px; checkGridShapeUnits_v001).
% Fraser, D.S. (2026)  v002
    arguments
        root (1,1) string = fileparts(fileparts(mfilename("fullpath")))
    end
    fn = fullfile(root, "src", "functions");
    src = fullfile(fn, "AccuracyMeasure", "Sinusoidal_Curves.mat");
    if ~isfile(src), error("gridShape:source", "%s", "Missing: " + src); end
    S = load(src, "path_all");
    p = S.path_all{6};                                   % nu = '2' (generateSyntheticData_v011 index 6)
    P = (100 * [p(1,:) + 2; -p(2,:)])';

    % pixelScale: parse the caller rather than restate it, so a change there fails loud here
    caller = fullfile(root, "src", "Toolchain_caller_v058.m");
    if ~isfile(caller), error("gridShape:caller", "%s", "Missing: " + caller); end
    tok = regexp(fileread(caller), "cfg\.pixelScale\s*=\s*([\d.]+)\s*/\s*([\d.]+)\s*;", "tokens", "once");
    if isempty(tok), error("gridShape:pixelScale", "%s", "cfg.pixelScale not found in " + caller); end
    pixelScale = str2double(tok{1}) / str2double(tok{2});

    % Cross-check against the cached baseline the grid actually loaded (resampled spline)
    c = dir(fullfile(fn, "baselineShp6_*Hz.mat"));
    if isempty(c), error("gridShape:cache", "%s", "No baselineShp6_*Hz.mat in " + fn); end
    ext = zeros(numel(c), 2);
    for i = 1:numel(c)
        B = load(fullfile(c(i).folder, c(i).name), "pathXYresample", "K");
        ext(i, :) = range(B.pathXYresample) / 2;
        if i == 1, Q = B.pathXYresample; K = B.K(:); end
    end
    if max(abs(ext - ext(1, :)), [], "all") > 1e-6
        error("gridShape:cacheMismatch", "%s", "baselineShp6 caches disagree across fs");
    end
    if max(abs(ext(1, :) - range(P) / 2)) > 0.5         % spline vs raw vertices, px
        error("gridShape:rawVsCache", "%s", sprintf("Raw %.2f x %.2f vs cache %.2f x %.2f px", ...
            range(P) / 2, ext(1, :)));
    end
    if norm(Q(end, :) - Q(1, :)) < 1e-9, Q(end, :) = []; K(end) = []; end   % open the loop
    if numel(K) ~= size(Q, 1), error("gridShape:K", "%s", "Curvature and path lengths differ"); end

    G = struct();
    G.pathPx     = Q;
    G.kappaPx    = K;
    G.aPx        = ext(1, 1);   G.bPx = ext(1, 2);   G.ab = G.aPx / G.bPx;
    G.pixelScale = pixelScale;
    G.aMM        = G.aPx / pixelScale;   G.bMM = G.bPx / pixelScale;
    fit          = fitEllipseAxes(Q);
    G.aFitPx     = fit(1);   G.bFitPx = fit(2);   G.abFit = G.aFitPx / G.bFitPx;
    G.aFitMM     = G.aFitPx / pixelScale;   G.bFitMM = G.bFitPx / pixelScale;
    G.perimPx    = sum(vecnorm(diff([Q; Q(1, :)]), 2, 2));
    G.source     = string(src);
    G.cacheCheck = string({c.name});
end

function ax = fitEllipseAxes(P)
% Axis-aligned least-squares ellipse about the centroid: (x/a)^2 + (y/b)^2 = 1.
    X = P - mean(P, 1);
    w = [X(:, 1).^2, X(:, 2).^2] \ ones(size(X, 1), 1);
    if any(w <= 0), error("gridShape:fit", "%s", "Ellipse fit not positive definite"); end
    ax = 1 ./ sqrt(w');
end
