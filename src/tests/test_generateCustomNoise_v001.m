function test_generateCustomNoise_v001()
% test_generateCustomNoise_v001  Verify generateCustomNoise_v003 spectral slope accuracy.
%
% For each alpha in {0, 0.5, 1, 2, 3}, generates N_REPS realisations and
% estimates the spectral exponent via pmtm log-log regression. Pass criterion:
% median |alpha_out - alpha_in| < SLOPE_TOL across all tested colours.
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
N        = 20000;       % samples per realisation (long for stable spectral estimate)
FS       = 100;         % Hz
AMP      = 1.0;         % unit amplitude (std dev)
N_REPS   = 20;          % realisations per alpha
FMIN     = 0.5;         % Hz -- lower bound for slope regression
FMAX     = 30;          % Hz -- upper bound (well below Nyquist)
SLOPE_TOL = 0.25;       % |alpha_out - alpha_in| tolerance

alphaVals = [0, 0.5, 1.0, 2.0, 3.0];

fprintf('=== test_generateCustomNoise_v001 ===\n');
fprintf('N=%d, fs=%d Hz, reps=%d per alpha, fit [%.1f %.1f] Hz, tol=%.2f\n\n', ...
    N, FS, N_REPS, FMIN, FMAX);
fprintf('  %-8s  %-10s  %-10s  %-10s  %s\n', 'alpha_in', 'median_out', 'mean_out', 'max_err', 'verdict');
fprintf('  %s\n', repmat('-', 1, 60));

allPass = true;

for ai = 1:numel(alphaVals)
    alphaIn  = alphaVals(ai);
    alphaEst = nan(N_REPS, 1);

    for r = 1:N_REPS
        try
            noise = generateCustomNoise_v003(N, alphaIn, AMP, FS);
            alphaEst(r) = estimateSpectralAlpha(noise, FS, FMIN, FMAX);
        catch ME
            warning('test_generateCustomNoise:genFail', 'rep %d: %s', r, ME.message);
        end
    end

    medErr  = abs(median(alphaEst, 'omitnan') - alphaIn);
    medOut  = median(alphaEst, 'omitnan');
    meanOut = mean(alphaEst, 'omitnan');
    maxErr  = max(abs(alphaEst - alphaIn), [], 'omitnan');

    pass = (medErr < SLOPE_TOL);
    if ~pass, allPass = false; end

    fprintf('  %-8.2f  %-10.3f  %-10.3f  %-10.3f  %s\n', ...
        alphaIn, medOut, meanOut, maxErr, verdict(pass));
end

fprintf('\n');
if allPass
    fprintf('PASS: all spectral slopes within tolerance %.2f\n', SLOPE_TOL);
else
    fprintf('FAIL: one or more spectral slopes exceeded tolerance %.2f\n', SLOPE_TOL);
end

end

% =========================================================================
function alphaHat = estimateSpectralAlpha(noise, fs, fMin, fMax)
% Estimate spectral exponent alpha from log-log pmtm slope.
% PSD ~ f^(-alpha), so slope of log(PSD) vs log(f) is -alpha.

[pxx, f] = pmtm(noise, 4, [], fs);

idx = (f > fMin) & (f <= fMax);
if sum(idx) < 5
    error('estimateSpectralAlpha:tooFewPoints', 'Fewer than 5 frequency bins in fit range');
end

p        = polyfit(log(f(idx)), log(pxx(idx)), 1);
alphaHat = -p(1);   % PSD ~ f^(-alpha) => slope = -alpha
end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
