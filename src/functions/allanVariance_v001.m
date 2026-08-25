function [avar, taus, alpha, alphaR2] = allanVariance_v001(x, fs, opts)
% allanVariance_v001 Overlapping Allan variance for noise colour estimation.
%
% Computes the overlapping Allan variance at logarithmically spaced
% averaging times tau, then fits log(AVAR) vs log(tau) to recover the
% noise spectral exponent alpha from the scaling relationship:
%
%   AVAR(tau) ~ tau^(alpha - 1)
%
% so alpha = slope + 1 in log-log space.
%
% White noise (alpha=0): slope = -1
% Pink noise  (alpha=1): slope =  0
% Red noise   (alpha=2): slope = +1
%
% USAGE:
%   [avar, taus, alpha, R2] = allanVariance_v001(x, fs)
%   [avar, taus, alpha, R2] = allanVariance_v001(x, fs, NTaus=30)
%
% INPUTS:
%   x    - 1D signal (column vector)
%   fs   - sampling rate (Hz)
%
% OUTPUTS:
%   avar    - Allan variance at each tau
%   taus    - averaging times (seconds)
%   alpha   - estimated noise spectral exponent
%   alphaR2 - R^2 of the log-log fit
%
% REFERENCES:
%   Allan (1966). Proceedings of the IEEE, 54(2), 221-230.
%   IEEE Standard 1139-2008.
%   El-Sheimy, Hou & Niu (2008). Sensors, 8(5), 3508-3547.
%
% Fraser, D.S. (2025)

    arguments
        x (:,1) double
        fs (1,1) double
        opts.NTaus (1,1) double = 30   % number of tau values
        opts.MinTau (1,1) double = 0   % 0 = auto (2/fs)
        opts.MaxTau (1,1) double = 0   % 0 = auto (N/5 / fs)
    end

    N = numel(x);
    x = x(:);

    % Tau range
    minM = 2;
    maxM = floor(N / 5);  % need at least 5 non-overlapping windows
    if opts.MinTau > 0, minM = max(minM, round(opts.MinTau * fs)); end
    if opts.MaxTau > 0, maxM = min(maxM, round(opts.MaxTau * fs)); end

    % Logarithmically spaced averaging lengths (in samples)
    mVals = unique(round(logspace(log10(minM), log10(maxM), opts.NTaus)));
    nTaus = numel(mVals);

    avar = NaN(nTaus, 1);
    taus = mVals(:) / fs;

    % Overlapping Allan variance
    for k = 1:nTaus
        m = mVals(k);
        % Cumulative sum for fast block averaging
        cs = [0; cumsum(x)];
        % Block averages of length m, overlapping by m-1
        nBlocks = N - m + 1;
        if nBlocks < 2, break; end
        blockAvg = (cs(m+1:N+1) - cs(1:N-m+1)) / m;
        % Allan variance = 0.5 * mean of squared consecutive differences
        avar(k) = 0.5 * mean(diff(blockAvg).^2);
    end

    % Remove any NaN taus (from break above)
    valid = ~isnan(avar);
    avar = avar(valid);
    taus = taus(valid);

    % Fit log-log slope
    if numel(taus) >= 3
        logTau  = log10(taus);
        logAvar = log10(avar);
        X = [ones(numel(logTau), 1), logTau(:)];
        b = X \ logAvar(:);
        slope = b(2);
        alpha = slope + 1;  % AVAR ~ tau^(alpha-1)
        ssRes = sum((logAvar(:) - X*b).^2);
        ssTot = sum((logAvar(:) - mean(logAvar)).^2);
        alphaR2 = 1 - ssRes / max(ssTot, eps);
    else
        alpha = NaN;
        alphaR2 = NaN;
    end

end
