function results = estimateNoiseGSJ_v001(trials, opts)
% estimateNoiseGSJ_v001 Difference-based noise sigma estimation (universal).
%
% Estimates noise standard deviation from any trajectory using the
% Gasser-Sroka-Jennen-Steinmetz (1986) second-order difference estimator:
%
%   sigma^2 = 0.5 * mean(diff(x, 2).^2)
%
% Works on any trajectory shape (periodic, aperiodic, single-trial) with
% no parameters to choose. Assumes noise is high-frequency relative to
% the signal, which is well satisfied for kinematic data at 60-240 Hz.
%
% Also computes the first-order Rice (1984) estimator for comparison:
%
%   sigma^2 = 0.5 * mean(diff(x, 1).^2)
%
% USAGE:
%   R = estimateNoiseGSJ_v001(trials)
%   R = estimateNoiseGSJ_v001(trials, Verbose=false)
%
% REFERENCES:
%   Gasser, Sroka & Jennen-Steinmetz (1986). Biometrics, 42(1), 87-101.
%   Rice (1984). Annals of Statistics, 12(4), 1215-1230.
%   Hall, Kay & Titterington (1990). Biometrika, 77(3), 521-528.
%
% PREREG REF: Section 7.3 (noise characterisation)
% Fraser, D.S. (2025)

    arguments
        trials (:,1) struct
        opts.Verbose (1,1) logical = true
    end

    nTrials = numel(trials);

    gsjSigmaX  = NaN(nTrials, 1);
    gsjSigmaY  = NaN(nTrials, 1);
    riceSigmaX = NaN(nTrials, 1);
    riceSigmaY = NaN(nTrials, 1);
    trialIDs   = strings(nTrials, 1);
    subjIDs    = strings(nTrials, 1);

    for ii = 1:nTrials
        t = trials(ii);
        x = double(t.x(:));
        y = double(t.y(:));
        trialIDs(ii) = string(t.trialID);
        subjIDs(ii)  = string(t.subjectID);

        % GSJ second-order difference estimator
        % Annihilates linear trends; robust to signal curvature
        d2x = diff(x, 2);
        d2y = diff(y, 2);
        gsjSigmaX(ii) = sqrt(0.5 * mean(d2x.^2));
        gsjSigmaY(ii) = sqrt(0.5 * mean(d2y.^2));

        % Rice first-order difference estimator (for comparison)
        d1x = diff(x, 1);
        d1y = diff(y, 1);
        riceSigmaX(ii) = sqrt(0.5 * mean(d1x.^2));
        riceSigmaY(ii) = sqrt(0.5 * mean(d1y.^2));
    end

    gsjSigmaMean  = mean([gsjSigmaX,  gsjSigmaY],  2);
    riceSigmaMean = mean([riceSigmaX, riceSigmaY], 2);

    results = table(trialIDs, subjIDs, ...
        gsjSigmaX, gsjSigmaY, gsjSigmaMean, ...
        riceSigmaX, riceSigmaY, riceSigmaMean, ...
        'VariableNames', ["trialID", "subjectID", ...
        "gsjSigmaX", "gsjSigmaY", "gsjSigmaMean", ...
        "riceSigmaX", "riceSigmaY", "riceSigmaMean"]);

    if opts.Verbose
        db = string(trials(1).database);
        fprintf("\n=== GSJ NOISE ESTIMATION: %s (%d trials) ===\n", db, nTrials);
        fprintf("GSJ sigma (2nd-order):  %.4f +/- %.4f (data units)\n", ...
            mean(gsjSigmaMean), std(gsjSigmaMean));
        fprintf("Rice sigma (1st-order): %.4f +/- %.4f (data units)\n", ...
            mean(riceSigmaMean), std(riceSigmaMean));
        fprintf("Rice/GSJ ratio:         %.2f (>1 means signal curvature contributes)\n", ...
            mean(riceSigmaMean) / mean(gsjSigmaMean));
    end

end
