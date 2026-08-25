function exclusionMask = computeReversalExclusion_v001(position3D, samplingRate)
% COMPUTEREVERSALEXCLUSION_V001  Reversal-point exclusion via the
% separation-ratio criterion.
%
% Extracted verbatim (logic unchanged, only renamed/re-styled per this
% project's MATLAB coding guideline) from computeReversalExclusion_local
% in checkCurvatureBetaPropagationUnified_v001.m, so the same validated
% criterion can be reused elsewhere without duplicating the logic by hand.
% Validated empirically in checkCurvatureAtReversals_v001.m: the outbound-
% to-inbound path separation ratio at each reversal, normalised by
% movement amplitude, predicts curvature instability almost
% deterministically (ULTRA-MoCap, r=-0.958; Santos gait, r=-0.791).
% README_CoherenceGap_v003.md Section 4.2 documents the full derivation
% and the rejected alternative fixes.
%
% Window scales with sampling rate (150 ms each side) rather than being a
% fixed sample count, so it behaves consistently across datasets with
% different sampling rates. Threshold is the bottom quartile of that
% trial's own distribution of separation ratios, not a fixed absolute
% value across trials or datasets.
%
%% inputs
% position3D   - N-by-3 position trajectory, [x y z]
% samplingRate - Hz
%
%% outputs
% exclusionMask - N-by-1 logical, true where a sample falls within the
%                 exclusion window of one of the most cusp-like quartile
%                 of reversal points in this trial
%
% Created 2026-07-21 (extraction of existing validated logic)
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

arguments
    position3D   (:,3) double
    samplingRate (1,1) double
end

numberOfSamples = size(position3D, 1);
windowSamples = max(3, round(0.15 * samplingRate));
positionY = position3D(:, 2);
velocity1D = gradient(positionY) * samplingRate;
signChanges = find(diff(sign(velocity1D)) ~= 0);
amplitude = max(positionY) - min(positionY);

exclusionMask = false(numberOfSamples, 1);
if numel(signChanges) < 3
    return
end

separationRatios = nan(numel(signChanges) - 2, 1);
for reversalIndex = 2:(numel(signChanges) - 1)
    reversalSample = signChanges(reversalIndex);
    lowSample = max(1, reversalSample - windowSamples);
    highSample = min(numberOfSamples, reversalSample + windowSamples);
    outboundPath = position3D(lowSample:reversalSample, :);
    inboundPath = position3D(reversalSample:highSample, :);
    minimumDistances = nan(size(outboundPath, 1), 1);
    for outboundIndex = 1:size(outboundPath, 1)
        distances = sqrt(sum((inboundPath - outboundPath(outboundIndex, :)) .^ 2, 2));
        minimumDistances(outboundIndex) = min(distances);
    end
    separationRatios(reversalIndex - 1) = median(minimumDistances) / amplitude;
end

exclusionThreshold = prctile(separationRatios, 25);
for reversalIndex = 2:(numel(signChanges) - 1)
    if separationRatios(reversalIndex - 1) <= exclusionThreshold
        reversalSample = signChanges(reversalIndex);
        lowSample = max(1, reversalSample - windowSamples);
        highSample = min(numberOfSamples, reversalSample + windowSamples);
        exclusionMask(lowSample:highSample) = true;
    end
end

end
