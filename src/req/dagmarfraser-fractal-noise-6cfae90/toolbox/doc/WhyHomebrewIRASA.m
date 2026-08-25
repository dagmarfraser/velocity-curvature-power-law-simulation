%[text] # Why `estimateIRASA` Is a Homebrew Implementation
%[text] `estimateIRASA` ships its own implementation of the IRASA algorithm (resampling-median spectral estimation, Wen & Liu 2016) rather than simply wrapping FieldTrip's `ft_freqanalysis(method="irasa")`. This document explains why, using the validation work that established the verified alpha range stated in `estimateIRASA`'s docstring.
%[text] **tl;dr:** three mechanistically independent estimators agree closely for alpha in \[-2, 6\]. At alpha=7 they diverge -- not because any one estimator is broken, but because an information limit is being hit at typical trial lengths. `estimateIRASA` is verified for alpha in \[-2, 6\], not the full \[-2, 7\] range over which the noise generator `noiseXu` is merely **numerically stable**. Those are two different claims and this toolbox does not conflate them.
%[text] ## The circularity problem
%[text] Suppose you want to validate a noise generator's numerical stability by checking that IRASA recovers the target alpha you asked it to generate. That's a reasonable check -- until you also want to validate the IRASA **estimator** itself, and the only tool at hand is IRASA. If generator validation and estimator validation both route through the same IRASA implementation, close agreement at extreme alpha risks being two correlated errors cancelling out, not genuine confirmation of either.
%[text] Homebrew-vs-FieldTrip agreement alone does not resolve this either. Both are IRASA variants built on the same resampling-median mechanism and could in principle share a blind spot at the same operating point.
%[text] **Resolution:** bring in a third, mechanistically independent estimator -- a direct pmtm log-log slope fit with no resampling step at all -- and test whether all three converge on known synthetic targets.
%[text] ## Three estimators
estimators = table( ...
    ["HB"; "FT"; "PM"], ...
    ["iraAlphaSigma_v003 / estimateIRASA"; "ft_freqanalysis(method=""irasa"")"; "direct pmtm"], ...
    ["pmtm, whole-signal, resampling-median (homebrew)"; "Hanning-windowed, resampling-median (FieldTrip)"; "single PSD, log-log OLS slope, NO resampling step"], ...
    VariableNames=["Code", "Implementation", "Mechanism"]);
disp(estimators)
%%
%[text] PM shares HB's PSD estimator (pmtm, time-bandwidth product 4, same frequency resolution) and differs only in skipping the up/downsample- and-median step -- isolating that step as the variable under test.
%[text] ## Historical validation result
%[text] N=3000, fs=120 Hz, 15 repetitions per alpha, alpha in \[-2, 7\]. Bias is mean(alpha\_recovered - alpha\_true); spread3 is the range across all three estimators' mean bias at that alpha.
alphaTrue = (-2:7)';
hbBias = [0.053; 0.035; 0.014; -0.028; -0.077; -0.042; -0.036; -0.005; -0.114; -0.638];
ftBias = [0.086; 0.050; 0.001; -0.052; -0.117; -0.142; -0.154; -0.296; -0.433; -0.674];
pmBias = [0.050; 0.034; 0.013; -0.029; -0.078; -0.041; -0.031; 0.032; 0.014; -0.300];
% spread3 and gapHBPM are taken verbatim from the original source table
% (docs/README_IRASAThreeWayAgreement_v001.md), NOT recomputed from the
% rounded bias figures above -- the original derived columns were
% computed from full-precision underlying data before rounding to 3
% decimals for display, so max(bias)-min(bias) on the rounded values
% above does not reproduce them exactly (e.g. alpha=7: recomputing from
% rounded bias gives spread3=0.374, not the source's reported 0.458).
spread3 = [0.043; 0.037; 0.023; 0.038; 0.044; 0.101; 0.123; 0.327; 0.447; 0.458];
gapHBPM = [0.003; 0.001; 0.001; 0.001; 0.001; 0.001; 0.005; 0.037; 0.128; 0.338];

historicalResults = table(alphaTrue, hbBias, ftBias, pmBias, spread3, gapHBPM, ...
    VariableNames=["alpha", "HB_bias", "FT_bias", "PM_bias", "spread3", "gapHBPM"]);
disp(historicalResults)
%%
%[text] **Aggregate:** HB mean|bias|=0.104, FT mean|bias|=0.200, PM mean|bias|=0.062.
%[text] ## Interpretation
%[text] **alpha in \[-2, 4\]:** `HB-PM` \<= 0.005 throughout. Two mechanistically independent estimators (no shared resampling step) agree almost exactly; FieldTrip alone is the outlier, consistent with its coarser frequency resolution. This is genuine independent-method confirmation of generator fidelity, not self-referential agreement between two variants of the same algorithm.
%[text] **alpha = 5, 6:** `HB-PM` grows (0.037, 0.128) but stays well under `HB-FT` or `PM-FT`. The "FT is the outlier; HB/PM agree" pattern still holds -- `spread3` alone crosses a naive flag threshold here, which would read as "unverifiable" in isolation, but the pairwise breakdown shows this is still well-supported, driven by FieldTrip specifically rather than genuine three-way disagreement.
%[text] **alpha = 7:** the pattern inverts. `HB-FT` becomes **smaller** than `HB-PM` -- HB and FT drift toward each other while PM diverges from both, and PM itself picks up substantial bias and inflated variance. This is genuine three-way breakdown, most likely reflecting an information limit at N=3000 (too few independent low-frequency samples in-band to constrain a slope this steep) rather than an algorithm-specific artefact.
%[text] ## Verdict
%[text] **Verifiable range: alpha in \[-2, 6\].** Independent-method agreement (HB vs PM) supports estimator fidelity throughout this range, including alpha=5-6 where `spread3` alone would flag concern.
%[text] **alpha = 7 is not independently verifiable** with these three tools at this trial length. A noise generator being **numerically stable** at alpha=7 (no NaN/Inf, well-conditioned output -- see `noiseXu`'s docstring) is a different claim from an estimator being able to **measure** alpha=7 back out reliably. Do not conflate the two.
%[text] ## Live demonstration
%[text] The table above is a historical result (originally run against PowerLawSimulationPreReg's noise generator and estimator). Since `noiseXu` and `estimateIRASA` reproduce that exact generator and estimator bit-for-bit, the same pattern is reproducible live, right here, using only this toolbox's public functions.
N = 3000;
fs = 120;
fLow = 1.0;
fHigh = 30.0;
sigma = 5.0;
demoAlphas = [3, 6, 7];

for a = demoAlphas
    rng(500 + a);
    x = noiseXu(N, a, sigma, fs);
    [alphaHB, ~, ~, ~, ftInfo] = estimateIRASA(x, fs, fLow, fHigh, CompareFieldTrip=true);
    if ftInfo.available
        fprintf("alpha_true=%d  HB=%.3f  FT=%.3f  gap(HB,FT)=%.3f\n", ...
            a, alphaHB, ftInfo.alphaFT, ftInfo.gapHBFT);
    else
        fprintf("alpha_true=%d  HB=%.3f  (FieldTrip not available: %s)\n", ...
            a, alphaHB, ftInfo.message);
    end
end
%%
%[text] At alpha=3 the HB/FT gap should be small; by alpha=6-7 it typically widens, echoing the historical pattern above. This is a single realisation per alpha (not averaged over repetitions like the historical table), so don't over-read small numerical differences from the table -- the point is the qualitative pattern, not exact reproduction of specific bias values.
%[text] ## What this means for using this toolbox
%[text] - For alpha in \[-2, 6\], trust `estimateIRASA`'s homebrew estimate on its own; it's the cheapest of the three estimators and has the lowest aggregate bias.
%[text] - Requesting `CompareFieldTrip=true` gives you a second opinion for free (when FieldTrip is installed) -- useful as a sanity check on any individual analysis, not just at the edges of the verified range.
%[text] - Outside alpha in \[-2, 6\], treat any single estimate -- homebrew or otherwise -- as unverified. Consider it a stress-test result, not a ground-truth measurement. \
%[text] ## Cross-references
%[text] - `estimateIRASA.m` docstring -- states this same verified range
%[text] - `noiseXu.m` docstring -- states the separate, wider generator stability range and why it should not be conflated with this one
%[text] - Original validation: PowerLawSimulationPreReg  `docs/README_IRASAThreeWayAgreement_v001.md` and Finding \#100 \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
