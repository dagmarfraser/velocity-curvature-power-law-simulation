%[text] # `checkSpectralKnee` Worked Example: Real Motivation, Synthetic Ground Truth
%[text] This example combines two things: the actual real-data result that motivated `checkSpectralKnee` (Preston, Smith & Voytek's 2026 Lorentzian knee diagnostic, applied to five kinematic datasets from PowerLawSimulationPreReg -- summary statistics only, no raw trial data shipped with this toolbox), and a live demonstration against **synthetic data with an exactly known ground-truth knee frequency and alpha**, which the real data alone cannot provide.
%[text] ## Background: why check for a knee at all
%[text] IRASA (used by `estimateIRASA`) assumes its input is scale-free -- a pure power-law spectrum with no preferred timescale. Real electrophysiological and, potentially, kinematic signals can instead follow a **Lorentzian** form: flat below a "knee frequency" `f_knee`, decaying as a power law above it. If a knee falls inside the fitting band, IRASA's alpha estimate is biased by the flat-to-decaying transition it's straddling.
%[text] `checkSpectralKnee` answers two questions for a single residual: (1) is there statistical support for a knee at all (via AIC comparison against a plain fixed power-law model), and (2), if so, does that knee fall inside the band you're fitting alpha over?
%[text] ## The five real datasets (per-trial results, PowerLawSimulationPreReg Finding \#99)
%[text] Per-trial fitting (each trial's residual fitted independently, then summarised) is canonical here, not a single pooled fit across all trials -- pooling PSDs across trials with different periodic-orbit structure was found to distort the result. These are aggregate summary statistics from the original analysis; raw trial data is not shipped with this toolbox (see README "Origin").
dataset = ["Cook CTRL"; "Cook ASD"; "Hickman PLAC"; "Hickman HALO"; "Zarandi"];
nTrials = [94; 102; 359; 338; 133];
alphaMedian = [3.765; 4.132; 4.509; 4.570; 2.656];
deltaAICMedian = [2.00; 2.00; -5.66; -4.41; 2.00];
pctDeltaAICSupportsFixed = [85; 64; 39; 43; 100]; % % of trials with deltaAIC >= -2
kneeHzMedian = [0.000; 0.000; 0.289; 0.268; 0.000];
pctKneeInIRASABand = [0; 2; 4; 6; 0]; % % of trials with knee inside [1, fs/2.2] Hz

results = table(dataset, nTrials, alphaMedian, deltaAICMedian, ...
    pctDeltaAICSupportsFixed, kneeHzMedian, pctKneeInIRASABand, ...
    VariableNames=["Dataset", "N_trials", "alpha_median", "deltaAIC_median", ...
    "pct_fixed_adequate", "knee_Hz_median", "pct_knee_in_band"]);
disp(results)
%%
%[text] ## Interpretation of the real result: three flavours of "no problem here"
%[text] **Cook CTRL, Cook ASD, Zarandi -- pure power law.** Median deltaAIC = +2.00 exactly -- the full AIC penalty for the extra knee-model parameter with zero RSS improvement. `checkSpectralKnee`'s optimiser collapses the fitted knee to 0 Hz trial after trial. These datasets are empirically scale-free from the wide-band floor upward.
%[text] **Hickman PLAC / HALO -- genuine sub-floor knee.** Both conditions show a real spectral knee around 0.27-0.29 Hz -- but that's well below the standard IRASA fit band floor of 1 Hz, so only 4-6% of trials show a knee actually inside the fitted band. Knee exists, but it's not biasing anything.
%[text] Notably absent from all five real datasets: a knee that actually falls **inside** the fit band, which would be the case that genuinely biases alpha. The live demonstration below constructs exactly that case, along with the two patterns above, against known ground truth.
%[text] ## Building synthetic data with an EXACT known knee
%[text] The real data can tell you a diagnostic's **qualitative** verdict was right, but not whether its **recovered numbers** (kneeHz, alpha) were close to a true underlying value -- the true value isn't observable in real data. This section builds surrogate residuals with an analytically exact Lorentzian spectrum, so the recovered numbers can be checked against a genuine ground truth.
%[text] The construction: compute the target Lorentzian magnitude spectrum `(kneeHz^2 + f^2)^(-alpha/4)` (square root of the Lorentzian PSD used by `checkSpectralKnee`'s own model) directly at every FFT bin frequency for a signal of length N at sampling rate fs, then hand that magnitude array to the toolbox's bundled `fftnoise` (Grinsted, 2011) to randomise its phase while preserving the exact magnitude -- precisely what `fftnoise` is for, just fed an analytic target spectrum instead of an empirical one. This example calls fftnoise directly, so it needs `external/fftnoise/` on the path in addition to `toolbox/` -- a plain `addpath('.../toolbox')` does not reach subfolders, so that's added explicitly below rather than assumed.
exampleDir = fileparts(mfilename("fullpath"));
addpath(fullfile(exampleDir, "..", "external", "fftnoise"));

fs = 133;
N = 4000;

lorentzianTemplate = @(alpha, kneeHz, sigma) buildLorentzianTemplate(N, fs, alpha, kneeHz, sigma);
%%
%[text] ## Case 1: sub-floor knee, like the real Hickman PLAC/HALO result
%[text] True alpha=4.5, true knee=0.28 Hz (below the standard 1 Hz fit floor)
rng(1);
subFloorSignal = lorentzianTemplate(4.5, 0.28, 3.0);
resultSubFloor = checkSpectralKnee(subFloorSignal, fs, WideBandLow=0.05);
fprintf("Sub-floor knee: verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f (true 0.28)  alphaWide=%.3f (true 4.50)  kneeInBand=%d\n", ...
    resultSubFloor.verdict, resultSubFloor.deltaAIC, resultSubFloor.kneeHz, resultSubFloor.alphaWide, resultSubFloor.kneeInBand);
%%
%[text] Matches the qualitative Hickman pattern: strong knee support, but `kneeInBand=false` -- and here, unlike the real data, the recovered `kneeHz` and `alphaWide` can be checked directly against the true values used to generate the signal. Recovery is close on both.
%[text] ## Case 2: in-band knee -- the case none of the five real datasets showed
%[text] True alpha=4.5, true knee=3 Hz (genuinely inside the default \[1,fHigh\] fit band). This is the actual bias-risk scenario the diagnostic exists to catch.
rng(2);
inBandSignal = lorentzianTemplate(4.5, 3.0, 3.0);
resultInBand = checkSpectralKnee(inBandSignal, fs);
fprintf("In-band knee:   verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f (true 3.00)  alphaWide=%.3f (true 4.50)  kneeInBand=%d\n", ...
    resultInBand.verdict, resultInBand.deltaAIC, resultInBand.kneeHz, resultInBand.alphaWide, resultInBand.kneeInBand);
fprintf("  alphaStd (naive fixed-band fit, straddling the knee): %.3f -- biased away from true 4.50\n", resultInBand.alphaStd);
%%
%[text] Here `kneeInBand=true`, correctly flagging the risk. And critically, `alphaStd` -- the naive power-law fit over the standard band, which straddles the knee -- comes out visibly below the true alpha=4.5, while `alphaWide` (the wide-band knee-aware fit) stays close to it. This is the actual mechanism `checkSpectralKnee` exists to catch: not just "a knee exists" but "a knee is biasing your alpha estimate", shown here with a quantified, known-ground-truth bias, not just an assertion.
%[text] ## Case 3: no knee at all, for contrast
rng(3);
flatSignal = noiseXu(N, 3.8, 3.0, fs);
resultFlat = checkSpectralKnee(flatSignal, fs);
fprintf("No knee:        verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f\n", ...
    resultFlat.verdict, resultFlat.deltaAIC, resultFlat.kneeHz);
%%
%[text] Matches the Cook/Zarandi pattern: deltaAIC near +2 (the pure parameter- penalty signature), fitted knee collapsed to (near) zero.
%[text] ## Takeaway for using this function
%[text] - A "none" verdict means treat the fixed power-law alpha estimate as reliable for that residual.
%[text] - A "marginal" or "strong" verdict with `kneeInBand=false` means a knee exists but isn't biasing your alpha fit -- informative, not alarming (the real Hickman result, and Case 1 above).
%[text] - A "marginal" or "strong" verdict with `kneeInBand=true` is the case that actually warrants caution -- Case 2 above quantifies exactly how much bias that can mean.
%[text] - Run this per-trial, not on a single pooled PSD across many trials -- pooling across trials with different structure can produce spurious or boundary-artefact knee fits (see the original Finding \#99 writeup's Section 7.1 for a worked example of exactly this failure mode). \
%[text] ## Cross-references
%[text] - `checkSpectralKnee.m` docstring
%[text] - Preston, M., Smith, S. & Voytek, B. (2026). Potential mechanisms and functional significance of aperiodic neural activity. Nature Human Behaviour.
%[text] - Original analysis: PowerLawSimulationPreReg  `docs/README_PrestonKneeAnalysis_v001.md`, Finding \#99
%[text] - `external/fftnoise/fftnoise.m`, Grinsted (2011) -- the phase randomiser used to build the synthetic templates above \
function x = buildLorentzianTemplate(N, fs, alpha, kneeHz, sigma)
%BUILDLORENTZIANTEMPLATE Synthetic residual with an exact known Lorentzian spectrum.
%   Example-local helper, not a toolbox function -- scoped to this demo.
%   Builds the target magnitude spectrum (kneeHz^2+f^2)^(-alpha/4)
%   directly at each FFT bin frequency (MATLAB's standard bin ordering:
%   DC, ascending positive frequencies, then negative), then hands it to
%   the bundled fftnoise (Grinsted, 2011) to randomise phase while
%   preserving that exact magnitude.
freqIdx = (0:N - 1)';
freqHz = freqIdx * fs / N;
freqHz(freqHz > fs / 2) = freqHz(freqHz > fs / 2) - fs;
magSpectrum = (kneeHz ^ 2 + freqHz .^ 2) .^ (-alpha / 4);
x = fftnoise(magSpectrum);
x = x - mean(x);
x = x * (sigma / std(x));
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
