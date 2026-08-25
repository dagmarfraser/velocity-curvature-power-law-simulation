%[text] # Comparing Noise Generation Options
%[text] This toolbox offers several ways to get 1/f^alpha-ish noise, each suited to a different situation. This example generates the same target alpha via each route and compares the recovered alpha and spectral shape, to help you choose.
%[text] ## The options
%[text] - `noiseXu` -- parametric generator. You specify alpha directly; it synthesises a signal to match. Use this when you want noise with a  **known** ground-truth alpha (simulation, testing, stress-testing a pipeline).
%[text] - `shapeNoise` with `PhaseSource="xu"` (its default) -- takes an  **empirical** residual's amplitude spectrum and imposes it onto a fresh  `noiseXu` phase carrier. Use this when you have real data and want a surrogate that matches its spectral shape without reusing its exact realisation.
%[text] - `shapeNoise` with other `PhaseSource` choices (`"white"`,  `"pinknoise"`, `"dsp"`, or a custom handle) -- same empirical-magnitude idea, different phase carrier. Useful for testing whether your result is sensitive to the specific phase-generation method, not just the magnitude spectrum.
%[text] - `external/fftnoise/fftnoise.m` -- a **phase-surrogate** generator (Aslak Grinsted, 2011): takes an existing signal's FFT and randomises its phase, keeping the exact magnitude spectrum. Different from  `noiseXu` in kind, not just implementation -- it needs a real FFT to start from, it doesn't synthesise one from a target alpha. Included as an alternative phase source for `shapeNoise`. \
%[text] ## Generating the same target alpha several ways
fs = 100;
N = 3000;
targetAlpha = 2.5;
fLow = 1.0;
fHigh = 25.0;

rng(10);
xuDirect = noiseXu(N, targetAlpha, 3.0, fs);

% For shapeNoise, first need an empirical-ish residual with roughly the
% target spectral shape to draw the magnitude from -- use a noiseXu
% draw for this demo (in real use this would be your actual data).
rng(11);
empiricalStandIn = noiseXu(N, targetAlpha, 3.0, fs);

rng(12);
shapedXu = shapeNoise(empiricalStandIn, fs, targetAlpha, PhaseSource="xu");
rng(13);
shapedWhite = shapeNoise(empiricalStandIn, fs, targetAlpha, PhaseSource="white");
%%
%[text] ## Recovering alpha from each
methodNames = ["noiseXu (direct)"; "shapeNoise (xu carrier)"; "shapeNoise (white carrier)"];
signals = {xuDirect, shapedXu, shapedWhite};
recoveredAlpha = zeros(numel(signals), 1);
for i = 1:numel(signals)
    recoveredAlpha(i) = estimateIRASA(signals{i}, fs, fLow, fHigh);
end

comparison = table(methodNames, recoveredAlpha, ...
    VariableNames=["Method", "alpha_recovered"]);
disp(comparison)
fprintf("Target alpha: %.2f\n", targetAlpha);
%%
%[text] ## Spectral comparison
figure;
tiledlayout(1, 2);

nexttile;
[pxxXu, fXu] = pmtm(xuDirect, 4, [], fs);
[pxxShapedXu, fShapedXu] = pmtm(shapedXu, 4, [], fs);
[pxxShapedWhite, fShapedWhite] = pmtm(shapedWhite, 4, [], fs);
loglog(fXu, pxxXu, DisplayName="noiseXu (direct)");
hold on;
loglog(fShapedXu, pxxShapedXu, DisplayName="shapeNoise (xu carrier)");
loglog(fShapedWhite, pxxShapedWhite, DisplayName="shapeNoise (white carrier)");
hold off;
xlabel("Frequency (Hz)");
ylabel("PSD");
title("Power spectra (log-log)");
legend(Location="southwest");
grid on;

nexttile;
bar(categorical(methodNames), recoveredAlpha);
yline(targetAlpha, "r--", "target alpha", LineWidth=1.5);
ylabel("Recovered alpha");
title("Alpha recovery by method");
xtickangle(20);
%%
%[text] All three methods should cluster near the target alpha -- the point of this comparison isn't that one is "more correct" (they're all valid ways to get 1/f^alpha-shaped output), it's that the **phase** carrier choice barely matters once the magnitude spectrum is fixed, which is itself informative: alpha, as IRASA estimates it, is a magnitude- spectrum property, largely independent of phase structure.
%[text] ## When would the phase carrier choice actually matter?
%[text] If your downstream analysis is sensitive to something **other** than the power spectrum -- e.g. higher-order statistics, waveform shape, specific autocorrelation structure at short lags -- the phase carrier choice can matter even though alpha recovery looks identical here. `shapeNoise`'s `PhaseSource` argument exists precisely so you can test that sensitivity directly, by swapping carriers and re-running your actual downstream analysis, not just an alpha check.
%[text] ## Cross-references
%[text] - `noiseXu.m`, `shapeNoise.m` docstrings
%[text] - `doc/GettingStarted.mlx` for a broader tour of the toolbox \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
