%[text] # Agreement Analysis Worked Example
%[text] A realistic method-comparison scenario: you have two ways of measuring the same quantity (two instruments, two raters, an old method and a new one) and want to know whether they can be used interchangeably. This walks through three scenarios that look similar on a scatter plot but tell very different stories once you compute CCC.
%[text] ## The scenario
%[text] Two wearable heart-rate monitors, worn simultaneously by 40 participants during a standard exercise protocol, each reporting a single summary heart rate (bpm) per participant. You want to know: can Device B replace the (more expensive, more validated) Device A?
n = 40;
rng(10);
deviceA = 70 + 60 * rand(n, 1);  % realistic HR range across a mixed sample
%%
%[text] ## Scenario 1: genuinely good agreement
%[text] Device B has honest measurement noise but no systematic bias.
deviceB_good = deviceA + randn(n, 1) * 3;

result1 = linCCC(deviceA, deviceB_good);
fprintf("Scenario 1 (good agreement): CCC=%.3f [%.3f, %.3f], r=%.3f, Cb=%.3f\n", ...
    result1.ccc, result1.LowerCI, result1.UpperCI, result1.r, result1.Cb);
%%
%[text] ## Scenario 2: systematic bias, hidden by a high correlation
%[text] Device B reads consistently ~8 bpm high -- a calibration offset. Pearson r stays high because the two still move together perfectly; only CCC penalises the offset.
deviceB_biased = deviceA + 8 + randn(n, 1) * 3;

result2 = linCCC(deviceA, deviceB_biased);
fprintf("Scenario 2 (systematic bias): CCC=%.3f [%.3f, %.3f], r=%.3f, Cb=%.3f, u=%.3f\n", ...
    result2.ccc, result2.LowerCI, result2.UpperCI, result2.r, result2.Cb, result2.u);
%%
%[text] Notice `r` barely moved from Scenario 1, but `ccc` dropped and `Cb` shows why: the location shift `u` is now nonzero and substantial. A report that only quoted Pearson r here would miss a real, fixable calibration problem.
%[text] ## Scenario 3: proportional bias (scale error)
%[text] Device B under-reads at low heart rates and over-reads at high heart rates -- a gain/calibration slope error rather than a fixed offset.
deviceB_scaled = 70 + (deviceA - 70) * 1.25 + randn(n, 1) * 3;

result3 = linCCC(deviceA, deviceB_scaled);
fprintf("Scenario 3 (proportional bias): CCC=%.3f [%.3f, %.3f], r=%.3f, Cb=%.3f, v=%.3f\n", ...
    result3.ccc, result3.LowerCI, result3.UpperCI, result3.r, result3.Cb, result3.v);
%%
%[text] Here the scale shift `v` (ratio of standard deviations) is the component driving the Cb penalty, not `u`. Distinguishing these two failure modes matters in practice: a fixed offset (Scenario 2) is often correctable with a simple additive calibration; a scale error (Scenario 3) usually needs a multiplicative correction, or points to a more fundamental measurement problem across the device's range.
%[text] ## Visual comparison
figure;
tiledlayout(1, 3);
scenarios = {deviceB_good, deviceB_biased, deviceB_scaled};
titles = ["Good agreement", "Systematic bias", "Proportional bias"];
results = {result1, result2, result3};
for i = 1:3
    nexttile;
    scatter(deviceA, scenarios{i}, 20, "filled");
    hold on;
    lims = [min(deviceA) - 5, max(deviceA) + 5];
    plot(lims, lims, "k--", DisplayName="Identity line");
    hold off;
    xlim(lims); ylim(lims);
    axis square;
    xlabel("Device A (bpm)"); ylabel("Device B (bpm)");
    title(sprintf("%s\nCCC=%.3f, r=%.3f", titles(i), results{i}.ccc, results{i}.r));
end
%%
%[text] All three scatter plots look like a strong linear relationship at a glance -- the visual difference between them is subtle. The CCC values make the distinction explicit and quantifiable in a way eyeballing a scatter plot does not reliably do.
%[text] ## Reporting the result
%[text] A defensible write-up states the CCC, its confidence interval, the sample size, and the McBride (2005) interpretation band -- not just a bare number:
fprintf("\nExample reportable sentence:\n");
fprintf("""Agreement between Device A and Device B was %s (CCC = %.3f, 95%% CI [%.3f, %.3f], n = %d).""\n", ...
    lower(string(interpretCCC(result2.ccc))), result2.ccc, result2.LowerCI, result2.UpperCI, result2.n);
%%
%[text] ## Cross-references
%[text] - `linCCC.m` docstring -- full formula and citation list
%[text] - `doc/GettingStarted.m` -- basic usage and output interpretation
%[text] - `doc/WhyValidatedAgainstR.m` -- validation against R's DescTools::CCC \
function label = interpretCCC(ccc)
%INTERPRETCCC McBride (2005) interpretation bands, for the reporting sentence above.
if ccc > 0.99
    label = "Almost perfect agreement";
elseif ccc > 0.95
    label = "Substantial agreement";
elseif ccc > 0.90
    label = "Moderate agreement";
else
    label = "Poor agreement";
end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
