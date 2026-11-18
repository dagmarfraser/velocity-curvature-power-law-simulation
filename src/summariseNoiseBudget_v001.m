function S = summariseNoiseBudget_v001()
%SUMMARISENOISEBUDGET_V001 Summary of noiseBudgetSensitivity_v001 from its saved runner outputs (2026-11-15).
%   The 12:35 run hung at Cook ASD (MATLAB restarted by Dagmar); Cook ASD was rerun separately. This reads every
%   loopClosureResults_<arm>_per_subject_endMatched_bud1[_ns<scale>]_nbs1_v027.mat in src/ and reports, per arm and scale,
%   the per-trial shift of betaGenStarMed (and of SG-IRLS) against scale 1: median, 90th percentile of |shift|, share beyond
%   0.003 (a tenth of the MDC) and beyond the MDC (0.03). Descriptive (Dagmar, 11:05: the budget is chosen without reference to
%   beta_gen*; this only reports its effect). Writes results/templateComparison/noiseBudgetSensitivity_summary_v002.csv.
here = fileparts(mfilename('fullpath'));
arms = ["FraserN48", "Hickman_PLAC", "Hickman_HALO", "Cook_CTRL", "Cook_ASD"];  scales = [1 1.01 1.02 1.05];  S = table();
for a = arms
    R = cell(1, 4);
    for k = 1:4
        tag = "";  if scales(k) ~= 1, tag = "_ns" + string(scales(k)); end
        f = fullfile(here, "loopClosureResults_" + a + "_per_subject_endMatched_bud1" + tag + "_nbs1_v027.mat");
        if isfile(f), R{k} = load(f, 'results').results; end
    end
    if isempty(R{1}), continue, end
    ids = string({R{1}.trialID});  b0 = [R{1}.betaGenStarMed];  s0 = arrayfun(@(r) r.betaGenStar(6), R{1});
    for k = 2:4
        if isempty(R{k}), continue, end
        [~, i0, i1] = intersect(ids, string({R{k}.trialID}));
        d = [R{k}(i1).betaGenStarMed] - b0(i0);  d = d(isfinite(d));
        ds = arrayfun(@(r) r.betaGenStar(6), R{k}(i1)) - s0(i0);  ds = ds(isfinite(ds));
        S = [S; table(a, scales(k), numel(d), median(d), prctile(abs(d), 90), mean(abs(d) > 0.003), mean(abs(d) > 0.03), median(ds), prctile(abs(ds), 90), ...
            'VariableNames', {'arm', 'scale', 'n', 'medShift', 'p90AbsShift', 'shareOver0003', 'shareOverMDC', 'medShiftSGIRLS', 'p90AbsShiftSGIRLS'})]; %#ok<AGROW>
    end
end
writetable(S, fullfile(here, '..', 'results', 'templateComparison', 'noiseBudgetSensitivity_summary_v002.csv'));  disp(S);
end
