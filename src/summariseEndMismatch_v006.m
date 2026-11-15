function T = summariseEndMismatch_v006()
%SUMMARISEENDMISMATCH_V006 Collate the endMismatchSensitivity_v006 summaries over the five arms (alpha-selected kinds).
%   T = summariseEndMismatch_v006() reads results/endMismatchSensitivity/endMismatch_<arm>_summary_v006.csv and the
%   per-trial CSVs, and writes results/endMismatchSensitivity/endMismatch_collated_v006.csv with one row per arm x kind:
%   median Lfrac (p10-p90 over trial-axes), the BWFD velocity and acceleration noise ratios against the scoring reference
%   (segment for e2e kinds, whole residual for ext kinds) and against the whole residual (Full), each as "maj ; min"
%   medians, the p10-p90 over trial-axes, the share of surrogates with no accepted draw, and FAIL counts.
%   Also prints the fairShare segment-length distribution per arm and the share of trial-axes where fairShare chose a
%   segment shorter than 75% / 90% of the trial.
here = fileparts(mfilename('fullpath'));  d = fullfile(here, '..', 'results', 'endMismatchSensitivity');
arms = ["FraserN48", "HickmanPLAC", "HickmanHALO", "CookCTRL", "CookASD"];
kinds = ["e2e75Sel", "e2e90Sel", "e2eFairSel", "extWrapFairSel", "extSpecFairSel", "extSpec75Sel"];
T = table();  L = table();
for a = arms
    R = readtable(fullfile(d, "endMismatch_" + a + "_v006.csv"), 'TextType', 'string');
    nFail = nnz(R.status ~= "ok");
    for k = kinds
        Rk = R(R.kind == k & R.status == "ok", :);  mj = Rk.axis == "maj";  mn = Rk.axis == "min";
        f = @(v, m) median(v(m), 'omitnan');  pr = @(v) sprintf('%.2f-%.2f', prctile(v, 10), prctile(v, 90));
        T = [T; table(a, k, height(Rk) / 2, f(Rk.Lfrac, true(height(Rk), 1)), string(pr(Rk.Lfrac)), ...
            string(sprintf('%.2f ; %.2f', f(Rk.velRatioBW, mj), f(Rk.velRatioBW, mn))), string(pr(Rk.velRatioBW)), ...
            string(sprintf('%.2f ; %.2f', f(Rk.accRatioBW, mj), f(Rk.accRatioBW, mn))), string(pr(Rk.accRatioBW)), ...
            string(sprintf('%.2f ; %.2f', f(Rk.velRatioBWFull, mj), f(Rk.velRatioBWFull, mn))), ...
            string(sprintf('%.2f ; %.2f', f(Rk.accRatioBWFull, mj), f(Rk.accRatioBWFull, mn))), ...
            mean(Rk.shareNoAccept, 'omitnan'), nFail, ...
            'VariableNames', {'arm','kind','nTrials','LfracMed','LfracP10P90','velBW','velBW_p10p90','accBW','accBW_p10p90','velBWFull','accBWFull','shareNoAccept','nFailRowsArm'})]; %#ok<AGROW>
    end
    Rf = R(R.kind == "e2eFairSel" & R.status == "ok", :);
    L = [L; table(a, height(Rf), median(Rf.Lfrac), prctile(Rf.Lfrac, 10), prctile(Rf.Lfrac, 90), min(Rf.Lfrac), mean(Rf.Lfrac < 0.75), mean(Rf.Lfrac < 0.90), ...
        'VariableNames', {'arm','nTrialAxes','LfracMed','p10','p90','min','shareBelow75','shareBelow90'})]; %#ok<AGROW>
end
writetable(T, fullfile(d, "endMismatch_collated_v006.csv"));  writetable(L, fullfile(d, "endMismatch_fairShareLength_v006.csv"));
disp(T);  disp(L);
end
