function out = analyseR29PostHoc_v001()
%ANALYSER29POSTHOC_V001  POST HOC description after R29 (spec 10.54; no pre-stated rule): beta_obs slope by model variant (V1 as published, V2 true f0, V3 no noise covariates) for the four arms,
%   and the median estimated alpha of the empirical (v015) trials against the synthetic arms, to see whether the structure effect passes through the lsa/al covariates and whether the level scaling changed colour.
%   Writes results/tempoSim/analyseR29PostHoc_v001.mat/.txt   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  ROOT = fileparts(srcDir);  addpath(srcDir);  simDir = fullfile(ROOT,'results','tempoSim');
txt = fullfile(simDir,'analyseR29PostHoc_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR29PostHoc_v001  %s  (POST HOC, no pre-stated rule)\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
ARMS = ["AXbF" "AXpF" "AXvF" "AXwF"];  rows = table();  al = table();
CASES = {"Fraser","SG-IRLS"; "Fraser","BWFD-OLS"; "Hickman_PLAC","SG-IRLS"};
for ds = ["Fraser" "Hickman_PLAC"]
    for arm = ARMS
        T = analyseTempoSlopeRecovery_v002(arm, Datasets=ds, Tag="PerDS_"+ds).slopes;
        for c = find(string(CASES(:,1))' == ds)
            p = string(CASES{c,2});
            for v = ["V1" "V2" "V3"]
                r = T(T.source=="synthetic" & T.variant==v & T.pipeline==p & T.outcome=="betaObs", :);
                if height(r) ~= 1, error('R29ph:pick','%s %s %s %s: %d rows',ds,arm,p,v,height(r)); end
                rows = [rows; table(ds, p, arm, v, r.b, r.SE, 'VariableNames', ["dataset","pipeline","arm","variant","b","SE"])]; %#ok<AGROW>
            end
        end
        R = load(fullfile(simDir,"loopClosureResults_SYN-"+arm+"-"+ds+"_all_shaped_xu_syn001.mat"),'results').results;
        al = [al; table(ds, arm, median([R.alphaIRA],'omitnan'), 'VariableNames', ["dataset","arm","medianAlpha"])]; %#ok<AGROW>
    end
    E = load(fullfile(srcDir,"loopClosureResults_"+ds+"_all_shaped_xu_v015.mat"),'results').results;
    al = [al; table(ds, "EMPIRICAL", median([E.alphaIRA],'omitnan'), 'VariableNames', ["dataset","arm","medianAlpha"])]; %#ok<AGROW>
end
W = unstack(rows(:,["dataset","pipeline","arm","variant","b"]), "b", "arm");
fprintf('beta_obs slope per doubling by variant and arm:\n'); disp(W); fprintf('median estimated alpha, synthetic arms against the empirical trials:\n'); disp(al);
out = struct('slopes', rows, 'wide', W, 'alpha', al);  save(fullfile(simDir,'analyseR29PostHoc_v001.mat'), '-struct', 'out');
fprintf('=== DONE ===\n');
end
