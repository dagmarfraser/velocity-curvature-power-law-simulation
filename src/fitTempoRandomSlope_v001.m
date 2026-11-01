function out = fitTempoRandomSlope_v001()
%FITTEMPORANDOMSLOPE_V001  Published tempo decomposition (random intercept) against a random-slope refit, for lfw.
%   out = fitTempoRandomSlope_v001()          from the project root or src; minutes on a laptop
%   Spec and pre-stated reading: docs/SPEC_T2_TempoSlopeRecovery_v001.md section 10.15 (written before the fit).
%   Data as checkTempoDecompositionDomain_v001 (Dhieb out; lfb re-centred; lsa, al re-z-scored). Models, per pipeline
%   (SG-IRLS, BWFD-OLS), variant (V1 as published; V3 without lsa, al) and outcome (betaObs, betaGen):
%     RI  outcome ~ 1 + lfw + lfb [+ lsa + al] + datasetC + (1|sessionC)          published
%     RS  outcome ~ 1 + lfw + lfb [+ lsa + al] + datasetC + (1 + lfw|sessionC)    random tempo slope per session (= per person)
%   Reports lfw: estimate, SE, t, Satterthwaite df, p, 95% CI; the random-slope SD; the LR test RS vs RI; fitlme warnings verbatim.
%   Writes results/tempoSim/fitTempoRandomSlope_v001.mat and .txt.   Fraser, D.S. (2026)  v001
root = fileparts(fileparts(mfilename('fullpath')));  simDir = fullfile(root,'results','tempoSim');
inMat = fullfile(root,'results','tempoDecomposition_v002.mat');
if ~isfile(inMat), error('rsTempo:input','%s', "FAILED PATH: " + inMat); end
if ~isfolder(simDir), mkdir(simDir); end
diary(fullfile(simDir,'fitTempoRandomSlope_v001.txt')); cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('fitTempoRandomSlope_v001 | %s | root %s\n', string(datetime('now')), root);
S = load(inMat,'T','PIPES');  PIPES = string(S.PIPES);  T = S.T(S.T.dataset ~= "Dhieb", :);
T.dataset = string(T.dataset);  T.session = string(T.session);  T.pipeline = string(T.pipeline);
z = @(x) (x - mean(x,'omitnan')) ./ std(x,'omitnan');                       % checkTempoDecompositionDomain_v001 L28-34
G = groupsummary(T(T.pipeline == PIPES(1),:), 'session', 'mean', 'lf');
T.lfb = T.lfMean - mean(G.mean_lf);  T.lsa = z(log2(T.sigma ./ T.a));  T.al = z(T.alpha);
T.datasetC = categorical(T.dataset);  T.sessionC = categorical(T.session);

R = table();
for p = PIPES
    V = T(T.pipeline == p,:);
    for variant = ["V1","V3"]
        cov = " + lsa + al";  if variant == "V3", cov = ""; end
        for oc = ["betaObs","betaGen"]
            W = V(isfinite(V.(oc)),:);
            [mi, wi] = fit_local(W, oc + " ~ 1 + lfw + lfb" + cov + " + datasetC + (1|sessionC)");
            [ms, ws] = fit_local(W, oc + " ~ 1 + lfw + lfb" + cov + " + datasetC + (1 + lfw|sessionC)");
            lr = compare(mi, ms, 'CheckNesting', true);
            [~, ~, st] = covarianceParameters(ms);  s1 = st{1};                                         % sessionC group
            q = string(s1.Type) == "std" & string(s1.Name1) == "lfw";
            if nnz(q) ~= 1, error('rsTempo:sd','Random-slope SD not found in the covariance table.'); end
            sdSlope = s1.Estimate(q);
            for tag = ["RI","RS"]
                m = mi;  w = wi;  if tag == "RS", m = ms;  w = ws; end
                [~, ~, fe] = fixedEffects(m, 'DFMethod', 'satterthwaite');  k = string(fe.Name) == "lfw";
                R = [R; table(p, variant, oc, tag, height(W), fe.Estimate(k), fe.SE(k), fe.tStat(k), fe.DF(k), fe.pValue(k), fe.Lower(k), fe.Upper(k), ...
                    ternary_local(tag == "RS", sdSlope, NaN), ternary_local(tag == "RS", lr.pValue(2), NaN), string(w), ...
                    'VariableNames', ["pipeline","variant","outcome","model","n","b","SE","t","df","p","CIlo","CIhi","sdSlope","pLRT","warning"])]; %#ok<AGROW>
            end
        end
    end
end
disp(R(:, setdiff(R.Properties.VariableNames, "warning", 'stable')))
if any(R.warning ~= ""), fprintf('\nFITLME WARNINGS (verbatim, never suppressed):\n'); disp(R(R.warning ~= "", ["pipeline","variant","outcome","model","warning"])); end

% regression anchor: RI estimates are the published slopes
A = load(fullfile(root,'results','checkTempoDecompositionDomain_v001.mat'),'D1');  D = A.D1{2};
for p = PIPES
    for pair = [["betaObs","beta_obs"]; ["betaGen","beta_gen*"]]'
        ref = D.b(D.pipeline == p & D.outcome == pair(2) & D.term == "lfw");
        got = R.b(R.pipeline == p & R.variant == "V1" & R.outcome == pair(1) & R.model == "RI");
        fprintf('anchor %-8s %-9s published %+.4f   RI here %+.4f\n', p, pair(2), ref, got);
        if abs(round(got,4) - ref) > 1e-4, error('rsTempo:anchor','RI %s %s does not reproduce the published slope.', p, pair(2)); end
    end
end
sg = R(R.pipeline == PIPES(1) & R.variant == "V1" & R.model == "RS", :);
fprintf('\nREADING (SG-IRLS V1, spec 10.15 R4, R5):\n');
for oc = ["betaObs","betaGen"]
    r = sg(sg.outcome == oc,:);  ok = r.warning == "" && isfinite(r.sdSlope) && r.sdSlope > 1e-6;
    fprintf('  %-8s RS b %+.4f SE %.4f t(%.0f) = %.2f p %.2g CI %+.4f..%+.4f  slope SD %.4f  LRT p %.2g  fit %s\n', oc, r.b, r.SE, r.df, r.t, r.p, r.CIlo, r.CIhi, r.sdSlope, r.pLRT, ternary_local(ok, 'CLEAN (use for t and p)', 'NOT CLEAN (quote bootstrap intervals only)'));
end
r = sg(sg.outcome == "betaGen",:);
fprintf('  R5: RS CI for the beta_gen* slope %s 0\n', ternary_local(r.CIlo > 0 || r.CIhi < 0, 'excludes', 'includes'));
out = struct('R',R,'runDate',datetime('now'),'matlab',version);
save(fullfile(simDir,'fitTempoRandomSlope_v001.mat'),'out');
fprintf('Saved %s\n', fullfile(simDir,'fitTempoRandomSlope_v001.mat'));
end

function [m, w] = fit_local(W, f)
% fit and keep any warning text verbatim
lastwarn('');  m = fitlme(W, f);  w = lastwarn;
end

function v = ternary_local(c, a, b), if c, v = a; else, v = b; end, end
