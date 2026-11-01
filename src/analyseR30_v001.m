function out = analyseR30_v001(DS)
%% analyseR30_v001  R30 (spec 10.55): which part of the residual's structure carries the pipeline-made beta_obs tempo slope?
% Needs the saved runs for AXcF, AXaF, AXsF, AXnF (src/runR30HPC_v001.m) for Fraser and Hickman PLAC, plus the saved AXr baselines (as analyseR23_v001).
% R30a: slopes from analyseTempoSlopeRecovery_v002 (per dataset, Tag PerDS_<dataset>), variants V1 and V3 separately; contrast d(arm) = s(arm) - s(AXr),
%   SE = sqrt(sum of arm variances), |z| >= 2 = interval excludes zero, SEM threshold 0.0108 per doubling. Labels exactly as 10.55 (not exclusive).
% R30b: within-session corr(true log2 f0, log2(sigma/a)) and corr(.., alpha) per arm; median alpha per arm vs the AXr run (if its run file is present) and vs the empirical v015 median.
% Writes results/tempoSim/analyseR30_v001.mat/.txt   Fraser, D.S. (2026) v001
srcDir = fileparts(mfilename('fullpath'));  ROOT = fileparts(srcDir);  addpath(srcDir);
simDir = fullfile(ROOT, 'results', 'tempoSim');
txt = fullfile(simDir, 'analyseR30_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR30_v001  %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
if nargin < 1, DS = ["Fraser" "Hickman_PLAC"]; end;  ARMS = ["AXcF" "AXaF" "AXsF" "AXnF"];  VAR = ["V1" "V3"];  SEM = 0.0108;
CASES = {"Fraser","SG-IRLS"; "Fraser","BWFD-OLS"; "Hickman_PLAC","SG-IRLS"};
rows = table();  corrRows = table();
for ds = DS
    Fit = struct();
    for arm = ARMS
        Fit.(arm) = analyseTempoSlopeRecovery_v002(arm, Datasets=ds, Tag="PerDS_"+ds).slopes;
    end
    if ds == "Fraser", bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrF.mat"); else, bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrPerDS_"+ds+".mat"); end
    if ~isfile(bf), error('R30:baseline', 'Missing AXr baseline fit: %s', bf); end
    Fit.AXr = load(bf,'out').out.slopes;
    for c = find(string(CASES(:,1))' == ds)
        p = string(CASES{c,2});
        for v = VAR
            g = @(arm) pick_local(Fit.(arm), p, v, ds, arm);
            sR = g("AXr");  sC = g("AXcF");  sA = g("AXaF");  sS = g("AXsF");  sN = g("AXnF");
            dz = @(a, b) deal(a.b - b.b, (a.b - b.b)/sqrt(a.SE^2 + b.SE^2));
            [dC, zC] = dz(sC, sR);  [dS, zS] = dz(sS, sR);  [dN, zN] = dz(sN, sR);  [dA, zA] = dz(sA, sR);
            [~, zAC] = dz(sA, sC);  [~, zRA] = dz(sR, sA);
            lab = strings(0,1);
            if abs(dC) < SEM && abs(zC) < 2, lab(end+1) = "COLOUR SUFFICIENT"; end %#ok<AGROW>
            if abs(zS) >= 2 || abs(zAC) >= 2, lab(end+1) = "CROSS-AXIS MATTERS"; end %#ok<AGROW>
            if abs(zRA) >= 2, lab(end+1) = "PHASE RESIDUE MATTERS"; end %#ok<AGROW>
            if abs(zN) >= 2, lab(end+1) = "COLOUR MATTERS"; end %#ok<AGROW>
            if isempty(lab), lab = "NOT RESOLVED"; end
            rows = [rows; table(ds, p, v, sR.b, sC.b, sA.b, sS.b, sN.b, dC, zC, dA, zA, dS, zS, dN, zN, zAC, zRA, strjoin(lab," + "), ...
                'VariableNames', ["dataset","pipeline","variant","sXr","sXc","sXa","sXs","sXn","dXc","zXc","dXa","zXa","dXs","zXs","dXn","zXn","zXa_Xc","zXr_Xa","label"])]; %#ok<AGROW>
        end
    end
    % R30b
    emp = load(fullfile(srcDir, "loopClosureResults_" + ds + "_all_shaped_xu_v015.mat"), 'results').results;
    medEmp = median(arrayfun(@(x) x.alphaIRA, emp), 'omitnan');   % same field the synthetic runs report
    medAXr = NaN;  fr = fullfile(simDir, "loopClosureResults_SYN-AXrF-" + ds + "_all_shaped_xu_syn001.mat");
    if ~isfile(fr), fr = fullfile(simDir, "loopClosureResults_SYN-AXr-" + ds + "_all_shaped_xu_syn001.mat"); end
    if isfile(fr)
        medAXr = median(arrayfun(@(x) x.alphaIRA, load(fr,'results').results), 'omitnan');
    else
        warning('R30:noAXrRun', '%s: AXr run file not on this machine; the alpha comparator is the empirical median only.', ds);
    end
    for arm = ARMS
        f = fullfile(simDir, "loopClosureResults_SYN-" + arm + "-" + ds + "_all_shaped_xu_syn001.mat");
        if ~isfile(f), error('R30:input', 'Missing: %s', f); end
        r = load(f, 'results').results;
        t = table(string({r.subjectID})', arrayfun(@(x) x.f0, r)', arrayfun(@(x) x.sigmaMM, r)', arrayfun(@(x) x.a_mm, r)', arrayfun(@(x) x.alphaIRA, r)', ...
            'VariableNames', ["session","f0","sigma","a","alpha"]);
        lf = log2(t.f0);  lsa = log2(t.sigma ./ t.a);  ok = isfinite(lf) & isfinite(lsa) & isfinite(t.alpha) & t.f0 > 0;
        t = t(ok,:);  lf = lf(ok);  lsa = lsa(ok);  gi = findgroups(t.session);  cen = @(x) x - cenMean_local(x, gi);
        S = load(fullfile(simDir, "tempoSim_" + arm + "_" + ds + ".mat"), 'dDiag');
        ma = median(t.alpha);
        corrRows = [corrRows; table(ds, arm, corr(cen(lf), cen(lsa)), corr(cen(lf), cen(t.alpha)), ma, ma - medAXr, ma - medEmp, S.dDiag.corrOwn, S.dDiag.corrReal, ...
            'VariableNames', ["dataset","arm","corrLfwLsa","corrLfwAlpha","medianAlpha","dAlphaVsAXr","dAlphaVsEmp","corrLevelOwn","corrLevelReal"])]; %#ok<AGROW>
    end
end
corrRows.colourPreserved = abs(corrRows.dAlphaVsAXr) <= 0.05;
corrRows.colourPreserved(corrRows.arm == "AXnF") = false;   % expected NOT preserved (white control); reported, not tested
fprintf('\nR30a (beta_obs slope per doubling; d = s(arm) - s(AXr); z by model SE; labels from 10.55):\n'); disp(rows);
fprintf('R30b (within-session correlations; alpha medians; dAlphaVsAXr is NaN where the AXr run file is absent; colourPreserved uses AXr only):\n'); disp(corrRows);
out = struct('R30a', rows, 'R30b', corrRows);
save(fullfile(simDir, 'analyseR30_v001.mat'), '-struct', 'out');
fprintf('=== DONE ===\n');
end

function r = pick_local(T, p, v, ds, arm)
r = T(T.source == "synthetic" & T.variant == v & T.pipeline == p & T.outcome == "betaObs", :);
if height(r) ~= 1, error('R30:pick', '%s %s %s %s: %d rows', ds, p, v, arm, height(r)); end
end

function m = cenMean_local(x, g)
mm = splitapply(@mean, x, g);  m = mm(g);
end
