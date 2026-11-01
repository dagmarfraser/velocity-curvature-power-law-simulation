function out = analyseR29SeamFree_v001(DS)
%% analyseR29SeamFree_v001  R29 robustness (2) of spec 10.56: the R29 question on SEAM-FREE arms Xd (=Xp), Xo (=Xv, own level), Xl (=Xw, donor level), baseline AXr (seam-free own residual).
% Rule (10.56, written before any run): contrasts Delta = s(Xd)-s(AXr), L = s(Xl)-s(AXr) (level decoupled), T = s(Xo)-s(AXr) (structure decoupled); |z| >= 2 by model SE
% (SE = sqrt of summed arm variances); labels as analyseR29_v001 (LEVEL / STRUCTURE / BOTH / PARTIAL / NOT RESOLVED / NOT APPLICABLE), V1 and V3 separately.
% R29b: within-session corr(true log2 f0, log2(sigma/a)) and corr(.., alpha) per arm, and the 10.51 hypothesis check with Xr as comparator.
% Needs the saved runs for AXdF, AXoF, AXlF (src/runOvernightHPC_v001.m) plus the saved AXr baseline fits. Writes results/tempoSim/analyseR29SeamFree_v001.mat/.txt   Fraser, D.S. (2026) v001
srcDir = fileparts(mfilename('fullpath'));  ROOT = fileparts(srcDir);  addpath(srcDir);
simDir = fullfile(ROOT, 'results', 'tempoSim');
txt = fullfile(simDir, 'analyseR29SeamFree_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR29SeamFree_v001  %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
if nargin < 1, DS = ["Fraser" "Hickman_PLAC"]; end;  ARMS = ["AXdF" "AXoF" "AXlF"];  VAR = ["V1" "V3"];
CASES = {"Fraser","SG-IRLS"; "Fraser","BWFD-OLS"; "Hickman_PLAC","SG-IRLS"};
rows = table();  corrRows = table();
for ds = DS
    Fit = struct();
    for arm = ARMS, Fit.(arm) = analyseTempoSlopeRecovery_v002(arm, Datasets=ds, Tag="PerDS_"+ds).slopes; end
    if ds == "Fraser", bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrF.mat"); else, bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrPerDS_"+ds+".mat"); end
    if ~isfile(bf), error('R29sf:baseline', 'Missing AXr baseline fit: %s', bf); end
    Fit.AXr = load(bf,'out').out.slopes;
    for c = find(string(CASES(:,1))' == ds)
        p = string(CASES{c,2});
        for v = VAR
            g = @(arm) pick_local(Fit.(arm), p, v, ds, arm);
            sR = g("AXr"); sD = g("AXdF"); sO = g("AXoF"); sL = g("AXlF");
            cz = @(a, b) deal(a.b - b.b, sqrt(a.SE^2 + b.SE^2));
            [Dl, seD] = cz(sD, sR);  [L, seL] = cz(sL, sR);  [Tt, seT] = cz(sO, sR);
            zD = Dl/seD;  zL = L/seL;  zT = Tt/seT;  l = L/Dl;  t = Tt/Dl;
            if abs(zD) < 2, lab = "NOT APPLICABLE (Delta interval includes zero)";
            elseif abs(zL) >= 2 && abs(zT) < 2 && l >= 0.65 && t <= 0.35, lab = "LEVEL carries it";
            elseif abs(zT) >= 2 && abs(zL) < 2 && t >= 0.65 && l <= 0.35, lab = "STRUCTURE carries it";
            elseif abs(zL) >= 2 && abs(zT) >= 2, lab = "BOTH";
            elseif abs(zL) < 2 && abs(zT) < 2, lab = "NOT RESOLVED";
            else, lab = "PARTIAL (one factor excludes zero but its share is outside the 0.65 / 0.35 rule)"; end
            rows = [rows; table(ds, p, v, sR.b, sD.b, sO.b, sL.b, Dl, zD, L, zL, Tt, zT, l, t, abs(l+t-1) > 0.35, lab, 'VariableNames', ...
                ["dataset","pipeline","variant","sXr","sXd","sXo","sXl","Delta","zDelta","L","zL","T","zT","lShare","tShare","nonAdditive","label"])]; %#ok<AGROW>
        end
    end
    for arm = ARMS
        f = fullfile(simDir, "loopClosureResults_SYN-" + arm + "-" + ds + "_all_shaped_xu_syn001.mat");
        if ~isfile(f), error('R29sf:input', 'Missing: %s', f); end
        r = load(f, 'results').results;
        t = table(string({r.subjectID})', arrayfun(@(x) x.f0, r)', arrayfun(@(x) x.sigmaMM, r)', arrayfun(@(x) x.a_mm, r)', arrayfun(@(x) x.alphaIRA, r)', 'VariableNames', ["session","f0","sigma","a","alpha"]);
        lf = log2(t.f0);  lsa = log2(t.sigma ./ t.a);  ok = isfinite(lf) & isfinite(lsa) & isfinite(t.alpha) & t.f0 > 0;
        t = t(ok,:);  lf = lf(ok);  lsa = lsa(ok);  gi = findgroups(t.session);  cen = @(x) x - cenMean_local(x, gi);
        S = load(fullfile(simDir, "tempoSim_" + arm + "_" + ds + ".mat"), 'dDiag', 'truth');
        corrRows = [corrRows; table(ds, arm, corr(cen(lf), cen(lsa)), corr(cen(lf), cen(t.alpha)), median(t.alpha), S.dDiag.nSingleton, sum(S.truth.wrapped), S.dDiag.corrOwn, S.dDiag.corrReal, ...
            'VariableNames', ["dataset","arm","corrLfwLsa","corrLfwAlpha","medianAlpha","nSkipped","nWrapped","corrLevelOwn","corrLevelReal"])]; %#ok<AGROW>
    end
end
fprintf('\nR29 seam-free (slope per doubling of beta_obs; Delta = Xd-Xr, L = Xl-Xr level decoupled, T = Xo-Xr structure decoupled):\n'); disp(rows);
fprintf('R29b seam-free (nWrapped must be 0):\n'); disp(corrRows);
out = struct('rows', rows, 'corr', corrRows);
save(fullfile(simDir, 'analyseR29SeamFree_v001.mat'), '-struct', 'out');
fprintf('=== DONE ===\n');
end

function r = pick_local(T, p, v, ds, arm)
r = T(T.source == "synthetic" & T.variant == v & T.pipeline == p & T.outcome == "betaObs", :);
if height(r) ~= 1, error('R29sf:pick', '%s %s %s %s: %d rows', ds, p, v, arm, height(r)); end
end

function m = cenMean_local(x, g)
mm = splitapply(@mean, x, g);  m = mm(g);
end
