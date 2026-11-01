function out = analyseR29_v001()
%% analyseR29_v001  R29 (spec 10.52): is it the noise LEVEL or the noise STRUCTURE whose coupling to tempo holds the pipeline-made beta_obs slope down?
% Needs the saved runs for AXbF, AXpF (R23) and AXvF, AXwF (R29, src/runR29HPC_v001.m) for Fraser and Hickman PLAC.
% R29a: slopes from analyseTempoSlopeRecovery_v002 (V1, per dataset, as analyseR23_v001); contrasts with SE = sqrt(sum of arm variances), |z| >= 2 = interval excludes zero.
% R29b: within-session corr(lfw, noise level / lsa / alpha) per arm, and the smoke-test alpha check from the saved runs.
% Writes results/tempoSim/analyseR29_v001.mat/.txt   Fraser, D.S. (2026) v001
srcDir = fileparts(mfilename('fullpath'));  ROOT = fileparts(srcDir);  addpath(srcDir);
simDir = fullfile(ROOT, 'results', 'tempoSim');
txt = fullfile(simDir, 'analyseR29_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR29_v001  %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
DS = ["Fraser" "Hickman_PLAC"];  ARMS = ["AXbF" "AXpF" "AXvF" "AXwF"];
PL = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
CASES = {"Fraser","SG-IRLS"; "Fraser","BWFD-OLS"; "Hickman_PLAC","SG-IRLS"};
rows = table();  corrRows = table();  alphaRows = table();
for ds = DS
    Fit = struct();
    for arm = ARMS
        Fit.(arm) = analyseTempoSlopeRecovery_v002(arm, Datasets=ds, Tag="PerDS_"+ds).slopes;
    end
    for c = find(string(CASES(:,1))' == ds)
        p = string(CASES{c,2});
        g = @(arm) pick_local(Fit.(arm), p, ds, arm);
        sB = g("AXbF"); sP = g("AXpF"); sV = g("AXvF"); sW = g("AXwF");
        cz = @(a, b) deal(a.b - b.b, sqrt(a.SE^2 + b.SE^2));
        [Dl, seD] = cz(sP, sB);  [L, seL] = cz(sW, sB);  [T, seT] = cz(sV, sB);
        zD = Dl/seD;  zL = L/seL;  zT = T/seT;
        l = L/Dl;  t = T/Dl;
        if abs(zD) < 2, lab = "NOT APPLICABLE (Delta interval includes zero)";
        elseif abs(zL) >= 2 && abs(zT) < 2 && l >= 0.65 && t <= 0.35, lab = "LEVEL carries it";
        elseif abs(zT) >= 2 && abs(zL) < 2 && t >= 0.65 && l <= 0.35, lab = "STRUCTURE carries it";
        elseif abs(zL) >= 2 && abs(zT) >= 2, lab = "BOTH";
        elseif abs(zL) < 2 && abs(zT) < 2, lab = "NOT RESOLVED";
        else, lab = "PARTIAL (one factor excludes zero but its share is outside the 0.65 / 0.35 rule)"; end
        nonAdd = abs(l + t - 1) > 0.35;
        rows = [rows; table(ds, p, sB.b, sP.b, sV.b, sW.b, Dl, zD, L, zL, T, zT, l, t, nonAdd, lab, ...
            'VariableNames', ["dataset","pipeline","sXb","sXp","sXv","sXw","Delta","zDelta","L","zL","T","zT","lShare","tShare","nonAdditive","label"])]; %#ok<AGROW>
    end
    % R29b: within-session correlations from the saved runs (SG-IRLS beta not needed; summaries are pipeline-independent)
    for arm = ARMS
        f = fullfile(simDir, "loopClosureResults_SYN-" + arm + "-" + ds + "_all_shaped_xu_syn001.mat");
        if ~isfile(f), error('R29:input', 'Missing: %s', f); end
        r = load(f, 'results').results;
        t = table(string({r.subjectID})', arrayfun(@(x) x.f0, r)', arrayfun(@(x) x.sigmaMM, r)', arrayfun(@(x) x.a_mm, r)', arrayfun(@(x) x.alphaIRA, r)', ...
            'VariableNames', ["session","f0","sigma","a","alpha"]);
        lf = log2(t.f0);  lsa = log2(t.sigma ./ t.a);  ok = isfinite(lf) & isfinite(lsa) & isfinite(t.alpha) & t.f0 > 0;
        t = t(ok,:);  lf = lf(ok);  lsa = lsa(ok);  gi = findgroups(t.session);  cen = @(x) x - cenMean_local(x, gi);
        corrRows = [corrRows; table(ds, arm, corr(cen(lf), cen(lsa)), corr(cen(lf), cen(t.alpha)), median(t.alpha), ...
            'VariableNames', ["dataset","arm","corrLfwLsa","corrLfwAlpha","medianAlpha"])]; %#ok<AGROW>
        S = load(fullfile(simDir, "tempoSim_" + arm + "_" + ds + ".mat"), 'dDiag');
        cr = NaN;  if isfield(S.dDiag, 'corrReal'), cr = S.dDiag.corrReal; end   % v002 corpora (AXbF, AXpF) did not record the realised level: NaN, disclosed
        alphaRows = [alphaRows; table(ds, arm, cr, 'VariableNames', ["dataset","arm","corrLevelReal"])]; %#ok<AGROW>
    end
end
fprintf('\nR29a (slope per doubling of beta_obs; Delta = Xp-Xb, L = Xw-Xb level decoupled, T = Xv-Xb structure decoupled):\n'); disp(rows);
fprintf('R29b (within-session correlations with true log2 f0; corrLevelReal = noise level actually inserted; AXbF has no dDiag value):\n'); disp(corrRows); disp(alphaRows);
hy = table();
for ds = DS
    a = @(arm) corrRows.corrLfwAlpha(corrRows.dataset == ds & corrRows.arm == arm);
    dv = a("AXvF") - a("AXbF");  dw = a("AXwF") - a("AXbF");
    if abs(dv) > 0.10 && abs(dw) <= 0.10, h = "SUPPORTED (donor structure re-couples alpha)";
    elseif abs(dw) > 0.10 && abs(dv) <= 0.10, h = "CONTRADICTED (level scaling re-couples alpha)";
    else, h = "NOT SEPARATED"; end
    m = @(arm) corrRows.medianAlpha(corrRows.dataset == ds & corrRows.arm == arm);
    colourOK = abs(m("AXvF") - m("AXpF")) <= 0.05 && abs(m("AXwF") - m("AXbF")) <= 0.05;
    hy = [hy; table(ds, dv, dw, h, colourOK, 'VariableNames', ["dataset","dAlphaCorrXv","dAlphaCorrXw","hypothesis10_51","scalingColourPreserving"])]; %#ok<AGROW>
end
disp(hy);
out = struct('R29a', rows, 'R29b', corrRows, 'level', alphaRows, 'hypothesis', hy);
save(fullfile(simDir, 'analyseR29_v001.mat'), '-struct', 'out');
fprintf('=== DONE ===\n');
end

function r = pick_local(T, p, ds, arm)
r = T(T.source == "synthetic" & T.variant == "V1" & T.pipeline == p & T.outcome == "betaObs", :);
if height(r) ~= 1, error('R29:pick', '%s %s %s: %d rows', ds, p, arm, height(r)); end
end

function m = cenMean_local(x, g)
mm = splitapply(@mean, x, g);  m = mm(g);
end
