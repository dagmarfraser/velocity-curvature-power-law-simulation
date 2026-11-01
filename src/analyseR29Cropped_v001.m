function out = analyseR29Cropped_v001()
%% analyseR29Cropped_v001  R29 robustness (1) of spec 10.56: the R29 contrasts refitted on CROPPED-ONLY trials (the trials whose Xp/Xv/Xw donor segment was NOT wrapped).
% Rule (10.56, written before this ran): same arms, contrasts, labels and thresholds as analyseR29_v001 (Delta = Xp-Xb, L = Xw-Xb, T = Xv-Xb; |z| >= 2; l,t shares 0.65/0.35), but on the
% trial subset whose Xp donor was cropped; the SAME subset (by source trial) is used for Xb, Xp, Xv, Xw. V1 and V3 reported. Also reported: trials kept, SD of within-session log2 f0 kept
% (wrapping tracks tempo, so the subset has a narrower tempo range: this is a SENSITIVITY check, the seam-free rebuild is the clean test).
% Writes results/tempoSim/analyseR29Cropped_v001.mat/.txt   Fraser, D.S. (2026) v001
srcDir = fileparts(mfilename('fullpath'));  ROOT = fileparts(srcDir);  addpath(srcDir); addpath(genpath(fullfile(srcDir,'functions'))); addpath(genpath(fullfile(srcDir,'req')));
simDir = fullfile(ROOT, 'results', 'tempoSim');
txt = fullfile(simDir, 'analyseR29Cropped_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR29Cropped_v001  %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
DS = ["Fraser" "Hickman_PLAC"];  ARMS = ["AXbF" "AXpF" "AXvF" "AXwF"];  VAR = ["V1" "V3"];
CASES = {"Fraser","SG-IRLS"; "Fraser","BWFD-OLS"; "Hickman_PLAC","SG-IRLS"};
rows = table();  keepRows = table();
for ds = DS
    tr = load(fullfile(simDir, "tempoSim_AXpF_" + ds + ".mat"), 'truth').truth;
    cropSrc = extractAfter(string(tr.trialID(~tr.wrapped)), "SYN-AXpF-");
    Tv = struct();
    for arm = ARMS
        o = analyseTempoSlopeRecovery_v002(arm, Datasets=ds, Tag="PerDS_" + ds + "_crop", KeepTable=true);
        T = o.Tv;  Tv.(arm) = T(ismember(string(T.trialSrc), cropSrc), :);
    end
    for c = find(string(CASES(:,1))' == ds)
        p = string(CASES{c,2});
        for v = VAR
            Fit = struct();
            for arm = ARMS, Fit.(arm) = fit_local(Tv.(arm), p, v, ds); end
            sB = Fit.AXbF; sP = Fit.AXpF; sV = Fit.AXvF; sW = Fit.AXwF;
            cz = @(a, b) deal(a.b - b.b, sqrt(a.SE^2 + b.SE^2));
            [Dl, seD] = cz(sP, sB);  [L, seL] = cz(sW, sB);  [Tt, seT] = cz(sV, sB);
            zD = Dl/seD;  zL = L/seL;  zT = Tt/seT;  l = L/Dl;  t = Tt/Dl;
            if abs(zD) < 2, lab = "NOT APPLICABLE (Delta interval includes zero)";
            elseif abs(zL) >= 2 && abs(zT) < 2 && l >= 0.65 && t <= 0.35, lab = "LEVEL carries it";
            elseif abs(zT) >= 2 && abs(zL) < 2 && t >= 0.65 && l <= 0.35, lab = "STRUCTURE carries it";
            elseif abs(zL) >= 2 && abs(zT) >= 2, lab = "BOTH";
            elseif abs(zL) < 2 && abs(zT) < 2, lab = "NOT RESOLVED";
            else, lab = "PARTIAL (one factor excludes zero but its share is outside the 0.65 / 0.35 rule)"; end
            rows = [rows; table(ds, p, v, sB.b, sP.b, sV.b, sW.b, Dl, zD, L, zL, Tt, zT, l, t, abs(l+t-1) > 0.35, lab, 'VariableNames', ...
                ["dataset","pipeline","variant","sXb","sXp","sXv","sXw","Delta","zDelta","L","zL","T","zT","lShare","tShare","nonAdditive","label"])]; %#ok<AGROW>
        end
    end
    Pall = Tv.AXpF(Tv.AXpF.pipeline == "SG-IRLS", :);  gm = splitapply(@mean, log2(Pall.f0), findgroups(Pall.session));  lfwKept = log2(Pall.f0) - gm(findgroups(Pall.session));
    keepRows = [keepRows; table(ds, height(Pall), numel(cropSrc), mean(~tr.wrapped), std(lfwKept), 'VariableNames', ["dataset","nKeptRows","nCropped","shareCropped","sdLfwKept"])]; %#ok<AGROW>
end
fprintf('\nR29 contrasts on CROPPED-ONLY trials (compare with analyseR29_v001 on all trials):\n'); disp(rows);
fprintf('Subset size and tempo range (SD of within-session log2 f0 among kept trials; the full-sample values are in spec 10.53):\n'); disp(keepRows);
out = struct('rows', rows, 'kept', keepRows);
save(fullfile(simDir, 'analyseR29Cropped_v001.mat'), '-struct', 'out');
fprintf('=== DONE ===\n');
end

function r = fit_local(T, p, variant, ds)
% the published model on one dataset (as analyseTempoSlopeRecovery_v002 fitAll_local), on the supplied trial subset
z = @(x) (x - mean(x,"omitnan")) ./ std(x,"omitnan");
W = T(T.pipeline == p & isfinite(T.betaObs), :);
W.lf = log2(W.f0);
G = groupsummary(W, "session", "mean", "lf");  [~, loc] = ismember(W.session, G.session);
W.lfw = W.lf - G.mean_lf(loc);  W.lfb = G.mean_lf(loc) - mean(G.mean_lf);
W.lsa = z(log2(W.sigma ./ W.a));  W.al = z(W.alpha);  W.sessionC = categorical(W.session);
if variant == "V3", rhs = " ~ 1 + lfw + lfb + (1|sessionC)"; else, rhs = " ~ 1 + lfw + lfb + lsa + al + (1|sessionC)"; end
m = fitlme(W, "betaObs" + rhs);  c = m.Coefficients;  k = string(c.Name) == "lfw";
if nnz(k) ~= 1, error('R29Cropped:fit', '%s %s %s: lfw coefficient not found', ds, p, variant); end
r = struct('b', c.Estimate(k), 'SE', c.SE(k), 'n', height(W));
end
