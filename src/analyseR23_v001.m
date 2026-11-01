function out = analyseR23_v001(opts)
%ANALYSER23_V001  R23 (spec SPEC_T2 10.44 as amended by 10.45): does decoupling noise from tempo remove the pipeline-alone tempo slope?
%   out = analyseR23_v001()   needs results/tempoSim/tempoSim_<arm>_<dataset>.mat and the runner output for arms AXpF and AXbF (src/runR23HPC_v001.m).
%   Fits each new arm per dataset (analyseTempoSlopeRecovery_v002, Tag 'PerDS_<dataset>', as R21) and compares AXpF (noise from another trial of the session)
%   with the coupled baselines: PRIMARY AXbF (own residual, circularly shifted; same seam load); SECONDARY AXrF for Fraser and the saved AXr per-dataset
%   fits for PLAC and HALO (own residual, seam-free). R23e compares AXbF with the secondary baseline (seam check).
%   Writes results/tempoSim/analyseR23_v001.mat and .txt.   Fraser, D.S. (2026)  v001
arguments
    opts.Datasets (1,:) string = ["Fraser","Hickman_PLAC","Hickman_HALO"]
end
root = fileparts(fileparts(mfilename('fullpath'))); simDir = fullfile(root,'results','tempoSim');
if get(0,'Diary') == "off", diary(fullfile(simDir,'analyseR23_v001.txt')); cleanup = onCleanup(@() diary('off')); end %#ok<NASGU>
fprintf('analyseR23_v001 | %s | root %s\n', string(datetime('now')), root);
Rows = table(); Chk = table(); Fit = struct();
for ds = opts.Datasets
    Fit.(ds).AXpF = analyseTempoSlopeRecovery_v002("AXpF", Datasets=ds, Tag="PerDS_"+ds).slopes;
    Fit.(ds).AXbF = analyseTempoSlopeRecovery_v002("AXbF", Datasets=ds, Tag="PerDS_"+ds).slopes;
    if ds == "Fraser", bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrF.mat"); else, bf = fullfile(simDir,"analyseTempoSlopeRecovery_AXrPerDS_"+ds+".mat"); end
    if ~isfile(bf), error('R23:baseline','Missing secondary baseline %s',bf); end
    Fit.(ds).AXr = load(bf,'out').out.slopes;
    % ---- R23a checks from the AXpF corpus
    C = load(fullfile(simDir, "tempoSim_AXpF_"+ds+".mat"), 'dDiag','truth');
    nT = height(C.truth); nS = numel(unique(C.truth.session)); skipFrac = C.dDiag.nSingleton / (nT + C.dDiag.nSingleton);
    tol = max(0.10, 1.96/sqrt(nT - nS));   % amendment 2 (10.45): not distinguishable from zero
    Chk = [Chk; table(ds, nT, nS, C.dDiag.nSingleton, skipFrac, C.dDiag.nWrapped/nT, C.dDiag.corrOwn, C.dDiag.corrDonor, tol, skipFrac > 0.05, abs(C.dDiag.corrDonor) < tol, abs(C.dDiag.corrOwn) >= tol, ...
        'VariableNames',["dataset","nTrials","nSessions","nSingletonSkipped","skipFraction","wrappedFraction","rOwn","rDonor","tolerance","void","decouplingAccepted","couplingPresent"])]; %#ok<AGROW>
    for p = ["SG-IRLS","BWFD-OLS"]
        for oc = ["betaObs","betaGen"]
            pick = @(T) T(T.source=="synthetic" & T.variant=="V1" & T.pipeline==p & T.outcome==oc, :);
            d = pick(Fit.(ds).AXpF); cP = pick(Fit.(ds).AXbF); cS = pick(Fit.(ds).AXr); bx = pick(Fit.(ds).AXbF);
            if any([height(d) height(cP) height(cS)] ~= 1), error('R23:pick','%s %s %s: row counts %d %d %d',ds,p,oc,height(d),height(cP),height(cS)); end
            Rows = [Rows; cmp_local(ds,p,oc,"AXpF vs AXbF (primary)",cP,d)]; %#ok<AGROW>
            Rows = [Rows; cmp_local(ds,p,oc,"AXpF vs AXr (secondary)",cS,d)]; %#ok<AGROW>
            Rows = [Rows; cmp_local(ds,p,oc,"AXbF vs AXr (R23e seam check)",cS,bx)]; %#ok<AGROW>
        end
    end
end
Rows.label = strings(height(Rows),1);
for i = 1:height(Rows)
    r = Rows(i,:); k = Chk(Chk.dataset==r.dataset,:);
    if startsWith(r.comparison,"AXbF vs AXr"), Rows.label(i) = seamLabel_local(r); continue; end
    if k.void, Rows.label(i) = "VOID (more than 5 per cent of trials without a donor)"; continue; end
    if ~k.decouplingAccepted, Rows.label(i) = "NOT READ (decoupling not accepted)"; continue; end
    if r.bC <= 0, Rows.label(i) = "NOT APPLICABLE (baseline slope not positive)"; continue; end
    aSup = r.bD <= 0.5*r.bC && r.z <= -2;  bSup = r.bD >= 0.8*r.bC && abs(r.z) < 2;
    if aSup, lab = "(a) noise-tempo coupling produces the slope"; elseif bSup, lab = "(b) the slope survives decoupling"; else, lab = "PARTIAL or UNRESOLVED"; end
    if ~k.couplingPresent, lab = lab + " [little noise-tempo coupling to remove]"; end
    Rows.label(i) = lab;
end
fprintf('\nR23a checks (10.44, amended 10.45):\n'); disp(Chk)
fprintf('R23b, beta_obs lfw slope, decoupled (D) vs coupled (C); Delta = D - C, z = Delta/SE_Delta:\n'); disp(Rows(Rows.outcome=="betaObs",:))
fprintf('R23c, null beta_gen* slope (inversion offset), read only where the primary coupled |offset| >= 0.01:\n');
big = Rows.outcome=="betaGen" & contains(Rows.comparison,"primary") & abs(Rows.bC) >= 0.01; disp(Rows(big,:))
fprintf('\nREADINGS (rules fixed in spec 10.44 and 10.45 before the run):\n');
for ds = opts.Datasets
    for p = ["SG-IRLS","BWFD-OLS"]
        if p == "BWFD-OLS" && ds ~= "Fraser", continue; end
        a = Rows(Rows.dataset==ds & Rows.pipeline==p & Rows.outcome=="betaObs" & contains(Rows.comparison,"primary"),:);
        s = Rows(Rows.dataset==ds & Rows.pipeline==p & Rows.outcome=="betaObs" & contains(Rows.comparison,"secondary"),:);
        q = Rows(Rows.dataset==ds & Rows.pipeline==p & Rows.outcome=="betaObs" & contains(Rows.comparison,"seam"),:);
        cat1 = extractBefore(a.label+" ", " "); cat2 = extractBefore(s.label+" ", " ");
        fin = a.label; if cat1 ~= cat2, fin = "SEAM-SENSITIVE (primary and secondary baselines disagree): UNRESOLVED"; end
        fprintf('  %-13s %-8s AXbF %+.4f -> AXpF %+.4f (ratio %.2f, z %+.2f) | vs AXr %+.4f (z %+.2f) | seam AXbF-AXr z %+.2f (%s)\n    -> %s\n', ...
            ds, p, a.bC, a.bD, a.ratio, a.z, s.bC, s.z, q.z, q.label, fin);
    end
end
out = struct('Rows',Rows,'Checks',Chk,'runDate',datetime('now'),'matlab',version);
save(fullfile(simDir,'analyseR23_v001.mat'),'out');
end

function T = cmp_local(ds, p, oc, name, c, d)
dl = d.b - c.b; se = sqrt(c.SE^2 + d.SE^2);
T = table(ds, p, oc, name, c.n, c.b, c.SE, d.n, d.b, d.SE, dl, se, dl/se, d.b/c.b, ...
    'VariableNames',["dataset","pipeline","outcome","comparison","nC","bC","seC","nD","bD","seD","delta","seDelta","z","ratio"]);
end

function s = seamLabel_local(r)
if abs(r.z) < 2, s = "no detectable seam effect (|z| < 2)"; else, s = "SEAM EFFECT detectable (|z| >= 2)"; end
end
