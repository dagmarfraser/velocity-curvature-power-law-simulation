function out = simulateTempoCorpus_v002(arm, opts)
%SIMULATETEMPOCORPUS_V002  = v001 plus the Xp noise suffix (spec SPEC_T2 10.44, R23): the trial's noise is ANOTHER trial's own residual from the same
%   session (random derangement, cropped or wrapped to length), i.e. noise independent of tempo within a session. Every other arm is built exactly as in v001.
%SIMULATETEMPOCORPUS_V001  Synthetic recordings with KNOWN beta_gen for T2/T2b.
%   out = simulateTempoCorpus_v001("A")  -- one arm, all six in-domain datasets.
%
%   Design: docs/SPEC_T2_TempoSlopeRecovery_v001.md. Each synthetic trial copies an
%   EMPIRICAL trial's coordinates (f0, a, b, theta, residual, alphaMaj/Min) from the
%   v015 shaped_xu corpus, replaces beta_gen with a KNOWN truth, and builds the record
%   exactly as the runner's forward map does (template at true f0 + shaped_xu noise
%   drawn from the empirical residual). The runner then analyses it blind, unchanged.
%
%   ARMS  A     null: dataset median + session intercept, no tempo term
%         B1    A + 0.0265 * within-session log2 f0 deviation
%         B2    A + 0.084  * same
%         E     A, but TRUE f0 constant within session (estimator coupling only)
%         Sm2 Sm1 Sp1 Sp2   A shifted by -2,-1,+1,+2 between-session SDs (plateau sweep)
%         SD01 SD05         A with session SD fixed at 0.01 / 0.05
%   (Gate-off arm G is NOT built: the runner returns NaN beta_gen* when gated out.)
%
%   Known approximations (disclosed, not hidden):
%    - Synthetic record length is the empirical TEMPLATE-SUBTRACTED residual length M
%      (shorter than the raw record by the window trim), so the runner trims again.
%    - Truth f0 is the empirical ESTIMATED f0 (treated as true).
%    - Trial-level truth scatter is opts.SdTrial (default 0); compare recovered
%      within-session SD with the empirical value in out.diag before trusting an arm.
%
%   Fails loud on missing v015 file, unmatched trial IDs, or a noise length mismatch.
arguments
    arm (1,1) string
    opts.Datasets (1,:) string = ["Fraser","Cook CTRL","Cook ASD","Hickman HALO","Hickman PLAC","Zarandi"]
    opts.FraserSessions (1,1) double = 10     % whole sessions kept for Fraser
    opts.SdTrial (1,1) double = 0
    opts.Seed (1,1) double = 20261014
    opts.MaxTrials (1,1) double = Inf         % smoke tests only: whole sessions up to this many trials
    opts.OutDir (1,1) string = ""
end
srcDir = fileparts(mfilename('fullpath'));
addpath(srcDir); addpath(genpath(fullfile(srcDir,'functions'))); addpath(genpath(fullfile(srcDir,'req')));
if opts.OutDir == "", opts.OutDir = fullfile(srcDir,'..','results','tempoSim'); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

% trailing F on the arm name = ALL Fraser sessions (spec 10.24); the full name still feeds the seed and the trial IDs
armBase = arm;  if endsWith(arm,"F") && strlength(arm) > 1, armBase = extractBefore(arm, strlength(arm)); opts.FraserSessions = Inf; end
P = armParams_local(armBase);
out = struct('arm',arm,'files',strings(0,1),'diag',table());
D = table();
for ds = opts.Datasets
    rng(opts.Seed + sum(char(ds)) + sum(char(arm)), 'twister');
    [importFn, sigToMM, FS] = datasetCfg_local(ds);
    tag = strrep(ds,' ','_');
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v015.mat', tag));
    if ~isfile(f), error('simulateTempoCorpus:noV015','Missing %s',f); end
    R = load(f,'results').results;
    rid = string({R.trialID})'; sid = string({R.subjectID})';
    bgm = arrayfun(@(r) r.betaGenStarMed, R)';
    f0e = arrayfun(@(r) r.f0, R)';
    keep = true(size(rid));
    if ds == "Fraser"
        us = unique(sid); us = us(randperm(numel(us)));
        keep = ismember(sid, us(1:min(opts.FraserSessions,numel(us))));
    end
    if isfinite(opts.MaxTrials)
        uk = unique(sid(keep)); nk = arrayfun(@(s) sum(sid==s & keep), uk);
        [~,o] = sort(nk); cs = cumsum(nk(o)); pick = uk(o(cs <= max(opts.MaxTrials, nk(o(1)))));
        keep = keep & ismember(sid,pick);
    end
    % empirical variance components (all trials, not the Fraser subsample)
    [sdB, medB, sdW] = varComp_local(bgm, sid);
    if ~isnan(P.sdSess), sdB = P.sdSess; end
    trAll = importFn(); tid = string({trAll.trialID})';
    idx = find(keep);
    % Trials whose empirical IRASA alpha failed have no shaped_xu noise model (the runner itself
    % returned NaN for them), so they cannot be simulated. Skipped LOUDLY and counted.
    aOK = arrayfun(@(r) isfinite(r.alphaMaj) && isfinite(r.alphaMin), R(idx));
    nSkip = sum(~aOK);
    if nSkip > 0
        warning('simulateTempoCorpus:noAlpha','%s: %d of %d trials skipped (empirical alphaMaj/Min not finite).',ds,nSkip,numel(idx));
    end
    idx = idx(aOK);
    [tf, loc] = ismember(rid(idx), tid);
    if ~all(tf), error('simulateTempoCorpus:unmatched','%s: %d empirical trials not matched to import', ds, sum(~tf)); end
    % --- Xp (spec 10.44): noise = ANOTHER trial's own residual from the SAME session, so noise is independent of tempo within a session ---
    RES = {}; donorOf = []; dStream = []; origPos = (1:numel(idx))'; idx0 = idx; nSingle = 0;
    if strcmp(P.noise,'realperm')
        RES = cell(numel(idx),1);
        for k = 1:numel(idx)
            r = R(idx(k)); tr = trAll(loc(k));
            [xr, yr] = templateSubtract_local(double(tr.x), double(tr.y), FS, r.f0, 4, 4);
            RES{k} = [xr*cos(r.theta)+yr*sin(r.theta), -xr*sin(r.theta)+yr*cos(r.theta)]; %#ok<AGROW>
        end
        dStream = RandStream('twister','Seed',mod(opts.Seed + sum(char(ds)) + sum(char(arm)) + 7919, 2^32-1));   % own stream: the global rng sequence is untouched
        [~,~,g0] = unique(sid(idx));
        donorOf = (1:numel(idx))'; isSingle = false(numel(idx),1);
        for s0 = 1:max(g0)
            m0 = find(g0==s0);
            if numel(m0) < 2, isSingle(m0) = true; continue; end
            pm = m0(randperm(dStream,numel(m0)));        % random order; each trial takes the residual of its predecessor: a derangement
            donorOf(pm) = circshift(pm,1);
        end
        nSingle = sum(isSingle);
        if nSingle > 0, warning('simulateTempoCorpus:singleton','%s: %d of %d trials skipped (their session has one trial, so no donor).',ds,nSingle,numel(idx)); end
        origPos = find(~isSingle);  idx = idx0(origPos);  loc = loc(origPos);
    end
    wrapped = false(numel(idx),1);
    % session structure for the kept trials
    sK = sid(idx); [us2,~,g] = unique(sK);
    u = sdB*randn(numel(us2),1);
    lf = log2(f0e(idx)); lfSessMean = accumarray(g,lf,[],@mean);
    if P.constTempo, f0True = 2.^lfSessMean(g); else, f0True = f0e(idx); end
    lfTrue = log2(f0True); lfSm = accumarray(g,lfTrue,[],@mean);
    lfw = lfTrue - lfSm(g);
    bTrue = medB + P.shift*sdB + u(g) + P.slope*lfw + opts.SdTrial*randn(numel(idx),1);
    bTrue = min(max(bTrue,0.02),0.73);
    trials = struct('trialID',{},'subjectID',{},'x',{},'y',{},'fs',{});
    for k = 1:numel(idx)
        r = R(idx(k)); tr = trAll(loc(k));
        [xr, yr] = templateSubtract_local(double(tr.x), double(tr.y), FS, r.f0, 4, 4);
        M = numel(xr);
        th = r.theta; rm = xr*cos(th)+yr*sin(th); rn = -xr*sin(th)+yr*cos(th);
        if strcmp(P.noise,'realperm'), [nMaj, nMin, wrapped(k)] = donorNoise_local(RES{donorOf(origPos(k))}, M, dStream);   % Xp: donor residual, same session
        elseif strcmp(P.noise,'real'), nMaj = rm(:); nMin = rn(:);   % Xr: the trial's own residual, unshifted
        else, [nMaj, nMin] = generateLoopClosureNoise_v003(P.noise, rm(:), rn(:), FS, r.alphaMaj, r.alphaMin); end
        if numel(nMaj) ~= M || any(~isfinite(nMaj)) || any(~isfinite(nMin))
            error('simulateTempoCorpus:noise','%s trial %s: noise length %d vs M %d or non-finite',ds,r.trialID,numel(nMaj),M);
        end
        aN = r.a_mm/sigToMM; bN = r.b_mm/sigToMM;
        [xT, yT] = template_local(aN,bN,th,bTrue(k),f0True(k),M,FS);
        trials(k).trialID   = char("SYN-"+arm+"-"+r.trialID);
        trials(k).subjectID = char(r.subjectID);
        trials(k).x = xT + nMaj*cos(th) - nMin*sin(th);
        trials(k).y = yT + nMaj*sin(th) + nMin*cos(th);
        trials(k).fs = FS;
    end
    dDiag = struct('nSingleton',nSingle,'nWrapped',NaN,'corrOwn',NaN,'corrDonor',NaN);
    if strcmp(P.noise,'realperm')   % coupling index (spec 10.44 R23a ii): within-session correlation of noise level with true tempo deviation
        lvl = @(Z) 0.5*log(var(Z(:,1)) + var(Z(:,2)));
        Lown = cellfun(lvl, RES(origPos)); Ldon = cellfun(lvl, RES(donorOf(origPos)));
        dDiag.nWrapped = sum(wrapped);
        dDiag.corrOwn = withinCorr_local(Lown, lfw, g); dDiag.corrDonor = withinCorr_local(Ldon, lfw, g);
        fprintf('[%s] %-13s Xp: %d wrapped, %d singleton-skipped | within-session corr(noise level, tempo): own residual %+.3f, donor %+.3f\n', arm, ds, dDiag.nWrapped, nSingle, dDiag.corrOwn, dDiag.corrDonor);
    end
    truth = table(string({trials.trialID})', sK, f0True, f0e(idx), bTrue, lfw, ...
        'VariableNames',["trialID","session","f0True","f0Emp","betaTrue","lfwTrue"]);
    if strcmp(P.noise,'realperm'), truth.donorTrialID = rid(idx0(donorOf(origPos))); truth.wrapped = wrapped; end
    % Noise covariates ESTIMATED FROM THE SYNTHETIC RECORD by the same characterisation routine as the
    % empirical bioResults (template residual IRASA alpha, template sigma), NOT copied from the source trial.
    bioResults = characteriseBiologicalNoise_v004(trials, Method="irasa", Verbose=false);
    if ~all(isfinite(bioResults.ira_alphaMean)), warning('simulateTempoCorpus:alphaNaN','%s: %d synthetic trials with non-finite estimated alpha (runner will gate them out).',ds,sum(~isfinite(bioResults.ira_alphaMean))); end
    outFile = fullfile(opts.OutDir, sprintf('tempoSim_%s_%s.mat',arm,tag));
    maxTrials = opts.MaxTrials; %#ok<NASGU>
    covariatesEstimated = true; %#ok<NASGU>
    save(outFile,'trials','bioResults','truth','sigToMM','FS','P','sdB','medB','sdW','nSkip','maxTrials','covariatesEstimated','dDiag','-v7.3');
    out.files(end+1,1) = outFile;
    D = [D; table(ds, numel(idx), numel(us2), medB, sdB, sdW, nSkip, dDiag.nSingleton, dDiag.nWrapped, dDiag.corrOwn, dDiag.corrDonor, 'VariableNames', ...
        ["dataset","nTrials","nSessions","medianTruth","sdSessionUsed","sdWithinEmp","nSkippedNoAlpha","nSkippedSingleton","nWrapped","corrNoiseTempoOwn","corrNoiseTempoDonor"])]; %#ok<AGROW>
    fprintf('[%s] %-13s %4d trials %3d sessions  median %.3f  sdSess %.3f\n',arm,ds,numel(idx),numel(us2),medB,sdB);
end
out.diag = D;
end

% ------------------------------------------------------------------------
function P = armParams_local(arm)
P = struct('slope',0,'shift',0,'sdSess',NaN,'constTempo',false,'noise','shaped_xu');
% reseeded replicate '<base>s<k>' (As2, B1s2): same design as <base>; the arm name feeds the seed so the draw differs
tok = regexp(char(arm),'^(.+?)s(\d+)$','tokens','once'); if ~isempty(tok), arm = string(tok{1}); end
% generator suffix on any base arm: Xb = bootstrap (circularly shifted empirical residual), Xu = unshaped Xu,
% Xr = the empirical residual itself, UNSHIFTED (real noise, no wrap seam): the closest null to the empirical situation
sfx = ["Xb","bootstrap"; "Xu","xu"; "Xr","real"; "Xp","realperm"];
for q = 1:size(sfx,1), if endsWith(arm, sfx(q,1)) && strlength(arm) > 2, P.noise = char(sfx(q,2)); arm = extractBefore(arm, strlength(arm)-1); end, end
switch arm
    case "A"
    case "B1",  P.slope = 0.0265;
    case "B2",  P.slope = 0.084;
    case "E",   P.constTempo = true;
    case "Sm2", P.shift = -2;
    case "Sm1", P.shift = -1;
    case "Sp1", P.shift = 1;
    case "Sp2", P.shift = 2;
    case "SD01", P.sdSess = 0.01;
    case "SD05", P.sdSess = 0.05;
    case "B3",   P.slope = 0.05;          % dose-response point between B1 and B2
    otherwise, error('simulateTempoCorpus:arm','Unknown or unbuilt arm "%s" (gate-off G is not built).',arm);
end
end

function [sdB, med, sdW] = varComp_local(b, sid)
% shrunken between-session SD: var(session means) - mean(within var / n); floor 0.005
ok = isfinite(b); [~,~,g] = unique(sid(ok)); bb = b(ok);
m = accumarray(g,bb,[],@mean); n = accumarray(g,1); w = bb - m(g);
sdW = std(w); sw2 = var(w)*numel(w)/max(numel(w)-numel(m),1);
sdB = sqrt(max(var(m(n>1)) - mean(sw2./n(n>1)), 0.005^2));
med = median(bb);
end

function [x, y] = template_local(a,b,theta,beta,f0,M,FS)
% identical construction to runner processOneTrial_local L882-897
phi = linspace(0,2*pi,10000)';
dsdp = sqrt((a.*sin(phi)).^2 + (b.*cos(phi)).^2);
den = (a^2.*sin(phi).^2 + b^2.*cos(phi).^2).^1.5;
kp = (a*b)./max(den,eps);
cumT = cumsum(dsdp.*kp.^beta); cumT = cumT/cumT(end);
tN = min(mod((0:M-1)'/FS*f0,1),1-1e-9);
phiA = interp1(cumT,phi,tN,'linear','extrap');
xE = a*cos(phiA); yE = b*sin(phiA);
x = xE*cos(theta)-yE*sin(theta); y = xE*sin(theta)+yE*cos(theta);
end

function [importFn, sigToMM, FS] = datasetCfg_local(ds)
switch ds   % mirrors runLoopClosureFftnoise_v013 L170-215
    case "Cook CTRL",    sigToMM = 0.248; FS = 133; importFn = @() importDB_cook_v002(Group="CTRL",Tasks=7,Verbose=false);
    case "Cook ASD",     sigToMM = 0.248; FS = 133; importFn = @() importDB_cook_v002(Group="ASD",Tasks=7,Verbose=false);
    case "Hickman PLAC", sigToMM = 0.248; FS = 133; importFn = @() importDB_hickman_v003(Study=2,Group="PLAC",Verbose=false);
    case "Hickman HALO", sigToMM = 0.248; FS = 133; importFn = @() importDB_hickman_v003(Study=2,Group="HALO",Verbose=false);
    case "Zarandi",      sigToMM = 10.0;  FS = 100; importFn = @() importDB_zarandi_v001(Verbose=false);
    case "Fraser",       sigToMM = 1/10.41793; FS = 240; importFn = @() importDB_fraser_v001(Verbose=false);
    otherwise, error('simulateTempoCorpus:dataset','Unknown dataset %s',ds);
end
end

function [resX, resY] = templateSubtract_local(x, y, fs, f0, nH, nCW)
% copy of runner templateSubtract_local (v013 L1000-1035), residual only
N = numel(x); win = round(nCW/f0*fs); hop = round(win/2);
win = min(win,N); win = max(win,round(2/f0*fs));
resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1); s = 1;
while s+win-1 <= N
    e = s+win-1; tw = (0:win-1)'/fs; han = 0.5*(1-cos(2*pi*(0:win-1)'/(win-1)));
    Dm = [ones(win,1), tw];
    for h = 1:nH, Dm(:,end+1)=cos(2*pi*h*f0*tw); Dm(:,end+1)=sin(2*pi*h*f0*tw); end %#ok<AGROW>
    Dw = Dm.*han; bX = Dw\(x(s:e).*han); bY = Dw\(y(s:e).*han);
    resX(s:e) = resX(s:e)+(x(s:e)-Dm*bX).*han; resY(s:e) = resY(s:e)+(y(s:e)-Dm*bY).*han;
    ww(s:e) = ww(s:e)+han; s = s+hop;
end
ok = ww>0; resX(ok)=resX(ok)./ww(ok); resY(ok)=resY(ok)./ww(ok);
cl = max(round(win/4),1); cs = cl; ce = min(max(N-cl,cs+100),N);
resX = resX(cs:ce); resY = resY(cs:ce);
end

function [nMaj, nMin, didWrap] = donorNoise_local(Z, M, stream)
% Z = donor residual (Md x 2, along and across the DONOR's axes); cropped from a random start if long enough, else wrapped circularly (counted)
Md = size(Z,1);
if Md >= M, off = randi(stream, Md-M+1) - 1; didWrap = false; else, off = randi(stream, Md) - 1; didWrap = true; end
ii = mod(off + (0:M-1)', Md) + 1;
nMaj = Z(ii,1); nMin = Z(ii,2);
end

function r = withinCorr_local(L, x, g)
mL = accumarray(g,L,[],@mean); mx = accumarray(g,x,[],@mean);
r = corr(L - mL(g), x - mx(g));
end
