% visualiseLoopClosurePipeline_v002.m
%
% v002 changes from v001:
%   - SG-LMLS (pp=4) designated primary pipeline throughout (was SG-IRLS pp=6).
%     Matches skeleton v009 designation (zero geometric gap; consistent across
%     noise colours; Finding #85).
%   - Noise model (shaped_xu) was already correct in v001; no change needed.
%     shaped_xu = empirical |FFT| amplitude on Xu fractional phase.
%     Consequence: synthetic traces do NOT drift from the template —
%     the real residual's spectral amplitude (including the low-frequency shelf)
%     bounds the surrogate's low-frequency energy, removing the systematic drift
%     concern that motivated the orbit-centring visualisation.  The centreByOrbit
%     local function is retained for the Figure 2b comparison panel only.
%   - PATCH (in-place, 2026-07-14): Figure 2's global speed colour scale
%     (spdGlobalMin/Max) was pooled only from the 6 grid-point templates
%     (showIdx). The beta_gen* panel plots templateBetaStar, an off-grid
%     interpolated template whose speed range was never included in that
%     pool, so its colours could be silently clipped to the colourbar's
%     extreme if its speed range fell outside the 6 grid templates' range.
%     Fixed by folding templateBetaStar's speed into spdAll before taking
%     min/max, whenever USE_EXACT_STAR is true. No other behaviour changed.
%
% 8-panel walkthrough of the loop closure pipeline for one Cook CTRL trial.
%
%  Panel 1  Raw 2D trajectory, coloured by time
%  Panel 2  x(t): raw signal, OLA template, OLA residual
%  Panel 3  Lab-frame residuals xResid, yResid overlaid
%  Panel 4  PCA ellipse with major / minor axes
%  Panel 5  PCA-rotated residuals: resid_maj and resid_min
%  Panel 6  Log-log PSD of both PCA-frame residuals + IRASA slopes
%  Panel 7  Shaped-Xu surrogate PSD vs real residual PSD
%  Panel 8  Forward map: median(betaRec) vs beta_gen, all 6 pipelines,
%           with betaObs and betaGenStar marked
%
% Requires: project src/ as working directory.
%           importDB_cook_v002, differentiateKinematicsEBR, regressDataEBR,
%           curvatureKinematicEBR, iraAlphaSigma_v001, generateLoopClosureNoise_v002.
%           For Panel 8: loopClosureResults_Cook_CTRL_per_subject_shaped_xu_v007.mat
%
% Usage (from src/):
%   visualiseLoopClosurePipeline_v001
%
% Fraser, D.S.  2026

clc;
fprintf('=== visualiseLoopClosurePipeline_v001 ===\n');

%% -----------------------------------------------------------------------
%% CONFIG
%% -----------------------------------------------------------------------
TRIAL_IDX   = 3;       % which trial to visualise (1-indexed into per_subject set)
FS          = 133;
SIG_TO_MM   = 0.248;
EDGE_CLIP   = 20;
N_BETA      = 25;
BETA_RANGE  = [0 0.75];
N_REPS      = 20;      % for Panel 8 forward-map rebuild
N_OLA_H     = 4;       % OLA harmonics
N_OLA_CW    = 4;       % OLA cycle windows
IRA_FLOW    = 1.0;
IRA_FHIGH   = min(20, FS/2 - 1);
IRA_HSET    = 1.1:0.05:1.9;
VIS_ORBIT_CENTRE = true;  % VISUALISATION ONLY — do not apply in v007 forward map
N_VIS_REPS  = 12;                     % surrogates shown in all vis figures
SUR_CLR     = [0.85 0.60 0.15 0.28];  % gold, alpha=0.28 for 12-trace density

DERIV_CFG   = struct('label',  {"BWFD","SG"}, ...
                     'filterType',   {2, 6}, ...
                     'filterParams', {[2 10 1],[4 17]});
REG_CFG     = struct('label', {"OLS","LMLS","IRLS"}, 'type', {3,4,5});
PIPE_LABELS = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
PP_PRIMARY  = 4;   % SG-LMLS — primary pipeline (v002; was 6=SG-IRLS)
N_PIPES     = 6;
PIPE_CLRS   = [0.2 0.5 0.8;   % BWFD-OLS  blue
               0.1 0.7 0.4;   % SG-OLS    green
               0.8 0.4 0.1;   % BWFD-LMLS orange
               0.6 0.2 0.7;   % SG-LMLS   purple
               0.6 0.6 0.6;   % BWFD-IRLS grey
               0.9 0.2 0.2];  % SG-IRLS   red

betaGenVec  = linspace(BETA_RANGE(1), BETA_RANGE(2), N_BETA);
nD          = 2;   % number of differentiator configs (BWFD, SG)
nR          = 3;   % number of regression configs (OLS, LMLS, IRLS)
rng(1729, 'twister');

srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
addpath(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));
addpath(genpath(fullfile(srcDir, 'req')));

%% -----------------------------------------------------------------------
%% LOAD DATA  (skipped if tr already in workspace)
%% -----------------------------------------------------------------------
if ~exist('tr', 'var')
    fprintf('Loading Cook CTRL...\n');
    trials     = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
    bioMat     = load(fullfile(srcDir, 'noiseCharacterisation_cook.mat'), 'bioResults');
    bio        = bioMat.bioResults;
    alphaVec   = bio.ira_alphaMean;
    sigmaVec   = bio.sigmaMean * SIG_TO_MM;

    % Per-subject median trial selection (same logic as v007)
    subjIDs   = unique(bio.subjectID);
    selIdx    = zeros(numel(subjIDs), 1);
    for s = 1:numel(subjIDs)
        mask  = bio.subjectID == subjIDs(s);
        aV    = alphaVec(mask); sV = sigmaVec(mask);
        valid = ~isnan(aV) & ~isnan(sV);
        muA   = median(aV(valid)); muS = median(sV(valid));
        d     = (aV - muA).^2 / max(var(aV(valid)), eps) + ...
                (sV - muS).^2 / max(var(sV(valid)), eps);
        d(~valid) = Inf;
        bioRows = find(mask);
        [~, best] = min(d);
        selIdx(s) = bioRows(best);
    end
    trialIDs   = string({trials.trialID}');
    selTIDs    = bio.trialID(selIdx);
    matchIdx   = zeros(numel(selIdx), 1);
    for s = 1:numel(selIdx)
        m = find(trialIDs == selTIDs(s), 1);
        if ~isempty(m), matchIdx(s) = m; end
    end
    matchIdx   = matchIdx(matchIdx > 0);
    tr         = trials(matchIdx(TRIAL_IDX));
    alphaIRA   = alphaVec(selIdx(TRIAL_IDX));
    x          = double(tr.x);
    y          = double(tr.y);
    N          = numel(x);
    t          = (0:N-1)' / FS;
    fprintf('Trial: %s   alpha_IRA=%.3f\n', tr.trialID, alphaIRA);
else
    fprintf('Using cached trial: %s\n', tr.trialID);
end

%% -----------------------------------------------------------------------
%% PIPELINE STEPS  (skipped if xResid already in workspace)
%% -----------------------------------------------------------------------
if ~exist('xResid', 'var')

%% f0
f0 = estimateF0_lcl(x, y, FS);
fprintf('f0 = %.3f Hz\n', f0);

%% OLA template subtraction
[xResid, yResid, xFit, yFit, csOLA] = templateSubtract_lcl(x, y, FS, f0, N_OLA_H, N_OLA_CW);
M    = numel(xResid);
tRes = t(csOLA : csOLA+M-1);   % absolute time for residuals and fit

%% PCA geometry (in mm)
xMM  = x*SIG_TO_MM; xMM = xMM - mean(xMM);
yMM  = y*SIG_TO_MM; yMM = yMM - mean(yMM);
xFMM = xFit*SIG_TO_MM;
yFMM = yFit*SIG_TO_MM;
[V, D] = eig(cov(xMM, yMM));
lams   = diag(D); [lams, ord] = sort(lams, 'descend'); V = V(:, ord);
a_mm   = sqrt(2*max(lams(1),0));
b_mm   = sqrt(2*max(lams(2),0));
theta  = atan2(V(2,1), V(1,1));
fprintf('PCA: a=%.2fmm  b=%.2fmm  theta=%.1f deg\n', a_mm, b_mm, rad2deg(theta));

%% PCA-frame noise rotation
resid_maj =  xResid*cos(theta) + yResid*sin(theta);
resid_min = -xResid*sin(theta) + yResid*cos(theta);
fftMaj    = fft(resid_maj);
fftMin    = fft(resid_min);

%% IRASA per-axis alpha
fprintf('Running IRASA...\n');
IRA_FLOW  = 1.0; IRA_FHIGH = min(20, FS/2-1); IRA_HSET = 1.1:0.05:1.9;
alphaMaj  = iraAlphaSigma_v001(resid_maj, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
alphaMin  = iraAlphaSigma_v001(resid_min, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
fprintf('IRASA: alphaMaj=%.3f  alphaMin=%.3f  alphaIRA_bio=%.3f\n', ...
    alphaMaj, alphaMin, alphaIRA);

%% Observed beta from real data
betaObs_all = NaN(1, N_PIPES);
for dIdx = 1:nD
    [dx, dy] = differentiateKinematicsEBR(x, y, ...
        DERIV_CFG(dIdx).filterType, DERIV_CFG(dIdx).filterParams, FS);
    vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2); vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);
    ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3); ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);
    sp = hypot(vx, vy);
    kp = curvatureKinematicEBR(vx, vy, ax, ay);
    lm = [1, -1/3];
    for rIdx = 1:nR
        ppIdx = (rIdx-1)*nD + dIdx;
        try
            [b, v] = regressDataEBR(sp, kp, REG_CFG(rIdx).type, lm, 0, 0);
            betaObs_all(ppIdx) = b;
            if rIdx == 1, lm = [v, b]; end
        catch; end
    end
end
%% Shaped-xu surrogate for Panel 7
fprintf('Generating shaped_xu surrogate...\n');
[nMaj_sxu, nMin_sxu] = generateLoopClosureNoise_v002('shaped_xu', ...
    resid_maj, resid_min, FS, alphaMaj, alphaMin); %#ok<NASGU>

else
    fprintf('Using cached pipeline results.\n');
end % guard 2

% Always re-derive primary beta_obs (cheap; PP_PRIMARY may differ from cached session)
betaObs_primary = betaObs_all(PP_PRIMARY);
fprintf('beta_obs SG-LMLS (primary) = %.4f\n', betaObs_primary);

%% -----------------------------------------------------------------------
%% FORWARD MAP  (skipped if betaRec already in workspace)
%% -----------------------------------------------------------------------
if ~exist('betaRec', 'var')
    fprintf('Building forward map (N_REPS=%d, N_BETA=%d)...\n', N_REPS, N_BETA);
a_nat = a_mm / SIG_TO_MM; b_nat = b_mm / SIG_TO_MM;
betaRec = NaN(N_REPS, N_PIPES, N_BETA);
nPhi  = 10000;
phi   = linspace(0, 2*pi, nPhi)';
dsdp  = sqrt((a_nat.*sin(phi)).^2 + (b_nat.*cos(phi)).^2);
den   = (a_nat^2.*sin(phi).^2 + b_nat^2.*cos(phi).^2).^1.5;
kp_phi = (a_nat*b_nat) ./ max(den, eps);
tVec   = (0:M-1)' / FS;
templates = zeros(M, 2, N_BETA);
for bi = 1:N_BETA
    wt   = dsdp .* kp_phi.^betaGenVec(bi);
    cumT = cumsum(wt); cumT = cumT / cumT(end);
    tN   = min(mod(tVec*f0, 1), 1-1e-9);
    phiA = interp1(cumT, phi, tN, 'linear', 'extrap');
    xE   = a_nat*cos(phiA); yE = b_nat*sin(phiA);
    templates(:,1,bi) = xE*cos(theta) - yE*sin(theta);
    templates(:,2,bi) = xE*sin(theta) + yE*cos(theta);
end
for bi = 1:N_BETA
    xTpl = templates(:,1,bi); yTpl = templates(:,2,bi);
    for rep = 1:N_REPS
        [nMaj, nMin] = generateLoopClosureNoise_v002('shaped_xu', ...
            resid_maj, resid_min, FS, alphaMaj, alphaMin);
        xSyn = xTpl + nMaj*cos(theta) - nMin*sin(theta);
        ySyn = yTpl + nMaj*sin(theta) + nMin*cos(theta);
        for dIdx = 1:nD
            [dx, dy] = differentiateKinematicsEBR(xSyn, ySyn, ...
                DERIV_CFG(dIdx).filterType, DERIV_CFG(dIdx).filterParams, FS);
            vx = dx(EDGE_CLIP:end-EDGE_CLIP,2); vy = dy(EDGE_CLIP:end-EDGE_CLIP,2);
            ax = dx(EDGE_CLIP:end-EDGE_CLIP,3); ay = dy(EDGE_CLIP:end-EDGE_CLIP,3);
            sp = hypot(vx,vy); kp = curvatureKinematicEBR(vx,vy,ax,ay);
            lm = [1, -1/3];
            for rIdx = 1:nR
                ppIdx = (rIdx-1)*nD + dIdx;
                try
                    [b, v] = regressDataEBR(sp, kp, REG_CFG(rIdx).type, lm, 0, 0);
                    betaRec(rep, ppIdx, bi) = b;
                    if rIdx == 1, lm = [v, b]; end
                catch; end
            end
        end
    end
end

% Inversion: SG-LMLS betaGenStar (primary pipeline, v002)
medC      = squeeze(median(betaRec(:, PP_PRIMARY, :), 1, 'omitnan'));
[mS, ord] = sort(medC, 'ascend');
bgS       = betaGenVec(ord)';
fin       = isfinite(mS) & isfinite(bgS);
else
    fprintf('Using cached forward map.\n');
end % guard 3

% Always re-derive primary beta_gen* from betaRec (cheap; handles cache hits)
medC_p = squeeze(median(betaRec(:, PP_PRIMARY, :), 1, 'omitnan'));
[mS_p, ord_p] = sort(medC_p, 'ascend');
bgS_p  = betaGenVec(ord_p)';
fin_p  = isfinite(mS_p) & isfinite(bgS_p);
betaGenStar_primary = NaN;
if sum(fin_p) >= 2 && betaObs_primary >= min(mS_p(fin_p)) && ...
        betaObs_primary <= max(mS_p(fin_p))
    betaGenStar_primary = interp1(mS_p(fin_p), bgS_p(fin_p), betaObs_primary, 'linear');
end
if isfinite(betaGenStar_primary)
    fprintf('beta_gen* SG-LMLS (primary, re-derived) = %.4f\n', betaGenStar_primary);
end

%% Load stored betaGenStar for this trial (exact interpolated value from full N_REPS=50 run)
% Used for the Figure 2 "best match" panel — avoids grid-snapping artefact.
% Falls back to the locally re-derived value if the mat is unavailable.
betaGenStarStored = NaN;
lcMat = fullfile(srcDir, 'loopClosureResults_Cook_CTRL_all_shaped_xu_v007.mat');
if exist(lcMat, 'file')
    lcD = load(lcMat, 'results');
    tidMatch = find(string({lcD.results.trialID}) == string(tr.trialID), 1);
    if ~isempty(tidMatch)
        betaGenStarStored = lcD.results(tidMatch).betaGenStarMed;
        fprintf('Stored beta_gen* for %s: %.4f\n', tr.trialID, betaGenStarStored);
    else
        fprintf('Trial %s not found in stored mat — using re-derived value.\n', tr.trialID);
    end
else
    fprintf('loopClosureResults mat not found — using re-derived beta_gen*.\n');
end
% Use stored value if available; fall back to re-derived
betaGenStarDisplay = betaGenStarStored;
if ~isfinite(betaGenStarDisplay)
    betaGenStarDisplay = betaGenStar_primary;
    fprintf('Fallback: using re-derived beta_gen* = %.4f\n', betaGenStarDisplay);
end

%% Generate exact "best match" template at betaGenStarDisplay (no grid snapping)
if isfinite(betaGenStarDisplay)
    nPhi_  = 10000;
    phi_   = linspace(0, 2*pi, nPhi_)';
    dsdp_  = sqrt((a_nat*sin(phi_)).^2 + (b_nat*cos(phi_)).^2);
    den_   = (a_nat^2*sin(phi_).^2 + b_nat^2*cos(phi_).^2).^1.5;
    kp_    = (a_nat*b_nat) ./ max(den_, eps);
    wt_    = dsdp_ .* kp_.^betaGenStarDisplay;
    cumT_  = cumsum(wt_); cumT_ = cumT_ / cumT_(end);
    tN_    = min(mod((0:M-1)'/FS * f0, 1), 1-1e-9);
    phiA_  = interp1(cumT_, phi_, tN_, 'linear', 'extrap');
    xE_    = a_nat*cos(phiA_); yE_ = b_nat*sin(phiA_);
    templateBetaStar = zeros(M, 2);
    templateBetaStar(:,1) = xE_*cos(theta) - yE_*sin(theta);
    templateBetaStar(:,2) = xE_*sin(theta) + yE_*cos(theta);
    fprintf('Exact template generated at beta_gen*=%.4f\n', betaGenStarDisplay);
else
    templateBetaStar = [];
    fprintf('WARNING: betaGenStarDisplay is NaN — best-match panel will use nearest grid point.\n');
end

%% -----------------------------------------------------------------------
%% FIGURE
%% -----------------------------------------------------------------------
fprintf('\nBuilding figure...\n');
fig = findobj('Type','figure','Tag','vlc_fig1');
if isempty(fig), fig = figure('Tag','vlc_fig1'); else clf(fig); end
set(fig, 'Name','Loop Closure Pipeline', ...
    'Position', [80 40 1400 920], 'Color', 'w');
tl  = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Loop Closure Pipeline — %s   (\\alpha_{IRA}=%.2f, f_0=%.2fHz)', ...
    strrep(tr.trialID, '_', '\_'), alphaIRA, f0), ...
    'FontSize', 12, 'FontWeight', 'bold');

CLR_RAW  = [0.15 0.45 0.75];
CLR_TPL  = [0.90 0.40 0.10];
CLR_RES  = [0.20 0.65 0.30];
CLR_MAJ  = [0.10 0.40 0.80];
CLR_MIN  = [0.75 0.20 0.60];
CLR_SXU  = [0.95 0.55 0.10];
ALPHA_LN = 0.6;

% ---- Panel 1: Raw 2D trajectory coloured by speed ----------------------
ax1 = nexttile;
% Speed from central finite differences (mm/s)
dxdt = gradient(xMM) * FS;
dydt = gradient(yMM) * FS;
spd  = hypot(dxdt, dydt);
spdN = (spd - min(spd)) / max(max(spd) - min(spd), eps);  % 0-1 normalised
cmap1 = turbo(256);
hold(ax1, 'on');
for k = 1:N-1
    cidx = max(1, round(spdN(k) * 255) + 1);
    plot(ax1, xMM(k:k+1), yMM(k:k+1), '-', ...
        'Color', cmap1(cidx,:), 'LineWidth', 1.2);
end
axis(ax1, 'equal'); grid(ax1, 'on');
xlabel(ax1, 'x (mm)'); ylabel(ax1, 'y (mm)');
title(ax1, '1 · Raw trajectory (colour = speed)');
colormap(ax1, turbo); clim(ax1, [min(spd) max(spd)]);
cb = colorbar(ax1); cb.Label.String = 'speed (mm/s)';

% ---- Panel 2: x(t) raw / template / residual ---------------------------
ax2 = nexttile; hold(ax2, 'on');
% Map OLA-trimmed indices back to raw time axis for overlay
nF   = numel(xFit);  % == M (same window)
tF   = tRes;         % OLA fit and residual share the same trimmed window

plot(ax2, t,   x*SIG_TO_MM,          '-',  'Color', [CLR_RAW, ALPHA_LN],  'LineWidth', 1.2);
plot(ax2, tF,  xFMM,                  '--', 'Color', CLR_TPL,              'LineWidth', 1.5);
plot(ax2, tRes, xResid*SIG_TO_MM,    '-',  'Color', CLR_RES,              'LineWidth', 1.2);
legend(ax2, {'raw x','OLA template','residual'}, 'Location','northwest', 'FontSize', 8);
xlabel(ax2, 'time (s)'); ylabel(ax2, 'x (mm)');
title(ax2, '2 · OLA template subtraction');
grid(ax2, 'on');

% ---- Panel 3: Real trajectory vs PCA template vs one synthetic trial ---
ax3 = nexttile; hold(ax3, 'on');
% Real trajectory (trimmed to OLA window, in mm)
xRealMM = x(csOLA:csOLA+M-1) * SIG_TO_MM;
yRealMM = y(csOLA:csOLA+M-1) * SIG_TO_MM;
xRealMM = xRealMM - mean(xRealMM);
yRealMM = yRealMM - mean(yRealMM);
% PCA-ellipse template at beta_gen = 1/3
bi13 = round(interp1(betaGenVec, 1:N_BETA, 1/3, 'linear', 'extrap'));
bi13 = max(1, min(N_BETA, round(bi13)));
xTpl13 = templates(:,1,bi13) * SIG_TO_MM;
yTpl13 = templates(:,2,bi13) * SIG_TO_MM;
xTpl13 = xTpl13 - mean(xTpl13);
yTpl13 = yTpl13 - mean(yTpl13);
% One shaped_xu synthetic trial at beta_gen = 1/3
[nMaj_ex, nMin_ex] = generateLoopClosureNoise_v002('shaped_xu', ...
    resid_maj, resid_min, FS, alphaMaj, alphaMin);
xSyn13 = templates(:,1,bi13) + nMaj_ex*cos(theta) - nMin_ex*sin(theta);
ySyn13 = templates(:,2,bi13) + nMaj_ex*sin(theta) + nMin_ex*cos(theta);
xSyn13 = xSyn13*SIG_TO_MM - mean(xSyn13*SIG_TO_MM);
ySyn13 = ySyn13*SIG_TO_MM - mean(ySyn13*SIG_TO_MM);

plot(ax3, xRealMM,  yRealMM,  '-',  'Color', [CLR_RAW, 0.55], 'LineWidth', 1.0);
plot(ax3, xTpl13,   yTpl13,   '--', 'Color', CLR_TPL,          'LineWidth', 1.8);
plot(ax3, xSyn13,   ySyn13,   '-',  'Color', [CLR_SXU, 0.75],  'LineWidth', 0.9);
legend(ax3, {'real (trimmed)','PCA template (\beta=1/3)','synthetic (shaped-Xu)'}, ...
    'Location','best', 'FontSize', 8);
axis(ax3, 'equal'); grid(ax3, 'on');
xlabel(ax3, 'x (mm)'); ylabel(ax3, 'y (mm)');
title(ax3, '3 · Real vs template vs synthetic');

% ---- Panel 4: PCA ellipse with axes ------------------------------------
ax4 = nexttile; hold(ax4, 'on');
% Background trajectory coloured by speed (same spd vector as P1)
cmap4 = parula(256);
xC = xMM - mean(xMM); yC = yMM - mean(yMM);
for k = 1:N-1
    cidx = max(1, round(spdN(k)*255)+1);
    plot(ax4, xC(k:k+1), yC(k:k+1), '-', 'Color', cmap4(cidx,:), 'LineWidth', 0.8);
end
% PCA ellipse
th_ = linspace(0,2*pi,300)';
xEl = a_mm*cos(th_);
yEl = b_mm*sin(th_);
xElR = xEl*cos(theta) - yEl*sin(theta);
yElR = xEl*sin(theta) + yEl*cos(theta);
plot(ax4, xElR, yElR, '-', 'Color', CLR_RAW, 'LineWidth', 1.8);
% Major axis arrow
quiver(ax4, 0, 0, a_mm*cos(theta), a_mm*sin(theta), 0, ...
    'Color', CLR_MAJ, 'LineWidth', 2.0, 'MaxHeadSize', 0.5);
quiver(ax4, 0, 0, b_mm*cos(theta+pi/2), b_mm*sin(theta+pi/2), 0, ...
    'Color', CLR_MIN, 'LineWidth', 2.0, 'MaxHeadSize', 0.5);
text(ax4, a_mm*cos(theta)*1.12, a_mm*sin(theta)*1.12, ...
    sprintf('major\n%.1fmm', a_mm), 'Color', CLR_MAJ, 'FontSize', 8, 'HorizontalAlignment','center');
text(ax4, b_mm*cos(theta+pi/2)*1.2, b_mm*sin(theta+pi/2)*1.2, ...
    sprintf('minor\n%.1fmm', b_mm), 'Color', CLR_MIN, 'FontSize', 8, 'HorizontalAlignment','center');
axis(ax4, 'equal'); grid(ax4, 'on');
xlabel(ax4, 'x (mm)'); ylabel(ax4, 'y (mm)');
title(ax4, '4 · PCA axes (theta = {%.1f}°)', rad2deg(theta));
title(ax4, sprintf('4 · PCA axes (\\theta=%.1f°)', rad2deg(theta)));

% ---- Panel 5: PCA-rotated residuals over time --------------------------
ax5 = nexttile; hold(ax5, 'on');
plot(ax5, tRes, resid_maj*SIG_TO_MM, '-', 'Color', CLR_MAJ, 'LineWidth', 1.0);
plot(ax5, tRes, resid_min*SIG_TO_MM, '-', 'Color', CLR_MIN, 'LineWidth', 1.0);
yline(ax5,  std(resid_maj)*SIG_TO_MM, '--', 'Color', CLR_MAJ, 'Alpha', 0.5);
yline(ax5, -std(resid_maj)*SIG_TO_MM, '--', 'Color', CLR_MAJ, 'Alpha', 0.5);
yline(ax5,  std(resid_min)*SIG_TO_MM, ':',  'Color', CLR_MIN, 'Alpha', 0.5);
yline(ax5, -std(resid_min)*SIG_TO_MM, ':',  'Color', CLR_MIN, 'Alpha', 0.5);
legend(ax5, {sprintf('resid_{maj} (\\alpha=%.2f)', alphaMaj), ...
             sprintf('resid_{min} (\\alpha=%.2f)', alphaMin)}, ...
    'Location', 'northeast', 'FontSize', 8);
xlabel(ax5, 'time (s)'); ylabel(ax5, 'amplitude (mm)');
title(ax5, '5 · PCA-frame residuals');
grid(ax5, 'on');

% ---- Panel 6: Log-log PSD + IRASA slopes -------------------------------
ax6 = nexttile; hold(ax6, 'on');
nfft  = 2^nextpow2(M);
fAx   = (0:nfft/2)' * FS / nfft;
keep  = fAx > IRA_FLOW & fAx < IRA_FHIGH;
powMaj = (abs(fft(resid_maj, nfft)).^2 / (FS*M));
powMin = (abs(fft(resid_min, nfft)).^2 / (FS*M));
powMaj = powMaj(1:nfft/2+1) * 2; powMaj(1) = powMaj(1)/2;
powMin = powMin(1:nfft/2+1) * 2; powMin(1) = powMin(1)/2;

plot(ax6, fAx(keep), powMaj(keep), '-', 'Color', [CLR_MAJ, 0.55], 'LineWidth', 1.0);
plot(ax6, fAx(keep), powMin(keep), '-', 'Color', [CLR_MIN, 0.55], 'LineWidth', 1.0);

% IRASA slope lines — anchored at fFit(1) (low end of fit band)
fFit  = fAx(keep);
anchorF  = fFit(1);
anchorPM = powMaj(find(keep, 1));
anchorPm = powMin(find(keep, 1));
slopeLineMaj = anchorPM * (fFit / anchorF).^(-alphaMaj);
slopeLineMin = anchorPm * (fFit / anchorF).^(-alphaMin);
plot(ax6, fFit, slopeLineMaj, '--', 'Color', CLR_MAJ, 'LineWidth', 2.0);
plot(ax6, fFit, slopeLineMin, '--', 'Color', CLR_MIN, 'LineWidth', 2.0);

set(ax6, 'XScale','log', 'YScale','log');
legend(ax6, {'PSD maj','PSD min', ...
    sprintf('slope \\alpha=%.2f', alphaMaj), ...
    sprintf('slope \\alpha=%.2f', alphaMin)}, ...
    'Location','southwest', 'FontSize', 8);
xlabel(ax6, 'frequency (Hz)'); ylabel(ax6, 'PSD');
title(ax6, '6 · Log-log PSD + IRASA slope');
grid(ax6, 'on'); xlim(ax6, [IRA_FLOW IRA_FHIGH]);

% ---- Panel 7: Shaped-xu surrogate PSD vs real residual PSD -------------
ax7 = nexttile; hold(ax7, 'on');
% Surrogates first (behind)
surAlphas = NaN(5,1);
for k2 = 1:5
    [nM_k, ~] = generateLoopClosureNoise_v002('shaped_xu', ...
        resid_maj, resid_min, FS, alphaMaj, alphaMin);
    p2 = (abs(fft(nM_k, nfft)).^2 / (FS*M));
    p2 = p2(1:nfft/2+1)*2;
    plot(ax7, fAx(keep), p2(keep), '-', ...
        'Color', [CLR_SXU, 0.35], 'LineWidth', 0.9, 'HandleVisibility','off');
    % Quick alpha check on this surrogate
    try
        surAlphas(k2) = iraAlphaSigma_v001(nM_k, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
    catch; end
end
% Dummy handle for surrogates legend entry
plot(ax7, NaN, NaN, '-', 'Color', CLR_SXU, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('surrogates \\alpha_{mean}=%.2f', mean(surAlphas,'omitnan')));
% Real residual on top
plot(ax7, fAx(keep), powMaj(keep), '-',  'Color', CLR_MAJ, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('real resid_{maj} \\alpha=%.2f', alphaMaj));
% IRASA slope for reference
plot(ax7, fFit, slopeLineMaj, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'DisplayName', 'IRASA slope');

set(ax7, 'XScale','log', 'YScale','log');
legend(ax7, 'Location','southwest', 'FontSize', 8);
xlabel(ax7, 'frequency (Hz)'); ylabel(ax7, 'PSD');
gapStr = '';
meanSurAlpha = mean(surAlphas, 'omitnan');
if isfinite(meanSurAlpha)
    gapStr = sprintf('  (R7 gap=%.2f)', abs(meanSurAlpha - alphaMaj));
end
title(ax7, ['7 · Shaped-Xu: surrogate vs real' gapStr]);
grid(ax7, 'on'); xlim(ax7, [IRA_FLOW IRA_FHIGH]);

% ---- Panel 8: Forward map ----------------------------------------------
ax8 = nexttile; hold(ax8, 'on');
% All 6 pipelines (translucent fills + median line)
for pp = 1:N_PIPES
    medFwd  = squeeze(median(betaRec(:, pp, :), 1, 'omitnan'));
    p5Fwd   = squeeze(prctile(betaRec(:, pp, :),  5, 1));
    p95Fwd  = squeeze(prctile(betaRec(:, pp, :), 95, 1));
    fill(ax8, [betaGenVec, fliplr(betaGenVec)], [p5Fwd', fliplr(p95Fwd')], ...
        PIPE_CLRS(pp,:), 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax8, betaGenVec, medFwd, '-', ...
        'Color', PIPE_CLRS(pp,:), 'LineWidth', 1.4, 'DisplayName', PIPE_LABELS(pp));
end
% Identity line
plot(ax8, BETA_RANGE, BETA_RANGE, 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');
% beta_obs markers (horizontal dashed lines)
for pp = 1:N_PIPES
    yline(ax8, betaObs_all(pp), '--', 'Color', [PIPE_CLRS(pp,:), 0.5], ...
        'LineWidth', 0.8, 'HandleVisibility','off');
end
% betaGenStar marker (SG-LMLS, primary pipeline)
if isfinite(betaGenStar_primary)
    xline(ax8, betaGenStar_primary, '-', 'Color', PIPE_CLRS(PP_PRIMARY,:), ...
        'LineWidth', 1.8, 'Label', sprintf('\\beta^*_{gen}=%.3f (SG-LMLS)', betaGenStar_primary), ...
        'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end
legend(ax8, 'Location','northwest', 'FontSize', 7);
xlabel(ax8, '\beta_{gen}'); ylabel(ax8, '\beta_{rec} (median)');
title(ax8, '8 · Forward map — shaped-Xu  N_{reps}=' + string(N_REPS));
xlim(ax8, BETA_RANGE); grid(ax8, 'on');
% Annotation: beta_obs for SG-LMLS (primary)
if isfinite(betaObs_all(PP_PRIMARY))
    text(ax8, BETA_RANGE(2)*0.98, betaObs_all(PP_PRIMARY), ...
        sprintf('\\beta_{obs}=%.3f  (SG-LMLS)', betaObs_all(PP_PRIMARY)), ...
        'Color', PIPE_CLRS(PP_PRIMARY,:), 'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', 'FontSize', 8);
end

fprintf('\nFigure complete.\n');

%% -----------------------------------------------------------------------
%% FIGURE 2: Generated ellipses at representative beta_gen values
%% Shows how beta_gen controls the speed profile on the ellipse.
%% Each panel: PCA template (speed-coloured) + 3 shaped_xu surrogates + real.
%% -----------------------------------------------------------------------
% Base display grid
% Panel order, 2x3 grid:
%   Row 1: [1/12 = low]  [β_gen* exact, green]  [1/3 canonical, red]
%   Row 2: [1/2 = mid]   [2/3 = high]            [3/4 = ceiling]
% β_gen* and 1/3 sit side-by-side so the visual comparison is immediate.
showBetaTargets = [1/12, 1/3, 1/3, 1/2, 2/3, 3/4];   % slot 2 overridden by STAR_SLOT
showIdx = arrayfun(@(b) find(abs(betaGenVec - b) == min(abs(betaGenVec - b)), 1), ...
    showBetaTargets);

% Slot 6 is replaced by the exact beta_gen* panel (uses templateBetaStar, not a grid point)
% showIdx(6) is kept as a fallback if templateBetaStar is unavailable.
USE_EXACT_STAR = isfinite(betaGenStarDisplay) && ~isempty(templateBetaStar);
STAR_SLOT      = 2;   % slot 2: β_gen* exact (green), slot 3: 1/3 canonical (red) — side by side row 1
fprintf('Best-match panel: beta_gen*=%.4f (%s)\n', betaGenStarDisplay, ...
    ternaryStr(USE_EXACT_STAR, 'exact stored', 'nearest grid fallback'));

fig2 = findobj('Type','figure','Tag','vlc_fig2');
if isempty(fig2), fig2 = figure('Tag','vlc_fig2'); else clf(fig2); end
set(fig2, 'Name','Generated Ellipses by beta_gen', ...
    'Position', [120 60 1300 700], 'Color', 'w');
tl2  = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, sprintf('Generated ellipses at representative \\beta_{gen}  —  %s', ...
    strrep(tr.trialID, '_', '\_')), 'FontSize', 12, 'FontWeight', 'bold');

% Real trajectory for reference (centred, in mm, OLA-trimmed window)
xRealRef = x(csOLA:csOLA+M-1)*SIG_TO_MM; xRealRef = xRealRef - mean(xRealRef);
yRealRef = y(csOLA:csOLA+M-1)*SIG_TO_MM; yRealRef = yRealRef - mean(yRealRef);

cmapSpd = turbo(256);

% Global speed range across all shown beta_gen values so low-beta panels
% appear nearly uniform and high-beta panels show the full variation.
spdAll = [];
for k = 1:6
    xTk = templates(:,1,showIdx(k))*SIG_TO_MM;
    yTk = templates(:,2,showIdx(k))*SIG_TO_MM;
    spdAll = [spdAll; hypot(gradient(xTk)*FS, gradient(yTk)*FS)]; %#ok<AGROW>
end
% PATCH: fold in templateBetaStar's own speed range too. The beta_gen* slot
% (STAR_SLOT) does NOT plot one of the six showIdx grid templates above --
% it plots templateBetaStar, an off-grid interpolated template. Without
% this, templateBetaStar's speed values could fall outside [min(spdAll),
% max(spdAll)] and get silently clipped to the colourbar's extreme colour
% by the clim() call applied per-panel below.
if USE_EXACT_STAR
    xTStar = templateBetaStar(:,1)*SIG_TO_MM;
    yTStar = templateBetaStar(:,2)*SIG_TO_MM;
    spdAll = [spdAll; hypot(gradient(xTStar)*FS, gradient(yTStar)*FS)]; %#ok<AGROW>
end
spdGlobalMin = min(spdAll);
spdGlobalMax = max(spdAll);

for k = 1:6
    bi   = showIdx(k);
    bVal = betaGenVec(bi);

    % For the beta_gen* slot: use exact stored template; override bi/bVal
    isStarPanel = (k == STAR_SLOT) && USE_EXACT_STAR;
    if isStarPanel
        bVal = betaGenStarDisplay;
        xT   = templateBetaStar(:,1)*SIG_TO_MM; xT = xT - mean(xT);
        yT   = templateBetaStar(:,2)*SIG_TO_MM; yT = yT - mean(yT);
    else
        % Template in mm, centred
        xT = templates(:,1,bi)*SIG_TO_MM; xT = xT - mean(xT);
        yT = templates(:,2,bi)*SIG_TO_MM; yT = yT - mean(yT);
    end

    % Speed of the template (mm/s) — normalised against GLOBAL range
    vxT  = gradient(xT)*FS; vyT = gradient(yT)*FS;
    spdT = hypot(vxT, vyT);
    spdTN = (spdT - spdGlobalMin) / max(spdGlobalMax - spdGlobalMin, eps);

    axK = nexttile; hold(axK, 'on');

    % Real trajectory (faint reference)
    plot(axK, xRealRef, yRealRef, '-', 'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.7, 'HandleVisibility','off');

    % N_VIS_REPS shaped_xu surrogates — raw noise, no orbit centring
    % Run SG-IRLS regression on the last surrogate to get beta_rec for this panel
    betaRec_panel = NaN;
    for rep = 1:N_VIS_REPS
        [nMj, nMn] = generateLoopClosureNoise_v002('shaped_xu', ...
            resid_maj, resid_min, FS, alphaMaj, alphaMin);
        xS = templates(:,1,bi) + nMj*cos(theta) - nMn*sin(theta);
        yS = templates(:,2,bi) + nMj*sin(theta) + nMn*cos(theta);
        xS = xS*SIG_TO_MM - mean(xS*SIG_TO_MM);
        yS = yS*SIG_TO_MM - mean(yS*SIG_TO_MM);
        plot(axK, xS, yS, '-', 'Color', SUR_CLR, ...
            'LineWidth', 0.9, 'HandleVisibility','off');
        if rep == N_VIS_REPS   % regress on last surrogate (SG-LMLS, primary)
            try
                xSnat = xS / SIG_TO_MM; ySnat = yS / SIG_TO_MM;
                [dxS, dyS] = differentiateKinematicsEBR(xSnat, ySnat, ...
                    DERIV_CFG(2).filterType, DERIV_CFG(2).filterParams, FS);
                vxS = dxS(EDGE_CLIP:end-EDGE_CLIP,2);
                vyS = dyS(EDGE_CLIP:end-EDGE_CLIP,2);
                axS = dxS(EDGE_CLIP:end-EDGE_CLIP,3);
                ayS = dyS(EDGE_CLIP:end-EDGE_CLIP,3);
                spS = hypot(vxS,vyS);
                kpS = curvatureKinematicEBR(vxS,vyS,axS,ayS);
                lm  = [betaObs_primary, betaObs_all(2)];  % LMLS seed
                [betaRec_panel, ~] = regressDataEBR(spS, kpS, REG_CFG(2).type, lm, 0, 0);
            catch; end
        end
    end

    % Template coloured by speed on top (global clim so panels are comparable)
    for seg = 1:M-1
        cidx = max(1, round(spdTN(seg)*255) + 1);
        plot(axK, xT(seg:seg+1), yT(seg:seg+1), '-', ...
            'Color', cmapSpd(cidx,:), 'LineWidth', 1.8, 'HandleVisibility','off');
    end

    axis(axK, 'equal'); grid(axK, 'on');
    xlabel(axK, 'x (mm)'); ylabel(axK, 'y (mm)');
    clim(axK, [spdGlobalMin spdGlobalMax]);

    % Panel title
    if isStarPanel
        titleStr = sprintf('\\beta_{gen}* = %.4f  (stored LC, exact)', bVal);
        titleClr = [0.10 0.50 0.10];
    elseif abs(bVal - 1/3) < 0.02
        titleStr = sprintf('\\beta_{gen} = %.3f  (\\approx 1/3, Lacquaniti 1983)', bVal);
        titleClr = [0.80 0.10 0.10];
    elseif isfinite(betaGenStar_primary) && abs(bVal - betaGenStar_primary) < ...
            (betaGenVec(2)-betaGenVec(1))*1.5
        titleStr = sprintf('\\beta_{gen} = %.3f  (\\beta^*_{gen}, interpolated)', bVal);
        titleClr = [0.10 0.50 0.10];
    else
        titleStr = sprintf('\\beta_{gen} = %.3f', bVal);
        titleClr = [0.1 0.1 0.1];
    end
    title(axK, titleStr, 'Color', titleClr, 'FontSize', 9, 'FontWeight','bold');

    % Annotations: speed ratio, beta_obs real vs simulated
    spdRatio = max(spdT) / max(min(spdT), eps);
    annoStr  = sprintf('speed ratio %.1f:1\n\\beta_{obs} = %.3f  (SG-LMLS, empirical)', ...
        spdRatio, betaObs_primary);
    if isfinite(betaRec_panel)
        annoStr = sprintf('%s\n\\beta_{obs} = %.3f  (SG-LMLS, simulated)', annoStr, betaRec_panel);
        if isfinite(betaObs_primary)
            annoStr = sprintf('%s\n\\Delta = %.3f', annoStr, betaRec_panel - betaObs_primary);
        end
    end
    text(axK, 0.03, 0.97, annoStr, ...
        'Units','normalized', 'VerticalAlignment','top', 'FontSize', 7, ...
        'Color', [0.2 0.2 0.2]);
end

% Single shared colorbar for the whole figure
colormap(fig2, turbo);
cb2 = colorbar;
cb2.Layout.Tile = 'east';
cb2.Label.String = 'template speed (mm/s)';
cb2.Label.FontSize = 9;

% Shared legend via annotation
annotation(fig2, 'textbox', [0.01 0.01 0.98 0.04], ...
    'String', ['Line colour = instantaneous speed of PCA template (turbo: blue=slow, red=fast)   |   ' ...
               'Gold lines = 3 shaped-Xu surrogates (noise centred per orbit — models one-cycle correction bandwidth)   |   ' ...
               'Grey = real trimmed trajectory'], ...
    'EdgeColor','none', 'FontSize', 8, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');

fprintf('Figure 2 complete.\n');

%% -----------------------------------------------------------------------
%% FIGURE 2b: Orbit-centred vs raw — direct comparison at beta_gen ~ 1/3
%% -----------------------------------------------------------------------
bi13cmp  = showIdx(3);
xTcmp    = templates(:,1,bi13cmp)*SIG_TO_MM;
yTcmp    = templates(:,2,bi13cmp)*SIG_TO_MM;
xTcmp    = xTcmp - mean(xTcmp); yTcmp = yTcmp - mean(yTcmp);
N_CMP    = N_VIS_REPS;   % alias — use shared constant

fig2b = findobj('Type','figure','Tag','vlc_fig2b');
if isempty(fig2b), fig2b = figure('Tag','vlc_fig2b'); else clf(fig2b); end
set(fig2b, 'Name','Orbit-centred vs raw comparison', ...
    'Position',[140 50 900 420], 'Color','w');
tl2b  = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl2b, sprintf('Noise injection comparison — \\beta_{gen}=%.3f  (N_{sur}=%d)', ...
    betaGenVec(bi13cmp), N_CMP), 'FontSize', 11, 'FontWeight','bold');

panelTitles = {'Raw noise (no orbit centring)', ...
               'Orbit-centred noise  [VIS ONLY]'};
for mode = 1:2
    axC = nexttile; hold(axC,'on');
    plot(axC, xRealRef, yRealRef, '-', 'Color',[0.75 0.75 0.75], ...
        'LineWidth',0.8, 'DisplayName','real');
    for rep = 1:N_VIS_REPS
        [nMj, nMn] = generateLoopClosureNoise_v002('shaped_xu', ...
            resid_maj, resid_min, FS, alphaMaj, alphaMin);
        if mode == 2                                    % orbit-centred
            nMj = centreByOrbit_lcl(nMj, f0, FS);     % VIS_ONLY
            nMn = centreByOrbit_lcl(nMn, f0, FS);     % VIS_ONLY
        end
        xS = templates(:,1,bi13cmp) + nMj*cos(theta) - nMn*sin(theta);
        yS = templates(:,2,bi13cmp) + nMj*sin(theta) + nMn*cos(theta);
        xS = xS*SIG_TO_MM - mean(xS*SIG_TO_MM);
        yS = yS*SIG_TO_MM - mean(yS*SIG_TO_MM);
        plot(axC, xS, yS, '-', 'Color', SUR_CLR, ...
            'LineWidth', 0.9, 'HandleVisibility','off');
    end
    % Template on top
    plot(axC, xTcmp, yTcmp, 'k-', 'LineWidth',1.6, 'DisplayName','template');
    axis(axC,'equal'); grid(axC,'on');
    xlabel(axC,'x (mm)'); ylabel(axC,'y (mm)');
    title(axC, panelTitles{mode}, 'FontSize',9);
    legend(axC,'Location','best','FontSize',8);
end

fprintf('Figure 2b complete.\n');

%% -----------------------------------------------------------------------
%% NOISE FIDELITY DIAGNOSTIC
%% Does running the full pipeline on synthetic trials recover the real alpha/sigma?
%% Generates N_DIAG synthetic trials at beta_gen ~ 1/3, runs OLA + PCA + IRASA
%% on each, and compares recovered alpha/sigma to the real residual values.
%% -----------------------------------------------------------------------
N_DIAG   = 12;
bi13     = showIdx(3);   % beta_gen closest to 1/3
xTpl13   = templates(:,1,bi13);
yTpl13   = templates(:,2,bi13);

fprintf('\nNoise fidelity diagnostic: IRASA on %d synthetic trials at beta_gen=%.3f\n', ...
    N_DIAG, betaGenVec(bi13));
alphaSyn  = NaN(N_DIAG, 2);   % col1=maj, col2=min
sigmaSyn  = NaN(N_DIAG, 2);

for k = 1:N_DIAG
    [nMj, nMn] = generateLoopClosureNoise_v002('shaped_xu', ...
        resid_maj, resid_min, FS, alphaMaj, alphaMin);
    nMj = centreByOrbit_lcl(nMj, f0, FS);   % VIS_ONLY — not applied in v007
    nMn = centreByOrbit_lcl(nMn, f0, FS);   % VIS_ONLY — not applied in v007
    xS = xTpl13 + nMj*cos(theta) - nMn*sin(theta);
    yS = yTpl13 + nMj*sin(theta) + nMn*cos(theta);

    % Same pipeline: OLA -> PCA rotation -> IRASA
    [xRs, yRs] = templateSubtract_lcl_posonly(xS, yS, FS, f0, N_OLA_H, N_OLA_CW);
    rMaj_s =  xRs*cos(theta) + yRs*sin(theta);
    rMin_s = -xRs*sin(theta) + yRs*cos(theta);

    try
        alphaSyn(k,1) = iraAlphaSigma_v001(rMaj_s, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
        sigmaSyn(k,1) = std(rMaj_s) * SIG_TO_MM;
    catch; end
    try
        alphaSyn(k,2) = iraAlphaSigma_v001(rMin_s, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
        sigmaSyn(k,2) = std(rMin_s) * SIG_TO_MM;
    catch; end
end

% Real values for comparison
sigmaReal_maj = std(resid_maj) * SIG_TO_MM;
sigmaReal_min = std(resid_min) * SIG_TO_MM;

fprintf('              alpha_maj        alpha_min        sigma_maj(mm)   sigma_min(mm)\n');
fprintf('Real:         %.3f            %.3f            %.4f          %.4f\n', ...
    alphaMaj, alphaMin, sigmaReal_maj, sigmaReal_min);
fprintf('Syn mean:     %.3f (sd%.3f)  %.3f (sd%.3f)  %.4f (sd%.4f)  %.4f (sd%.4f)\n', ...
    mean(alphaSyn(:,1),'omitnan'), std(alphaSyn(:,1),'omitnan'), ...
    mean(alphaSyn(:,2),'omitnan'), std(alphaSyn(:,2),'omitnan'), ...
    mean(sigmaSyn(:,1),'omitnan'), std(sigmaSyn(:,1),'omitnan'), ...
    mean(sigmaSyn(:,2),'omitnan'), std(sigmaSyn(:,2),'omitnan'));

% Figure 3: diagnostic scatter
fig3 = findobj('Type','figure','Tag','vlc_fig3');
if isempty(fig3), fig3 = figure('Tag','vlc_fig3'); else clf(fig3); end
set(fig3, 'Name','Noise fidelity diagnostic','Position',[160 80 700 380],'Color','w');
tl3  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl3, sprintf('Noise fidelity: synthetic IRASA recovery  (N=%d, \\beta_{gen}=%.3f)', ...
    N_DIAG, betaGenVec(bi13)), 'FontSize', 11, 'FontWeight','bold');

% Alpha recovery
ax3a = nexttile; hold(ax3a,'on');
xj   = 0.08 * (rand(N_DIAG,1) - 0.5);
scatter(ax3a, ones(N_DIAG,1)+xj, alphaSyn(:,1), 38, CLR_MAJ, 'filled', 'MarkerFaceAlpha',0.7);
scatter(ax3a, 2*ones(N_DIAG,1)+xj, alphaSyn(:,2), 38, CLR_MIN, 'filled', 'MarkerFaceAlpha',0.7);
% Real values as horizontal reference lines
yline(ax3a, alphaMaj, '--', 'Color', CLR_MAJ, 'LineWidth', 1.8, ...
    'Label', sprintf('real \\alpha_{maj}=%.2f', alphaMaj), ...
    'LabelHorizontalAlignment','left');
yline(ax3a, alphaMin, '--', 'Color', CLR_MIN, 'LineWidth', 1.8, ...
    'Label', sprintf('real \\alpha_{min}=%.2f', alphaMin), ...
    'LabelHorizontalAlignment','right');
set(ax3a, 'XTick', [1 2], 'XTickLabel', {'maj','min'}, 'XLim',[0.5 2.5]);
ylabel(ax3a, '\alpha (IRASA)'); grid(ax3a,'on');
title(ax3a, '\alpha: synthetic vs real');

% Sigma recovery
ax3b = nexttile; hold(ax3b,'on');
scatter(ax3b, ones(N_DIAG,1)+xj, sigmaSyn(:,1), 38, CLR_MAJ, 'filled', 'MarkerFaceAlpha',0.7);
scatter(ax3b, 2*ones(N_DIAG,1)+xj, sigmaSyn(:,2), 38, CLR_MIN, 'filled', 'MarkerFaceAlpha',0.7);
yline(ax3b, sigmaReal_maj, '--', 'Color', CLR_MAJ, 'LineWidth', 1.8, ...
    'Label', sprintf('real \\sigma_{maj}=%.3fmm', sigmaReal_maj), ...
    'LabelHorizontalAlignment','left');
yline(ax3b, sigmaReal_min, '--', 'Color', CLR_MIN, 'LineWidth', 1.8, ...
    'Label', sprintf('real \\sigma_{min}=%.3fmm', sigmaReal_min), ...
    'LabelHorizontalAlignment','right');
set(ax3b, 'XTick', [1 2], 'XTickLabel', {'maj','min'}, 'XLim',[0.5 2.5]);
ylabel(ax3b, '\sigma (mm)'); grid(ax3b,'on');
title(ax3b, '\sigma: synthetic vs real (sigma-inflation check)');

fprintf('Figure 3 complete.\n');

%% -----------------------------------------------------------------------
%% LOCAL FUNCTIONS
%% -----------------------------------------------------------------------

function n = centreByOrbit_lcl(n, f0, fs)
% Remove the per-orbit mean from a noise signal.
% VISUALISATION ONLY — this function is NOT called in runLoopClosureFftnoise_v007.
% The forward map uses uncentred surrogates consistent with the preregistered analysis.
% Here it is applied solely to produce figures whose trajectories are visually
% comparable to real drawings. Orbit centring models the assumption that online
% motor correction resets positional drift once per movement cycle.
    n          = n(:);
    orbitLen   = round(fs / f0);
    nOrbits    = floor(numel(n) / orbitLen);
    for k = 1:nOrbits
        idx    = (k-1)*orbitLen + (1:orbitLen);
        n(idx) = n(idx) - mean(n(idx));
    end
end

function [resX, resY] = templateSubtract_lcl_posonly(x, y, fs, f0, nH, nCW)
% Lightweight wrapper: returns only residuals, no cs. Used in the diagnostic.
    [resX, resY, ~, ~, ~] = templateSubtract_lcl(x, y, fs, f0, nH, nCW);
end
function out = ternaryStr(cond, a, b)
    if cond, out = a; else, out = b; end
end

function f0 = estimateF0_lcl(x, y, fs)
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x  = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band)), f0 = f0x;
    else, f0 = f0y;
    end
end

function [resX, resY, fitX, fitY, cs] = templateSubtract_lcl(x, y, fs, f0, nH, nCW)
    N   = numel(x);
    win = round(nCW/f0*fs); hop = round(win/2);
    win = min(win, N); win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    s = 1;
    while s + win - 1 <= N
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e).*han); bY = Dw \ (y(s:e).*han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX).*han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY).*han;
        ww(s:e)   = ww(s:e) + han; s = s + hop;
    end
    ok = ww > 0;
    resX(ok) = resX(ok)./ww(ok); resY(ok) = resY(ok)./ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl; ce = min(max(N-cl, cs+100), N);
    fitX = x(cs:ce) - resX(cs:ce); fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce); resY = resY(cs:ce);
end
