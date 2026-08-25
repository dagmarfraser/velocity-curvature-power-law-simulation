function out = fresSRpowerLaw_v001(frk, Freq0, opts)
% FRESSRPOWERLAW_V001  Huh & Sejnowski (2015) spectral power-law regression.
%
% Applies the angular-frequency bandpass pipeline to Frenet-Serret
% kinematics and returns two OLS beta estimates in one pass:
%   betaNoInt : OLS forced through origin  (canonical Huh method)
%   betaInt   : OLS with intercept
% Both are returned as positive values (~1/3) consistent with the
% regressDataEBR sign convention: v = yGain * k^(-beta).
%
% SYNTAX
%   out = fresSRpowerLaw_v001(frk, Freq0)
%   out = fresSRpowerLaw_v001(frk, Freq0, ThetaStart=pi/2, LPCutoff=25)
%
% INPUTS
%   frk     struct from frenetKinematics_v001
%           Required fields: v, k, theta, pos, iValid
%   Freq0   shape angular frequency [cycles per 2pi revolution]
%           e.g. 2 = standard ellipse, 4/3 = figure-eight, 4 = complex shape
%           Must be > 0 (spiral Freq0=0 not supported here)
%
% PARAMETERS (name-value)
%   ThetaStart (1x1) skip first N rad of accumulated angle (default pi/2)
%   LPCutoff   (1x1) lowpass cutoff [cycles/rev]           (default 25)
%
% OUTPUTS  struct:
%   betaNoInt   (1x1)  beta, OLS no intercept  (positive, ~1/3)
%   betaCI      (1x2)  95% CI [lo hi] for betaNoInt  (positive)
%   betaInt     (1x1)  beta, OLS with intercept       (positive, ~1/3)
%   yGain       (1x1)  speed-curvature gain = exp(intercept)
%   beta0       (1x1)  Huh theoretical prediction |f0(Freq0)| (positive)
%   powerRatio  (1x1)  sqrt(power_v / power_k) at Freq0
%   logKV_lp    (Nx2)  lowpass log-signals   [log|k|, log v]  on theta_clean grid
%   logKV_gauss (Nx2)  bandpass log-signals  [log|k|, log v]  on theta_clean grid
%   theta_clean (Nx1)  stitched continuous theta for valid samples
%   traj_xy     (Nx2)  position for valid samples
%   nu          (Mx1)  angular frequency axis [cycles/rev] on spline grid
%   nSegs       (1x1)  number of valid segments stitched
%   nSamples    (1x1)  total samples in stitched signal
%
% ALGORITHM
%   1.  Apply iValid AND theta > ThetaStart; find runs spanning > pi rad
%   2.  Trim pi/4 from each run end; stitch with continuous theta
%   3.  Spline-resample to Nsp = 2*NN+1 uniform theta points; build nu axis
%   4.  Lowpass |nu| < LPCutoff applied to [k, v] spline
%   5.  log(|kv_lp|) -> back-interpolate to theta_clean (mean-centred)
%   6.  Bandpass Gaussian at Freq0 [Paper version, H&S 2015]:
%         BPF = exp(-5 * Freq0 * log2(|nu|/Freq0)^2)
%       applied to lowpass log-signals; back-interpolate (mean-centred)
%   7.  OLS no-intercept:   regress(log_v_gauss, log_k_gauss)
%       OLS with-intercept: fitlm(log_k_gauss, log_v_gauss)
%   8.  Power ratio and theoretical beta0 = |f0(Freq0)|
%
% NOTE ON SIGNS
%   regress returns slope b = -beta (negative).  Both betas stored positive.
%   Theoretical: f0(nu) = -2/3*(1+nu^2/2)/(1+nu^2+nu^4/15)  [negative];
%   beta0 = |f0(Freq0)|.
%
% DEPENDENCIES  spline_fitting, splinefit
%
% REFERENCE
%   Huh D & Sejnowski TJ (2015) Spectrum of power laws for curved hand
%   movements. PNAS 112(29): E3950-E3958.
%
% SEE ALSO  frenetKinematics_v001, spline_fitting, regressDataEBR
%
% PROVENANCE
%   Distilled from Analysis_subfunction_Redux_v12.m (Huh / Cook Lab / Fraser).
%   7th toolchain (Frenet-Serret spectral) for PowerLawSimulationPreReg.
%
% AUTHOR  D.S. Fraser, CHBH, University of Birmingham, 2026.

arguments
    frk    (1,1) struct
    Freq0  (1,1) double {mustBePositive}
    opts.ThetaStart (1,1) double = pi/2
    opts.LPCutoff   (1,1) double = 25
end

%% Unpack ------------------------------------------------------------------
theta  = frk.theta;
v      = frk.v;
k      = frk.k;
pos    = frk.pos;
iValid = frk.iValid;

%% Combined gate: iValid + initiation-phase skip --------------------------
iGate = iValid & (theta > opts.ThetaStart);

%% Find contiguous runs; keep those spanning > pi rad ---------------------
runBeg  = find([iGate(1); diff(iGate) > 0]);
runEnd  = find([diff(iGate) < 0; iGate(end)]);
runSpan = theta(runEnd) - theta(runBeg);

keep    = runSpan > pi;
runBeg  = runBeg(keep);
runEnd  = runEnd(keep);

if isempty(runBeg)
    error('fresSRpowerLaw_v001:noValidSegments', '%s', ...
        'No valid segments with theta span > pi. Check iValid / ThetaStart.');
end

%% Stitch segments into continuous theta ----------------------------------
iidx        = [];
thetaStitch = [];
thetaTail   = 0;

for jj = 1:length(runBeg)
    seg  = runBeg(jj):runEnd(jj);
    iBeg = find(theta(seg) > theta(runBeg(jj)) + pi/4, 1, 'first');
    iEnd = find(theta(seg) < theta(runEnd(jj)) - pi/4, 1, 'last');
    if isempty(iBeg) || isempty(iEnd) || iBeg >= iEnd
        continue
    end
    seg0 = seg(iBeg):seg(iEnd);
    iidx = [iidx, seg0]; %#ok<AGROW>

    tSeg    = theta(seg0);
    tShift  = -tSeg(1) + thetaTail + (tSeg(2) - tSeg(1));
    thetaStitch = [thetaStitch; tSeg + tShift]; %#ok<AGROW>
    thetaTail   = thetaStitch(end);
end

if isempty(iidx)
    error('fresSRpowerLaw_v001:noSegmentsAfterTrim', '%s', ...
        'No segments survived pi/4 edge trimming.');
end

thetaClean = thetaStitch;           % N x 1, monotone
kv         = [k(iidx), v(iidx)];   % N x 2  (signed k, speed)
traj_xy    = pos(iidx, :);          % N x 2

nSamples = length(iidx);
NN       = round(nSamples / 2);
Nsp      = 2*NN + 1;

%% Uniform theta grid and angular frequency axis --------------------------
thetaSpline = linspace(thetaClean(1), thetaClean(end), Nsp)';
dTheta      = abs(diff(thetaSpline([1 end])));
nu          = (2*pi / dTheta) .* [0:NN, -NN:-1]'; % cycles per revolution

%% Spline-resample kv to uniform grid -------------------------------------
kvSpline = spline_fitting(thetaClean, thetaSpline, kv, NN);
% kvSpline: Nsp x 2

%% Lowpass filter in nu domain --------------------------------------------
%% Log-transform first, then lowpass in log domain -------------------------
% Applying lowpass to linear k,v then logging (log∘LP) introduces Jensen's
% inequality bias because log is concave: log(LP(k)) ≠ LP(log(k)).
% Correct order: log first, then filter in log domain for both LP and Gauss.
logKvSpline = log(abs(kvSpline));           % Nsp x 2  — log before any filtering

lpFilt    = double(abs(nu) < opts.LPCutoff);
lpFilt_if = ifft(lpFilt);
lpFilt_if = lpFilt_if([NN+1:end, 1:NN]);

logKvLpSpline = [conv(logKvSpline(:,1), lpFilt_if, 'same'), ...
                 conv(logKvSpline(:,2), lpFilt_if, 'same')];  % Nsp x 2

%% Back-interpolate lowpass log-signals to theta_clean (mean-centred) -----
logKvLp = spline_fitting(thetaSpline, thetaClean, logKvLpSpline, NN, true);
% logKvLp: N x 2  [log|k|_lp, log v_lp]

%% OLS on lowpass signals only (no Gaussian bandpass) ---------------------
logKlp = logKvLp(:,1);
logVlp = logKvLp(:,2);

[b_lp, bInt_lp]  = regress(logVlp, logKlp);
betaNoInt_lp     = -b_lp;
betaCI_lp        = sort(-bInt_lp);

mdl_lp           = fitlm(logKlp, logVlp, 'Intercept', true);
coef_lp          = mdl_lp.Coefficients.Estimate;
betaInt_lp       = -coef_lp(2);
yGain_lp         = exp(coef_lp(1));

%% Bandpass Gaussian at Freq0 (Paper version) -----------------------------
bpFilt    = exp(-5 .* Freq0 .* (log2(abs(nu) ./ Freq0)).^2); % Nsp x 1
bpFilt    = bpFilt ./ max(bpFilt);
bpFilt_if = ifft(bpFilt);
bpFilt_if = bpFilt_if([NN+1:end, 1:NN]);

logKvGaussSpline = [conv(logKvLpSpline(:,1), bpFilt_if, 'same'), ...
                    conv(logKvLpSpline(:,2), bpFilt_if, 'same')]; % Nsp x 2

if any(any(isnan(logKvGaussSpline)))
    error('fresSRpowerLaw_v001:gaussNaN', '%s', ...
        sprintf('Bandpass filter produced NaN at Freq0=%.4g.', Freq0));
end

%% Back-interpolate Gauss signals to theta_clean (mean-centred) -----------
logKvGauss = spline_fitting(thetaSpline, thetaClean, logKvGaussSpline, NN, true);
% logKvGauss: N x 2  [log|k|_gauss, log v_gauss]

logKg = logKvGauss(:,1);
logVg = logKvGauss(:,2);

%% OLS no intercept (Huh canonical) --------------------------------------
[b, bInt] = regress(logVg, logKg);
betaNoInt = -b;                     % positive (~1/3)
betaCI    = sort(-bInt);            % [lo hi], positive

%% OLS with intercept -----------------------------------------------------
mdl   = fitlm(logKg, logVg, 'Intercept', true);
coef  = mdl.Coefficients.Estimate;
betaInt = -coef(2);                 % positive
yGain   = exp(coef(1));

%% Power ratio at Freq0 ---------------------------------------------------
fftLogKv  = fft(logKvLpSpline) ./ Nsp;       % Nsp x 2
bpFilt0   = exp(-5 .* Freq0 .* (log2(abs(nu) ./ Freq0)).^2); % un-normalised

pwrKvFreq = sqrt(sum(abs(fftLogKv .* bpFilt0).^2, 1) ./ Nsp); % 1 x 2
powerRatio = pwrKvFreq(2) ./ pwrKvFreq(1);   % v / k at Freq0

%% Theoretical beta (Huh eq. 7) ------------------------------------------
f0fun = @(x) -2/3 .* (1 + x.^2./2) ./ (1 + x.^2 + x.^4./15);
beta0 = abs(f0fun(Freq0));          % positive

%% Assemble output --------------------------------------------------------
out.betaNoInt   = betaNoInt;
out.betaCI      = betaCI;
out.betaInt     = betaInt;
out.yGain       = yGain;
out.betaNoInt_lp = betaNoInt_lp;
out.betaCI_lp    = betaCI_lp;
out.betaInt_lp   = betaInt_lp;
out.yGain_lp     = yGain_lp;
out.beta0       = beta0;
out.powerRatio  = powerRatio;
out.logKV_lp    = logKvLp;
out.logKV_gauss = logKvGauss;
out.theta_clean = thetaClean;
out.traj_xy     = traj_xy;
out.nu          = nu;
out.nSegs       = length(runBeg);
out.nSamples    = nSamples;

end
