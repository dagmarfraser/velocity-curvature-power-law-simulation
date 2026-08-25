function [nMaj, nMin] = generateLoopClosureNoise_v002(model, srcMaj, srcMin, FS, alphaMaj, alphaMin)
% generateLoopClosureNoise_v002  Surrogate noise for loop-closure forward map.
%
% Iterates from v001 by adding a third surrogate model, 'shaped_xu', behind
% the same dispatcher interface. v001 callers are unaffected: 'fftnoise' and
% 'xu' branches are bit-identical to v001.
%
% MODEL = 'fftnoise':
%   srcMaj/srcMin are the FFT of the PCA-axis residuals (fftMaj/fftMin in the
%   runner). Phase-scrambled, amplitude preserved. alpha args ignored.
%   Reproduces the empirical amplitude spectral density exactly (spectral
%   breaks, tremor peaks) BUT alpha is NOT preserved above ~3 (Finding #67):
%   at Cook/Hickman residual alpha~5 the surrogate has alpha~2.7. The
%   diff-diff regularisation conjecture (that double differentiation washes
%   this out before the curvature step) was FALSIFIED by
%   checkDiffDiffAlphaAttenuation_v001: the alpha gap is invariant under
%   differentiation (a fixed linear filter shifts real and surrogate alpha
%   identically), persisting at ~2.65 from position to acceleration.
%
% MODEL = 'xu':
%   srcMaj/srcMin are the PCA-axis residual TIME SERIES (real column vectors).
%   Xu (2019) coloured noise at the supplied per-axis alpha, scaled to the
%   residual std. Reproduces alpha faithfully to ~6 via genuine fractional
%   structure (GGM operator, time domain) BUT imposes a single power-law
%   slope: no spectral break, no tremor peak. The missing break shows up as
%   a ~0.4 alpha offset at the derivative levels (checkDiffDiffAlpha) and as
%   forward-map over-compression (median beta_gen* 0.13 vs fftnoise 0.41).
%   alphaMaj/alphaMin REQUIRED.
%
% MODEL = 'shaped_xu'  (NEW in v002):
%   srcMaj/srcMin are the PCA-axis residual TIME SERIES (as 'xu').
%   Construction: wear the EMPIRICAL amplitude spectrum on Xu's FRACTIONAL
%   PHASE.
%       F_shaped = |FFT(resid)| .* exp(1i*angle(FFT(Xu_base)))
%   where Xu_base = generateCustomNoise_v003(M, alpha, std(resid), FS).
%   Rationale (checkFftnoiseAlphaPreservation_v002 + checkDiffDiffAlpha):
%   the ONLY lever on IRASA alpha is PHASE, not magnitude. fftnoise uses
%   random phase (alpha collapses to ~2.7); Xu's phase carries the high-alpha
%   fractional scaling. shaped_xu keeps the empirical magnitude (break,
%   shoulder, tremor — like fftnoise) but replaces random phase with Xu's
%   fractional phase (high alpha — like xu). It is the minimal surrogate that
%   targets BOTH the empirical fine structure AND the correct noise colour.
%   Conjugate symmetry is automatic (empirical magnitude symmetric, Xu phase
%   antisymmetric), so the IFFT is real. alphaMaj/alphaMin REQUIRED.
%
%   NOTE: whether Xu's phase survives the magnitude swap to still deliver
%   high alpha is an empirical question, gated by an alpha-preservation check
%   BEFORE shaped_xu is used for any beta_gen* run. Do not assume validity.
%
% The three models span the design space:
%   fftnoise  : empirical fine structure YES, correct alpha NO (caps ~2.7)
%   xu        : empirical fine structure NO,  correct alpha YES (single slope)
%   shaped_xu : empirical fine structure YES, correct alpha YES (target)
%
% INPUTS:
%   model     - 'fftnoise' | 'xu' | 'shaped_xu'
%   srcMaj    - fftnoise: fft(resid_maj) [complex Mx1]
%               xu/shaped_xu: resid_maj [real Mx1]
%   srcMin    - fftnoise: fft(resid_min) [complex Mx1]
%               xu/shaped_xu: resid_min [real Mx1]
%   FS        - sampling rate (Hz) — used by xu/shaped_xu
%   alphaMaj  - xu/shaped_xu: IRASA alpha of major-axis residual
%   alphaMin  - xu/shaped_xu: IRASA alpha of minor-axis residual
%
% OUTPUTS:
%   nMaj, nMin - real surrogate column vectors, length M, one per PCA axis
%
% Fraser, D.S. (2026)

    arguments
        model    (1,:) char {mustBeMember(model, {'fftnoise','xu','shaped_xu','bootstrap'})}
        srcMaj   (:,1)
        srcMin   (:,1)
        FS       (1,1) double {mustBePositive} = 100
        alphaMaj (1,1) double = NaN
        alphaMin (1,1) double = NaN
    end

    switch model
        case 'fftnoise'
            % srcMaj/srcMin are already FFTs; phase-scramble, preserve amplitude
            nMaj = fftnoise(srcMaj);
            nMin = fftnoise(srcMin);

        case 'xu'
            if ~isfinite(alphaMaj) || ~isfinite(alphaMin)
                error('generateLoopClosureNoise_v002:alphaRequired', ...
                    '%s', 'xu model requires finite alphaMaj and alphaMin.');
            end
            M        = numel(srcMaj);
            ALPHA_LO = -2.0; ALPHA_HI = 5.95;  % generateCustomNoise_v003 valid range
            if alphaMaj < ALPHA_LO || alphaMaj > ALPHA_HI || ...
               alphaMin < ALPHA_LO || alphaMin > ALPHA_HI
                warning('generateLoopClosureNoise_v002:xuAlphaOutOfRange', '%s', ...
                    sprintf('xu: alphaMaj=%.3f alphaMin=%.3f outside [%.1f,%.2f] -- rep rejected.', ...
                    alphaMaj, alphaMin, ALPHA_LO, ALPHA_HI));
                nMaj = NaN(M, 1);
                nMin = NaN(M, 1);
                return
            end
            sigMaj  = std(srcMaj, 0, 1);   % residual std per axis = injection sigma
            sigMin  = std(srcMin, 0, 1);
            nMaj = generateCustomNoise_v003(M, alphaMaj, sigMaj, FS);
            nMin = generateCustomNoise_v003(M, alphaMin, sigMin, FS);

        case 'bootstrap'
            % Same circular shift on both axes — preserves cross-axis phase coherence.
            % srcMaj/srcMin are residual time series (real column vectors; same
            % convention as 'xu'/'shaped_xu'). alpha args not used.
            %
            % A circular shift is a linear phase rotation in the frequency domain:
            %   FFT(circshift(r, k))_w = FFT(r)_w * exp(i*w*k)
            % Applying the SAME k to both axes leaves all cross-axis phase DIFFERENCES
            % invariant. Every per-axis spectral statistic (alpha, sigma, spectral break)
            % is preserved exactly. No IRASA required; no alpha range check needed.
            M     = numel(srcMaj);
            shift = randi(M);
            nMaj  = circshift(srcMaj(:), shift);
            nMin  = circshift(srcMin(:), shift);

        case 'shaped_xu'
            if ~isfinite(alphaMaj) || ~isfinite(alphaMin)
                error('generateLoopClosureNoise_v002:alphaRequired', ...
                    '%s', 'shaped_xu model requires finite alphaMaj and alphaMin.');
            end
            nMaj = shapeXu_local(srcMaj, FS, alphaMaj);
            nMin = shapeXu_local(srcMin, FS, alphaMin);
    end
end

%% ===== LOCAL FUNCTIONS =====================================================

function n = shapeXu_local(resid, FS, alpha)
% Empirical amplitude spectrum on Xu fractional phase.
%   F = |FFT(resid)| .* exp(1i*angle(FFT(Xu_base)))
% Conjugate symmetry is inherited automatically: |FFT(resid)| is symmetric
% and angle(FFT(Xu_base)) is antisymmetric for real inputs, so F is a valid
% real-signal spectrum. real() guards residual numerical asymmetry.
    resid = resid(:);
    M     = numel(resid);
    sig   = std(resid, 0, 1);
    if sig == 0
        n = zeros(M, 1);
        return
    end
    Eemp  = abs(fft(resid));                              % empirical magnitude (break/tremor)
    if alpha < -2 || alpha > 5.95                          % generateCustomNoise_v003 valid range [-2,6]
        warning('generateLoopClosureNoise:alphaOutOfRange', '%s', ...
            sprintf('alpha %.3f out of [-2, 5.95] for Xu generator — rep rejected (NaN returned).', alpha));
        n = NaN(M, 1);
        return
    end
    base  = generateCustomNoise_v003(M, alpha, sig, FS);  % Xu fractional-phase carrier
    Fb    = fft(base);
    Fout  = Eemp .* exp(1i*angle(Fb));                   % empirical magnitude, Xu phase
    n     = real(ifft(Fout));
    sN    = std(n, 0, 1);
    if sN > 0
        n = n * (sig / sN);                              % restore injection sigma
    end
end
