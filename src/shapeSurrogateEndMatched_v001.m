function [u, ref] = shapeSurrogateEndMatched_v001(r, carrierFn, method, opts)
%SHAPESURROGATEENDMATCHED_V001 Magnitude/phase-swap surrogate of r with a published remedy for end-point mismatch.
%   [u, ref] = shapeSurrogateEndMatched_v001(r, carrierFn, method) returns a surrogate u carrying the Fourier amplitudes
%   of r on the phases of carrierFn(n) (an n-by-1 real series, e.g. Xu noise), and ref, the series u is a surrogate OF
%   (r itself, or the chosen sub-segment for "endToEnd"). No rescaling and no mean removal: the caller does both.
%   method:
%     "plain"    |FFT(r)| on the carrier's phases. The old shapeXu construction; no remedy (baseline).
%     "mirror"   symmetrisation: swap on [r; flipud(r)], which has no end-to-start jump, with a carrier of length 2M;
%                first M samples kept (Nakamura, Small & Hirata 2006, Phys Rev E 74:026205, as used by Guarin Lopez
%                et al. 2010, arXiv 1008.1804).
%     "tfts"     truncated Fourier transform surrogate: r keeps its own phases below CutoffHz and takes the carrier's
%                phases at and above it (Nakamura et al. 2006).
%     "detrend"  least-squares polynomial trend of order TrendOrder removed, swap on the remainder, trend added back
%                (Lucio, Valdes & Rodriguez 2012, Phys Rev E 85:056202; TFTD in TimeseriesSurrogates.jl uses order 1).
%     "endToEnd" swap on the sub-segment with the smallest end mismatch (Schreiber & Schmitz 2000, Physica D 142:346;
%                endToEndSegment_v001). u is SHORTER than r: numel(u) = numel(ref).
%   Candidate for the fractal-noise toolbox; written 2026-11-14 for endMismatchSensitivity_v001. Nakamura 2006 and Lucio
%   2012 had not been read in full when this was written: check each branch against its paper before citing.
arguments
    r (:,1) double {mustBeFinite}
    carrierFn (1,1) function_handle
    method (1,1) string {mustBeMember(method, ["plain", "mirror", "tfts", "detrend", "endToEnd"])}
    opts.Fs (1,1) double {mustBePositive} = 1
    opts.CutoffHz (1,1) double {mustBePositive} = 1
    opts.TrendOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.Segment double = []                    % [i0 i1] for "endToEnd"; empty = search with endToEndSegment_v001
end
M = numel(r);  ref = r;  swap = @(a, c) real(ifft(abs(fft(a)) .* exp(1i * angle(fft(c)))));
switch method
    case "plain"
        u = swap(r, carrier(carrierFn, M));
    case "mirror"
        e = swap([r; flipud(r)], carrier(carrierFn, 2 * M));  u = e(1:M);
    case "tfts"
        R = fft(r);  ph = angle(fft(carrier(carrierFn, M)));
        f = min(0:M - 1, M:-1:1)' * opts.Fs / M;  keep = f < opts.CutoffHz;       % symmetric in +/- frequency
        if all(keep), error('shapeSurrogate:cutoff', 'CutoffHz %g keeps every phase of r (Nyquist %g).', opts.CutoffHz, opts.Fs / 2); end
        ph(keep) = angle(R(keep));  u = real(ifft(abs(R) .* exp(1i * ph)));
    case "detrend"
        t = (0:M - 1)' / (M - 1);  p = polyval(polyfit(t, r, opts.TrendOrder), t);
        u = swap(r - p, carrier(carrierFn, M)) + p;
    case "endToEnd"
        sg = opts.Segment;
        if isempty(sg), [sg(1), sg(2)] = endToEndSegment_v001(r); end
        if numel(sg) ~= 2 || sg(1) < 1 || sg(2) > M || sg(2) - sg(1) < 3, error('shapeSurrogate:segment', 'Segment must be [i0 i1] inside 1..%d.', M); end
        ref = r(sg(1):sg(2));  u = swap(ref, carrier(carrierFn, numel(ref)));
end
end

function c = carrier(fn, n)
c = fn(n);
if ~isreal(c) || numel(c) ~= n || any(~isfinite(c)), error('shapeSurrogate:carrier', 'carrierFn(%d) must return %d finite real samples.', n, n); end
c = c(:);
end
