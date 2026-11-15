function [i0, i1, info] = endToEndSegment_v002(s, opts)
%ENDTOENDSEGMENT_V002 Longest sub-segment of s whose wrap-around step looks like an ordinary step.
%   [i0, i1, info] = endToEndSegment_v002(s) returns the segment s(i0:i1) chosen by Rule:
%
%   Rule "fairShare" (default; data-driven, no fixed floor). A Fourier surrogate treats the segment as periodic, so the
%   step from s(i1) back to s(i0) becomes one more first difference, and two more second differences, of the noise. With
%   random phases that energy is spread over the whole surrogate, so it adds to the velocity and acceleration noise.
%   The rule accepts a segment when each wrap term is no larger than Tol times the mean of the segment's own terms:
%       cJump = (s(i0) - s(i1))^2 / mean(diff(seg).^2)                                   <= Tol
%       cSlip = mean([(s(i0+1) - 2 s(i0) + s(i1))^2, (s(i0) - 2 s(i1) + s(i1-1))^2])
%               / mean(diff(seg, 2).^2)                                                  <= Tol
%   With Tol = 1 the wrap adds at most 1/(L-1) to the first-difference energy and 2/(L-2) to the second-difference
%   energy of a segment of length L (under 1% for the trials here), which is the bias it can put into the derivative
%   noise. Lengths are searched from the full series downwards (the stepwise stopping strategy of Nakamura, Small &
%   Hirata 2006, Phys Rev E 74:026205, Sec. V C, applied to segment length instead of the randomised band), and the
%   longest length with an accepted segment is kept; among its segments the one with the smallest max(cJump, cSlip).
%   The resulting length is a property of the data, reported in info.Lfrac.
%
%   Rule "minFraction" reproduces endToEndSegment_v001: segments of at least MinFraction of the series, minimising
%   JumpWeight*gJump + (1 - JumpWeight)*gSlip with the measures of Schreiber & Schmitz 2000, Eqs. 16-17 (checked
%   against the paper on 2026-11-14: docs/FormulaCheck_SchreiberSchmitzEq16-17_20261114_v001.md).
%
%   info: Lfrac, cJump, cSlip (fairShare terms of the chosen segment), gJump, gSlip (Eqs. 16-17), rule.
arguments
    s (:,1) double {mustBeFinite}
    opts.Rule (1,1) string {mustBeMember(opts.Rule, ["fairShare", "minFraction"])} = "fairShare"
    opts.Tol (1,1) double {mustBePositive} = 1
    opts.MinFraction (1,1) double {mustBeInRange(opts.MinFraction, 0, 1)} = 0.75
    opts.JumpWeight (1,1) double {mustBeInRange(opts.JumpWeight, 0, 1)} = 0.5
    opts.MinLength (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.MinLength, 5)} = 5
end
N = numel(s);
if N < opts.MinLength, error('endToEndSegment:short', 'series of %d samples is shorter than MinLength %d.', N, opts.MinLength); end
c1 = [0; cumsum(s)];  c2 = [0; cumsum(s.^2)];
d1 = diff(s);  e1 = [0; cumsum(d1.^2)];  d2 = diff(s, 2);  e2 = [0; cumsum(d2.^2)];
i0 = NaN;  i1 = NaN;
if opts.Rule == "minFraction"
    Lmin = max(4, ceil(opts.MinFraction * N));  w = opts.JumpWeight;  best = Inf;
    for L = Lmin:N
        a = (1:N - L + 1)';  b = a + L - 1;
        den = (c2(b + 1) - c2(a)) - (c1(b + 1) - c1(a)).^2 / L;
        gj = (s(a) - s(b)).^2 ./ den;  gs = ((s(a + 1) - s(a)) - (s(b) - s(b - 1))).^2 ./ den;
        [g, k] = min(w * gj + (1 - w) * gs);
        if g < best, best = g;  i0 = a(k);  i1 = b(k); end
    end
    if ~isfinite(best), error('endToEndSegment:degenerate', 'no segment has a finite mismatch (constant series?).'); end
else
    for L = N:-1:opts.MinLength
        a = (1:N - L + 1)';  b = a + L - 1;
        m1 = (e1(b) - e1(a)) / (L - 1);                       % mean of d1(a:b-1).^2
        m2 = (e2(b - 1) - e2(a)) / (L - 2);                   % mean of d2(a:b-2).^2
        cj = (s(a) - s(b)).^2 ./ m1;
        cs = ((s(a + 1) - 2 * s(a) + s(b)).^2 + (s(a) - 2 * s(b) + s(b - 1)).^2) / 2 ./ m2;
        c = max(cj, cs);  c(~isfinite(c)) = Inf;
        [cMin, k] = min(c);
        if cMin <= opts.Tol, i0 = a(k);  i1 = b(k);  break, end
    end
    if isnan(i0), error('endToEndSegment:noSegment', 'no segment of %d or more samples has wrap terms within Tol %g.', opts.MinLength, opts.Tol); end
end
seg = s(i0:i1);  L = numel(seg);  v = sum((seg - mean(seg)).^2);
info.rule = opts.Rule;  info.Lfrac = L / N;
info.cJump = (seg(1) - seg(end))^2 / mean(diff(seg).^2);
info.cSlip = ((seg(2) - 2 * seg(1) + seg(end))^2 + (seg(1) - 2 * seg(end) + seg(end - 1))^2) / 2 / mean(diff(seg, 2).^2);
info.gJump = (seg(1) - seg(end))^2 / v;  info.gSlip = ((seg(2) - seg(1)) - (seg(end) - seg(end - 1)))^2 / v;
end
