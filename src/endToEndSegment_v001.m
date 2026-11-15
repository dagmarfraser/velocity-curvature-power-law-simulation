function [i0, i1, gJump, gSlip] = endToEndSegment_v001(s, opts)
%ENDTOENDSEGMENT_V001 Sub-segment of s with the smallest end-point mismatch (Schreiber & Schmitz 2000, Sec. 4.5).
%   [i0, i1, gJump, gSlip] = endToEndSegment_v001(s) searches every segment s(i0:i1) of at least MinFraction of the
%   series and returns the one minimising  JumpWeight*gJump + (1 - JumpWeight)*gSlip, with
%     gJump = (s(i0) - s(i1))^2 / sum((seg - mean(seg)).^2)
%     gSlip = ((s(i0+1) - s(i0)) - (s(i1) - s(i1-1)))^2 / sum((seg - mean(seg)).^2)
%   as in the TISEAN routine endtoend (equal weights by default). The two formulas were taken from a summary of the
%   paper on 2026-11-14 and must be checked against its Eqs. 16-17 before the function is cited.
arguments
    s (:,1) double {mustBeFinite}
    opts.MinFraction (1,1) double {mustBeInRange(opts.MinFraction, 0, 1)} = 0.75
    opts.JumpWeight (1,1) double {mustBeInRange(opts.JumpWeight, 0, 1)} = 0.5
end
N = numel(s);  Lmin = max(4, ceil(opts.MinFraction * N));  w = opts.JumpWeight;
if N < Lmin, error('endToEndSegment:short', 'series of %d samples is shorter than the minimum segment (%d).', N, Lmin); end
c1 = [0; cumsum(s)];  c2 = [0; cumsum(s.^2)];  best = Inf;
for L = Lmin:N
    a = (1:N - L + 1)';  b = a + L - 1;
    den = (c2(b + 1) - c2(a)) - (c1(b + 1) - c1(a)).^2 / L;
    gj = (s(a) - s(b)).^2 ./ den;  gs = ((s(a + 1) - s(a)) - (s(b) - s(b - 1))).^2 ./ den;
    [g, k] = min(w * gj + (1 - w) * gs);
    if g < best, best = g;  i0 = a(k);  i1 = b(k);  gJump = gj(k);  gSlip = gs(k); end
end
if ~isfinite(best), error('endToEndSegment:degenerate', 'no segment has a finite mismatch (constant series?).'); end
end
