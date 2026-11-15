function [i0, i1, info] = endToEndSegment_v003(S, opts)
%ENDTOENDSEGMENT_V003 endToEndSegment_v002's fairShare rule for several channels at once (one segment for all of them).
%   [i0, i1, info] = endToEndSegment_v003(S) with S an N-by-C matrix (e.g. [resid_maj, resid_min]) returns the longest
%   segment S(i0:i1, :) whose wrap terms are within Tol of the segment's own mean terms in EVERY channel:
%       cJump(c) = (S(i0,c) - S(i1,c))^2 / mean(diff(S(i0:i1,c)).^2)                                      <= Tol
%       cSlip(c) = mean of the two wrapped second differences squared / mean(diff(S(i0:i1,c), 2).^2)        <= Tol
%   Among the segments of that length, the one with the smallest max over channels of max(cJump, cSlip).
%   One segment for both axes, because the synthetic trajectory that the noise is added to has one length.
%   With one column it is identical to endToEndSegment_v002(s, Rule="fairShare") (checked in checkEndToEndSegmentV003_v001).
%   Rationale, the bound it gives on the derivative noise, and its source (Nakamura, Small & Hirata 2006, Sec. V C, applied
%   to length): see endToEndSegment_v002 and docs/EndMismatch_FairShareRule_FullLength_20261114_v001.md.
%   info: Lfrac, cJump (1-by-C), cSlip (1-by-C), crit (the accepted max), Tol.
arguments
    S (:,:) double {mustBeFinite}
    opts.Tol (1,1) double {mustBePositive} = 1
    opts.MinLength (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.MinLength, 5)} = 5
end
[N, C] = size(S);
if N < opts.MinLength, error('endToEndSegment:short', 'series of %d samples is shorter than MinLength %d.', N, opts.MinLength); end
e1 = [zeros(1, C); cumsum(diff(S).^2)];  e2 = [zeros(1, C); cumsum(diff(S, 2).^2)];
i0 = NaN;  i1 = NaN;  crit = NaN;
for L = N:-1:opts.MinLength
    a = (1:N - L + 1)';  b = a + L - 1;
    m1 = (e1(b, :) - e1(a, :)) / (L - 1);  m2 = (e2(b - 1, :) - e2(a, :)) / (L - 2);
    cj = (S(a, :) - S(b, :)).^2 ./ m1;
    cs = ((S(a + 1, :) - 2 * S(a, :) + S(b, :)).^2 + (S(a, :) - 2 * S(b, :) + S(b - 1, :)).^2) / 2 ./ m2;
    c = max([cj, cs], [], 2);  c(~isfinite(c)) = Inf;
    [cMin, k] = min(c);
    if cMin <= opts.Tol, i0 = a(k);  i1 = b(k);  crit = cMin;  break, end
end
if isnan(i0), error('endToEndSegment:noSegment', 'no segment of %d or more samples has wrap terms within Tol %g in all %d channels.', opts.MinLength, opts.Tol, C); end
G = S(i0:i1, :);
info.Lfrac = (i1 - i0 + 1) / N;  info.Tol = opts.Tol;  info.crit = crit;
info.cJump = (G(1, :) - G(end, :)).^2 ./ mean(diff(G).^2, 1);
info.cSlip = ((G(2, :) - 2 * G(1, :) + G(end, :)).^2 + (G(1, :) - 2 * G(end, :) + G(end - 1, :)).^2) / 2 ./ mean(diff(G, 2).^2, 1);
end
