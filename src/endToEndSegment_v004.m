function [i0, i1, info] = endToEndSegment_v004(S, opts)
%ENDTOENDSEGMENT_V004 fairShare segments for several channels with one common LENGTH (starts may differ per channel).
%   [i0, i1, info] = endToEndSegment_v004(S) with S an N-by-C matrix returns 1-by-C start and end indices. All segments have
%   the same length L, the longest length at which every channel has at least one segment whose wrap terms are within Tol
%   of the segment's own mean terms (the fairShare test of endToEndSegment_v002, per channel); each channel takes its own
%   best segment of that length (smallest max(cJump, cSlip)).
%   Why not one shared segment (v003, Mode "sameStart"): on the real residuals that demand gave a median of only 38-40% of
%   the trial, and no segment at all for about half of the Hickman trials (checkSegmentBetaObs_v001, 2026-11-14). The noise
%   only needs one LENGTH, because the synthetic trajectory has one length; each axis's noise is drawn from its own segment.
%   Mode "sameStart" reproduces endToEndSegment_v003 (i0, i1 then equal in every channel).
%   info: L, Lfrac, cJump (1-by-C), cSlip (1-by-C), mode, Tol.
arguments
    S (:,:) double {mustBeFinite}
    opts.Tol (1,1) double {mustBePositive} = 1
    opts.MinLength (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.MinLength, 5)} = 5
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ["commonLength", "sameStart"])} = "commonLength"
end
[N, C] = size(S);
if N < opts.MinLength, error('endToEndSegment:short', 'series of %d samples is shorter than MinLength %d.', N, opts.MinLength); end
e1 = [zeros(1, C); cumsum(diff(S).^2)];  e2 = [zeros(1, C); cumsum(diff(S, 2).^2)];
i0 = NaN(1, C);  i1 = NaN(1, C);
for L = N:-1:opts.MinLength
    a = (1:N - L + 1)';  b = a + L - 1;
    m1 = (e1(b, :) - e1(a, :)) / (L - 1);  m2 = (e2(b - 1, :) - e2(a, :)) / (L - 2);
    cj = (S(a, :) - S(b, :)).^2 ./ m1;
    cs = ((S(a + 1, :) - 2 * S(a, :) + S(b, :)).^2 + (S(a, :) - 2 * S(b, :) + S(b - 1, :)).^2) / 2 ./ m2;
    c = max(cj, cs);  c(~isfinite(c)) = Inf;                       % per start, per channel
    if opts.Mode == "sameStart"
        [cMin, k] = min(max(c, [], 2));
        if cMin <= opts.Tol, i0(:) = a(k);  i1(:) = b(k);  break, end
    else
        [cMin, k] = min(c, [], 1);
        if all(cMin <= opts.Tol), i0 = a(k)';  i1 = b(k)';  break, end
    end
end
if any(isnan(i0)), error('endToEndSegment:noSegment', 'no length of %d or more samples has a segment within Tol %g in all %d channels (%s).', opts.MinLength, opts.Tol, C, opts.Mode); end
L = i1(1) - i0(1) + 1;  info.L = L;  info.Lfrac = L / N;  info.mode = opts.Mode;  info.Tol = opts.Tol;
info.cJump = NaN(1, C);  info.cSlip = NaN(1, C);
for ch = 1:C
    G = S(i0(ch):i1(ch), ch);
    info.cJump(ch) = (G(1) - G(end))^2 / mean(diff(G).^2);
    info.cSlip(ch) = ((G(2) - 2 * G(1) + G(end))^2 + (G(1) - 2 * G(end) + G(end - 1))^2) / 2 / mean(diff(G, 2).^2);
end
end
