function [i0, i1, info] = endToEndSegment_v005(S, opts)
%ENDTOENDSEGMENT_V005 endToEndSegment_v004 (commonLength) with an optional derivative-energy budget instead of a fixed Tol.
%   [i0, i1, info] = endToEndSegment_v005(S, Budget=0.01) accepts, for each channel, a segment of length L whose wrap terms
%   add at most Budget of the segment's own first- and second-difference energy:
%       (S(i0) - S(i1))^2                         <= Budget * sum(diff(seg).^2)       i.e. cJump <= Budget * (L - 1)
%       sum of the two wrapped 2nd differences^2   <= Budget * sum(diff(seg, 2).^2)    i.e. cSlip <= Budget * (L - 2) / 2
%   which is exactly the bias bound endToEndSegment_v002 quoted for Tol 1 (1/(L-1), 2/(L-2)), now held fixed across lengths
%   instead of tightening as L grows. Lengths are searched from the full series down (Nakamura, Small & Hirata 2006), the
%   longest length at which every channel has an accepted segment is kept, each channel takes its own best segment
%   (smallest max of cJump/limJ, cSlip/limS). With Budget empty or NaN the rule is v004's (Tol, default 1), so v005 then
%   reproduces v004 exactly.
%   Why (2026-11-15): with the template tail excluded, the fixed Tol 1 left about a quarter of Hickman trials without a
%   segment of the runner's 280-sample minimum; Tol 1 is far stricter than its own 1% justification for long segments.
%   info: L, Lfrac, cJump, cSlip (1-by-C), rule ("tol" or "budget"), Tol, Budget.
arguments
    S (:,:) double {mustBeFinite}
    opts.Tol (1,1) double {mustBePositive} = 1
    opts.Budget double = []
    opts.MinLength (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.MinLength, 5)} = 5
end
useB = ~isempty(opts.Budget) && isfinite(opts.Budget);
if useB, mustBePositive(opts.Budget); end
[N, C] = size(S);
if N < opts.MinLength, error('endToEndSegment:short', 'series of %d samples is shorter than MinLength %d.', N, opts.MinLength); end
e1 = [zeros(1, C); cumsum(diff(S).^2)];  e2 = [zeros(1, C); cumsum(diff(S, 2).^2)];
i0 = NaN(1, C);  i1 = NaN(1, C);
for L = N:-1:opts.MinLength
    a = (1:N - L + 1)';  b = a + L - 1;
    m1 = (e1(b, :) - e1(a, :)) / (L - 1);  m2 = (e2(b - 1, :) - e2(a, :)) / (L - 2);
    cj = (S(a, :) - S(b, :)).^2 ./ m1;
    cs = ((S(a + 1, :) - 2 * S(a, :) + S(b, :)).^2 + (S(a, :) - 2 * S(b, :) + S(b - 1, :)).^2) / 2 ./ m2;
    if useB, limJ = opts.Budget * (L - 1);  limS = opts.Budget * (L - 2) / 2; else, limJ = opts.Tol;  limS = opts.Tol; end
    c = max(cj / limJ, cs / limS);  c(~isfinite(c)) = Inf;     % <= 1 means accepted
    [cMin, k] = min(c, [], 1);
    if all(cMin <= 1), i0 = a(k)';  i1 = b(k)';  break, end
end
if any(isnan(i0)), error('endToEndSegment:noSegment', 'no length of %d or more samples has a segment within the %s in all %d channels.', opts.MinLength, string(ifelse(useB, "budget", "tolerance")), C); end
L = i1(1) - i0(1) + 1;  info.L = L;  info.Lfrac = L / N;  info.Tol = opts.Tol;  info.Budget = opts.Budget;
info.rule = string(ifelse(useB, "budget", "tol"));  info.cJump = NaN(1, C);  info.cSlip = NaN(1, C);
for ch = 1:C
    G = S(i0(ch):i1(ch), ch);
    info.cJump(ch) = (G(1) - G(end))^2 / mean(diff(G).^2);
    info.cSlip(ch) = ((G(2) - 2 * G(1) + G(end))^2 + (G(1) - 2 * G(end) + G(end - 1))^2) / 2 / mean(diff(G, 2).^2);
end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
