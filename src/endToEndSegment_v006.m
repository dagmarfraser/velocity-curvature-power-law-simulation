function [i0, i1, info] = endToEndSegment_v006(S, opts)
%ENDTOENDSEGMENT_V006 endToEndSegment_v005 plus a parameter-free "detectability" rule (2026-11-15).
%   Rules (one length L shared by all channels, each channel its own segment, lengths searched from the full series down,
%   the longest accepted length kept; Nakamura, Small & Hirata 2006):
%     Rule "tol"     wrap terms <= Tol times the segment's mean terms (v002-v004; Tol 1 by default)
%     Rule "budget"  the wrap adds at most Budget of the segment's first- and second-difference energy (v005)
%     Rule "detect"  the wrap step is no larger than the largest step expected among the segment's own steps, and each of the
%                    two wrapped second differences no larger than the largest expected among its own second differences:
%                        (S(i0) - S(i1))^2 / mean(diff(seg).^2)            <= q(L - 1)
%                        max(wrapped 2nd differences^2) / mean(diff(seg,2).^2) <= q(L - 2)
%                    with q(n) = chi2inv(1 - 1/n, 1), the level a squared standardised Gaussian step exceeds once in n steps
%                    (about 2 ln n). No number is chosen: the wrap must not stand out from the segment's own steps. Its
%                    energy share is then at most q(L-1)/(L-1) (about 4% at L 280, 1.4% at 1000, 0.5% at 3000).
%   Defaults reproduce v005: Rule "tol" with Tol 1 unless Budget is given (then "budget"). v006 with Rule "tol" equals
%   endToEndSegment_v004/v005 exactly.
%   Why: Dagmar (2026-11-15) asked for the budget to be grounded without reference to beta_gen*; "detect" needs no parameter.
%   info: L, Lfrac, cJump, cSlip, cSlipMax (1-by-C), rule, Tol, Budget, qJ, qS.
arguments
    S (:,:) double {mustBeFinite}
    opts.Rule (1,1) string {mustBeMember(opts.Rule, ["auto", "tol", "budget", "detect"])} = "auto"
    opts.Tol (1,1) double {mustBePositive} = 1
    opts.Budget double = []
    opts.MinLength (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.MinLength, 5)} = 5
end
rule = opts.Rule;
if rule == "auto", if ~isempty(opts.Budget) && isfinite(opts.Budget), rule = "budget"; else, rule = "tol"; end, end
if rule == "budget" && (isempty(opts.Budget) || ~isfinite(opts.Budget) || opts.Budget <= 0), error('endToEndSegment:budget', 'Rule "budget" needs a positive Budget.'); end
[N, C] = size(S);
if N < opts.MinLength, error('endToEndSegment:short', 'series of %d samples is shorter than MinLength %d.', N, opts.MinLength); end
e1 = [zeros(1, C); cumsum(diff(S).^2)];  e2 = [zeros(1, C); cumsum(diff(S, 2).^2)];
i0 = NaN(1, C);  i1 = NaN(1, C);  qJ = NaN;  qS = NaN;
for L = N:-1:opts.MinLength
    a = (1:N - L + 1)';  b = a + L - 1;
    m1 = (e1(b, :) - e1(a, :)) / (L - 1);  m2 = (e2(b - 1, :) - e2(a, :)) / (L - 2);
    cj = (S(a, :) - S(b, :)).^2 ./ m1;
    w1 = (S(a + 1, :) - 2 * S(a, :) + S(b, :)).^2;  w2 = (S(a, :) - 2 * S(b, :) + S(b - 1, :)).^2;
    switch rule
        case "tol",    c = max(cj / opts.Tol, (w1 + w2) / 2 ./ m2 / opts.Tol);
        case "budget", c = max(cj / (opts.Budget * (L - 1)), (w1 + w2) / 2 ./ m2 / (opts.Budget * (L - 2) / 2));
        case "detect", qJ = chi2inv(1 - 1 / (L - 1), 1);  qS = chi2inv(1 - 1 / (L - 2), 1);  c = max(cj / qJ, max(w1, w2) ./ m2 / qS);
    end
    c(~isfinite(c)) = Inf;
    [cMin, k] = min(c, [], 1);
    if all(cMin <= 1), i0 = a(k)';  i1 = b(k)';  break, end
end
if any(isnan(i0)), error('endToEndSegment:noSegment', 'no length of %d or more samples has an accepted segment in all %d channels (rule %s).', opts.MinLength, C, rule); end
L = i1(1) - i0(1) + 1;  info.L = L;  info.Lfrac = L / N;  info.rule = rule;  info.Tol = opts.Tol;  info.Budget = opts.Budget;  info.qJ = qJ;  info.qS = qS;
info.cJump = NaN(1, C);  info.cSlip = NaN(1, C);  info.cSlipMax = NaN(1, C);
for ch = 1:C
    G = S(i0(ch):i1(ch), ch);  m2c = mean(diff(G, 2).^2);
    v1 = (G(2) - 2 * G(1) + G(end))^2;  v2 = (G(1) - 2 * G(end) + G(end - 1))^2;
    info.cJump(ch) = (G(1) - G(end))^2 / mean(diff(G).^2);  info.cSlip(ch) = (v1 + v2) / 2 / m2c;  info.cSlipMax(ch) = max(v1, v2) / m2c;
end
end
