function a = dfaAlpha_v001(x, order)
%DFAALPHA_V001 DFA of a given detrending order, returned as a spectral exponent alpha = 2*Hfluct - 1.
%   MATLAB twin of dfaN in fractal-noise stage1/run_stage1_v004.R: profile = cumsum(x - mean), 15 log-spaced box
%   sizes from 10 to N/10, non-overlapping boxes, polynomial detrend of the given order per box, OLS of log F on
%   log s. Order m can only reach alpha < 2m + 1 (Kiyono 2015).
x = x(:);  y = cumsum(x - mean(x));  n = numel(y);
if n < 100, error('dfaAlpha:short', 'Series too short for DFA (N = %d).', n); end
s = unique(round(exp(linspace(log(10), log(floor(n/10)), 15))));
F = zeros(size(s));
for i = 1:numel(s)
    w = s(i);  nb = floor(n / w);  Y = reshape(y(1:nb*w), w, nb);
    t = (1:w)' / w;  P = t .^ (0:order);  [Q, ~] = qr(P, 0);
    R = Y - Q * (Q' * Y);  F(i) = sqrt(mean(R(:).^2));
end
p = polyfit(log(s), log(F), 1);  a = 2 * p(1) - 1;
end
