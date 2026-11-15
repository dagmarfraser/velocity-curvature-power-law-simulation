function checkEndToEndSegmentV003_v001()
%CHECKENDTOENDSEGMENTV003_V001 Acceptance checks for endToEndSegment_v003 and endMatchedNoise_v001. Errors on any failure.
%   1. One column: v003 picks the same segment as endToEndSegment_v002(Rule="fairShare") on 200 random series
%      (random walks and their integrals, lengths 300-3000); where neither finds a segment, both must raise the same error.
%   2. Two columns: the chosen segment satisfies Tol in both channels, and no longer segment does (brute force, 20 short series).
%   3. endMatchedNoise_v001: length, std per axis equal to the segment's, accepted draws within Tol, NaN when nothing can be
%      accepted (Tol 1e-9), identical output for an identical substream.
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
rs = RandStream("threefry4x64_20", "Seed", 20261114);  nNone = 0;
for k = 1:200
    n = randi(rs, [300 3000]);  s = cumsum(randn(rs, n, 1));  if mod(k, 2) == 0, s = cumsum(s); end
    e2 = "";  e3 = "";  a2 = NaN;  b2 = NaN;  a3 = NaN;  b3 = NaN;
    try, [a2, b2] = endToEndSegment_v002(s, Rule="fairShare"); catch ME, e2 = string(ME.identifier); end
    try, [a3, b3] = endToEndSegment_v003(s); catch ME, e3 = string(ME.identifier); end
    if e2 ~= e3 || ~isequaln([a2 b2], [a3 b3]), error('check:oneColumn', 'series %d: v002 [%d %d] "%s", v003 [%d %d] "%s"', k, a2, b2, e2, a3, b3, e3); end
    nNone = nNone + (e2 ~= "");
end
fprintf('1. one column: v003 = v002 fairShare on 200 series (%d with no acceptable segment in both)\n', nNone);
for k = 1:20
    n = randi(rs, [60 140]);  S = cumsum(randn(rs, n, 2));
    [a, b, info] = endToEndSegment_v003(S);  L = b - a + 1;
    if any(info.cJump > 1) || any(info.cSlip > 1), error('check:tol', 'series %d: chosen segment exceeds Tol', k); end
    for LL = n:-1:L + 1
        for i = 1:n - LL + 1
            G = S(i:i + LL - 1, :);
            cj = (G(1, :) - G(end, :)).^2 ./ mean(diff(G).^2, 1);
            cs = ((G(2, :) - 2 * G(1, :) + G(end, :)).^2 + (G(1, :) - 2 * G(end, :) + G(end - 1, :)).^2) / 2 ./ mean(diff(G, 2).^2, 1);
            if all(cj <= 1) && all(cs <= 1), error('check:longer', 'series %d: a longer segment (%d > %d) also satisfies Tol', k, LL, L); end
        end
    end
end
fprintf('2. two columns: chosen segment within Tol in both, none longer qualifies (20 series, brute force)\n');
FS = 133;  seg = cumsum(randn(rs, 1500, 2));  [g0, g1] = endToEndSegment_v003(seg);  seg = seg(g0:g1, :);  seg = seg - mean(seg);  Lg = size(seg, 1);
aMaj = alphaPWLW_v001(seg(:, 1), FS, [1 20]);  aMin = alphaPWLW_v001(seg(:, 2), FS, [1 20]);
s0 = rs;  s0.Substream = 7;  RandStream.setGlobalStream(s0);  [n1, m1, i1] = endMatchedNoise_v001(seg(:, 1), seg(:, 2), FS, aMaj, aMin);
s0 = RandStream("threefry4x64_20", "Seed", 20261114);  s0.Substream = 7;  RandStream.setGlobalStream(s0);  [n2, m2] = endMatchedNoise_v001(seg(:, 1), seg(:, 2), FS, aMaj, aMin);
if ~isequaln(n1, n2) || ~isequaln(m1, m2), error('check:seed', 'same substream gave different noise'); end
if numel(n1) ~= Lg || numel(m1) ~= Lg, error('check:length', 'wrong length'); end
if ~i1.acceptedMaj && ~i1.acceptedMin, error('check:accept', 'no axis accepted a draw on an end-matched test segment'); end
if i1.acceptedMaj && abs(std(n1) - std(seg(:, 1))) > 1e-9 * std(seg(:, 1)), error('check:std', 'major-axis std not matched'); end
if i1.acceptedMaj && abs(i1.dAlphaMaj) > 0.3, error('check:alpha', 'accepted draw outside Tol'); end
if i1.acceptedMin && (abs(i1.dAlphaMin) > 0.3 || abs(std(m1) - std(seg(:, 2))) > 1e-9 * std(seg(:, 2))), error('check:minor', 'minor-axis std or alpha wrong'); end
[n3, m3, i3] = endMatchedNoise_v001(seg(:, 1), seg(:, 2), FS, aMaj, aMin, Tol=1e-9, MaxDraw=3);
if ~all(isnan(n3)) || ~all(isnan(m3)) || i3.acceptedMaj || i3.acceptedMin, error('check:failLoud', 'unaccepted draw not returned as NaN'); end
fprintf('3. endMatchedNoise_v001: length, std, alpha tolerance, NaN on no acceptance, reproducible by substream (accepted maj %d, min %d)\n', i1.acceptedMaj, i1.acceptedMin);
fprintf('ALL CHECKS PASSED\n');
end
