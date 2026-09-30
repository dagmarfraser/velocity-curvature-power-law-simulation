function checkFraserConstellationMedian_v001(matFile)
% CHECKFRASERCONSTELLATIONMEDIAN_V001  Provenance of the Fraser 0.3292 (claude.md dataset table, Table 6b).
%   Shows that the headline is median(betaGenStarMed), the per-trial NaN-omitting
%   median of the six pipelines, and prints each pipeline's own median for contrast.
%   Read-only. Fails loud if the file layout differs from what is asserted.
if nargin < 1
    here = fileparts(mfilename('fullpath'));
    matFile = fullfile(here, 'loopClosureResults_Fraser_all_shaped_xu_v008.mat');
end
S = load(matFile, 'results', 'pipelineLabels');
r = S.results;  labels = string(S.pipelineLabels);
assert(numel(labels) == 6, 'Expected 6 pipelines, got %d', numel(labels));
B = nan(numel(r), 6);
for i = 1:numel(r)
    v = r(i).betaGenStar(:)';
    assert(numel(v) == 6, 'Trial %d: betaGenStar has %d entries, expected 6', i, numel(v));
    B(i, :) = v;
end
bm = [r.betaGenStarMed]';
d  = max(abs(bm - median(B, 2, 'omitnan')), [], 'omitnan');
assert(d < 1e-12, 'betaGenStarMed is not the row-wise NaN-omitting median (max diff %.3g)', d);
fprintf('File: %s\n', matFile);
fprintf('N=%d  median(betaGenStarMed)=%.4f  nNaN=%d\n', numel(r), median(bm, 'omitnan'), nnz(isnan(bm)));
fprintf('max|betaGenStarMed - rowwise nanmedian| = %.3g\n', d);
for k = 1:6
    fprintf('%-10s median=%.4f  nNaN=%d\n', labels(k), median(B(:, k), 'omitnan'), nnz(isnan(B(:, k))));
end
end
