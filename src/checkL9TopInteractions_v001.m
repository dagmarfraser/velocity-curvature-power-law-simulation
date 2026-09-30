function tops = checkL9TopInteractions_v001(matFile, nTop)
% Largest fixed-effect interaction terms of the registered L9 LMM. Read-only.
% regressionType_4 = LMLS, _5 = IRLS (OLS = 3, reference).
if nargin < 1, matFile = fullfile('results','L9ModelInspect_v002.mat'); end
if nargin < 2, nTop = 8; end
assert(isfile(matFile), 'checkL9TopInteractions:noFile', 'Not found: %s', matFile);
S = load(matFile, 'out'); F = S.out.fixed;
nm = string(F.name); est = double(F.est); I = find(contains(nm, ':'));
[~, ix] = sort(abs(est(I)), 'descend'); I = I(ix(1:nTop));
tops = table(nm(I), est(I), 'VariableNames', {'term','est'});
disp(tops);
end
