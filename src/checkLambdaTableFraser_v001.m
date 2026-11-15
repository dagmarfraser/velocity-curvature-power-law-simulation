function R = checkLambdaTableFraser_v001(tableFile)
% CHECKLAMBDATABLEFRASER_V001  Fraser rows of a checkAdditiveVsRetimed_v002 lambda table: trials evaluable and excluded (by status),
%   per-subject median lambda, and Spearman(per-subject median lambda, per-subject median f0), the figures behind Part 3's
%   generative-model paragraph ("could not be evaluated (N of Fraser's M trials ...); across Fraser subjects, lambda ...").
%   Anchor: on the 2026-09-22 20:51 table (N=45) it must give 1,152 excluded of 2,829 and Spearman .32 (p = .04), Finding #222.
%   R = checkLambdaTableFraser_v001(<path to a checkAdditiveVsRetimed_v002_*.mat>)   Fraser, D.S. (2026)  v001
T = load(tableFile, 'T').T;  t = T(T.dataset == "Fraser", :);
R.n = height(t);  R.nOk = nnz(t.status == "ok");  R.nTooSlow = nnz(endsWith(t.status, "olaWindowExceedsTrace"));
R.nUnsep = nnz(t.status == "anchorsUnseparated");
m = isfinite(t.lambda);  g = findgroups(t.subjectID(m));
s = splitapply(@median, t.lambda(m), g);  f = splitapply(@median, t.f0(m), g);
[R.rho, R.p] = corr(s, f, 'Type', 'Spearman');  R.nSubj = numel(s);  R.subjMedLambda = median(s);
fprintf('%s\n  Fraser %d trials: %d too slow for the OLA window, %d anchors unseparated, %d evaluable; %d subjects, median of subject medians %.3f\n', ...
    tableFile, R.n, R.nTooSlow, R.nUnsep, R.nOk, R.nSubj, R.subjMedLambda);
fprintf('  Spearman(subject lambda, subject f0) = %.2f (p = %.3f)\n', R.rho, R.p);
if R.n == 2829 && (R.nTooSlow ~= 1152 || round(R.rho, 2) ~= 0.32), error('lamTab:anchor', '%s', 'N=45 table does not reproduce Finding #222 (1,152; .32)'), end
end
