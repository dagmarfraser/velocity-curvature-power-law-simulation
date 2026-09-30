function out = checkRegisteredNoiseExpectations_v001(matFile)
% Prereg 8.1/6.1 noise-specific SEM expectations vs the grid (coordTable). Read-only.
% Prints every number quoted in manuscript Section 3c's rebuilt table.
if nargin < 1, matFile = fullfile('src','perCoordinateSEM_v2_001.mat'); end
assert(isfile(matFile), 'checkRegisteredNoiseExpectations:noFile', 'Not found: %s', matFile);
S = load(matFile, 'coordTable'); T = S.coordTable;
crit = 0.03/2.77; pipes = categories(T.pipeline)';
sig = unique(T.sigma)'; sig = sig(sig >= 0.1);
rows = cell(0,1);
for a = 0:3
  for s = sig
    for k = 1:numel(pipes)
      m = T.alpha==a & abs(T.sigma-s)<1e-9 & T.pipeline==pipes{k} & isfinite(T.sem);
      x = T.sem(m);
      rows{end+1,1} = table(a, s, string(pipes{k}), median(x), 100*mean(x<crit), 100*mean(x>=crit), 100*mean(x>=0.020), numel(x), ...
        'VariableNames',{'alpha','sigma','pipeline','medSEM','pctAdequate','pctGEcrit','pctGE020','n'}); %#ok<AGROW>
    end
  end
end
out = vertcat(rows{:});
fprintf('criterion MDC/2.77 = %.5f\n', crit);
fprintf('\n[alpha=0] median SEM range and %%coords >= crit / >= 0.020, sigma 0.1 vs 2 mm\n');
for s = [0.1 2]
  r = out(out.alpha==0 & abs(out.sigma-s)<1e-9,:);
  fprintf(' sigma %.1f: medSEM %.4f-%.4f; pct>=crit %.1f-%.1f\n', s, min(r.medSEM), max(r.medSEM), min(r.pctGEcrit), max(r.pctGEcrit));
end
r = out(out.alpha==0 & abs(out.sigma-2)<1e-9,:);
for k = 1:height(r), fprintf('  %-9s med %.4f  >=crit %.1f%%  >=0.020 %.1f%%\n', r.pipeline(k), r.medSEM(k), r.pctGEcrit(k), r.pctGE020(k)); end
fprintf('\nSigma (grid points) at which each pipeline median first reaches crit:\n');
for a = 1:3
  for k = 1:numel(pipes)
    r = out(out.alpha==a & out.pipeline==pipes{k},:); i = find(r.medSEM >= crit, 1);
    if isempty(i), fprintf(' alpha %d %-9s never (max sigma %g)\n', a, pipes{k}, max(r.sigma));
    else, fprintf(' alpha %d %-9s first >= crit at sigma %g (previous grid point %g)\n', a, pipes{k}, r.sigma(i), r.sigma(max(i-1,1))); end
  end
end
fprintf('\nBWFD-OLS median SEM at sigma = 2 mm by alpha:\n');
for a = 0:3, r = out(out.alpha==a & abs(out.sigma-2)<1e-9 & out.pipeline=="BWFD-OLS",:); fprintf(' alpha %d: %.4f\n', a, r.medSEM); end
r = out(out.alpha==3 & abs(out.sigma-20)<1e-9,:);
fprintf('\nalpha=3, sigma=20: pct adequate %.1f-%.1f\n', min(r.pctAdequate), max(r.pctAdequate));
end
