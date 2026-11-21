% compareK1vsK4_v001 (2026-11-21): side-by-side K=4 vs K=1 for every Table 6b cell and every scalar number of the variance decomposition output.
% K=4 sources are the EXISTING outputs (old tree results/buildTable6bV028_v001_v028.mat and src/loopClosureVarDecomp_v014_gated_all6_v028.mat);
% K=1 sources are the v015/v002 outputs in BRM/outputs_K1. Difference = K1 - K4. "crosses003" = |difference| >= 0.03 (the registered MDC); it is only
% meaningful where betaUnits is true (quantities on the beta scale), and is computed for all rows regardless.
B2B = fullfile(getenv('HOME'),'Library','CloudStorage','Dropbox','Brain2Bee');
old = fullfile(B2B,'PowerLawSimulationPreReg');  od = fullfile(B2B,'PowerLaw_Repro_2026','BRM','outputs_K1');
outCsv = fullfile(od, 'K1_vs_K4_Table6b_varDecomp_20261121_v001.csv');
if isfile(outCsv), error('compareK1vsK4:exists', '%s exists; not overwritten', outCsv); end
T4 = load(fullfile(old,'results','buildTable6bV028_v001_v028.mat')).T;
T1 = load(fullfile(od,'buildTable6bV028_v002_K1.mat')).T;
V4 = load(fullfile(old,'src','loopClosureVarDecomp_v014_gated_all6_v028.mat'));
V1 = load(fullfile(od,'loopClosureVarDecomp_v015_all6_v028_nH1_K1.mat'));
assert(V1.meta.nHarm == 1 && V1.meta.dataSource == "v028");
beta = ["flat" "cluster" "ciLo" "ciHi" "mean" "sd" "clusterMinusOneThird"];   % Table 6b columns on the beta scale
R = struct('section', {}, 'item', {}, 'quantity', {}, 'K4', {}, 'K1', {}, 'text4', {}, 'text1', {}, 'betaUnits', {});
% ---- Table 6b ----
assert(isequal(T4.dataset, T1.dataset));
num = ["flat" "cluster" "ciLo" "ciHi" "mean" "sd" "Nfin" "nFiniteMed" "clusterMinusOneThird" "ciIncludesOneThird"];
for i = 1:height(T4)
    for q = num
        R(end+1) = struct('section', "Table6b", 'item', T4.dataset(i), 'quantity', q, 'K4', double(T4.(q)(i)), 'K1', double(T1.(q)(i)), 'text4', "", 'text1', "", 'betaUnits', any(beta == q)); %#ok<SAGROW>
    end
    R(end+1) = struct('section', "Table6b", 'item', T4.dataset(i), 'quantity', "bandStatus", 'K4', NaN, 'K1', NaN, 'text4', T4.bandStatus(i), 'text1', T1.bandStatus(i), 'betaUnits', false); %#ok<SAGROW>
end
% ---- variance decomposition: every scalar/short vector numeric field of the per-dataset entries and the syntheses ----
betaRx = '^(pointEst|se.*|ci.*|median|medians|mPool.*|sePool.*|range|maxCiUpper|poolMarginBeta.*|marginBeta|sig.*|tau2)$';
ents = ["cluster_median" "cluster_primary" "flat_median" "flat_primary" "test_median" "test_primary" "dec_median" "dec_primary"];
lab4 = string({V4.out.label});  lab1 = string({V1.out.label});  assert(isequal(lab4, lab1));
for i = 1:numel(V4.out)
    for e = ents
        R = flat_local(R, "varDecomp", lab4(i), e, V4.out(i).(e), V1.out(i).(e), betaRx);
    end
end
for e = ents(1:2)
    R = flat_local(R, "varDecomp", "Hickman (combined)", e, V4.hickmanCombined.(e), V1.hickmanCombined.(e), betaRx);
end
for s = ["synthMed" "synthPrim" "synthMed_k3" "synthPrim_k3" "synthMed_k4" "synthPrim_k4" "synthMed_k6" "synthPrim_k6"]
    R = flat_local(R, "varDecomp", s, "synth", V4.(s), V1.(s), betaRx);
end
T = struct2table(R);
T.difference = T.K1 - T.K4;  T.absDifference = abs(T.difference);  T.crosses003 = T.absDifference >= 0.03;
T.textChanged = T.text4 ~= T.text1;
T = T(:, {'section','item','quantity','betaUnits','K4','K1','difference','absDifference','crosses003','text4','text1','textChanged'});
writetable(T, outCsv);
fprintf('rows %d; Table6b %d; varDecomp %d; beta-unit rows crossing 0.03: %d (Table6b %d, varDecomp %d)\n', height(T), sum(T.section=="Table6b"), sum(T.section=="varDecomp"), ...
    sum(T.crosses003 & T.betaUnits), sum(T.crosses003 & T.betaUnits & T.section=="Table6b"), sum(T.crosses003 & T.betaUnits & T.section=="varDecomp"));
fprintf('NaN differences (one side NaN): %d\n', sum(isnan(T.difference) & isnan(T.K4) ~= isnan(T.K1)));
disp(T(T.crosses003 & T.betaUnits, {'section','item','quantity','K4','K1','difference'}));

function R = flat_local(R, section, item, grp, a, b, betaRx)
fn = fieldnames(a);
for j = 1:numel(fn)
    f = fn{j};  x = a.(f);  y = b.(f);
    if isstruct(x), R = flat_local(R, section, item, grp + "." + f, x, y, betaRx); continue, end
    if ~(isnumeric(x) || islogical(x)), continue, end
    if numel(x) > 12 || numel(y) > 12, continue, end   % bootstrap draws and long vectors are not tabulated
    if numel(x) ~= numel(y), error('compare:size', '%s %s %s size differs (%d vs %d)', item, grp, f, numel(x), numel(y)); end
    for k = 1:numel(x)
        q = grp + "." + f;  if numel(x) > 1, q = q + "[" + k + "]"; end
        R(end+1) = struct('section', section, 'item', item, 'quantity', q, 'K4', double(x(k)), 'K1', double(y(k)), 'text4', "", 'text1', "", ...
            'betaUnits', ~isempty(regexp(f, betaRx, 'once')) && isempty(regexp(f, '^(frac|p[A-Z]|poolP|poolZ|poolT)|Excludes$', 'once'))); %#ok<AGROW>
    end
end
end
