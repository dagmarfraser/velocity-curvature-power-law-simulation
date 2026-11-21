% v002: v001 of this test used isequal, which is false whenever a field holds NaN; replaced by isequaln (found by diagnosis, 2026-11-21).
% testLoadV028_v002_vs_v001_v002 (2026-11-21): for K=4 arms, loadV028_v002 must return exactly what loadV028_v001 returned (isequaln, NaN==NaN, on results and shared meta fields).
old = fullfile(getenv('HOME'),'Library','CloudStorage','Dropbox','Brain2Bee','PowerLawSimulationPreReg','src','functions');
here = fileparts(mfilename('fullpath'));
ds = ["Fraser" "Zarandi" "Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO" "Dhieb"];
arms = ["primary" "budget1" "shapeXu"];  bases = ["window" "whole"];
rows = {};  nFail = 0;
for a = arms, for b = bases, for d = ds
    addpath(old);  [r1, m1] = loadV028_v001(d, Arm=a, ObsBasis=b);  rmpath(old);
    addpath(here); [r2, m2] = loadV028_v002(d, Arm=a, ObsBasis=b);  rmpath(here);
    fn = fieldnames(m1);  okm = true;  for i = 1:numel(fn), okm = okm && isequaln(m1.(fn{i}), m2.(fn{i})); end
    ok = isequaln(r1, r2) && okm;  nFail = nFail + ~ok;
    rows(end+1,:) = {d, a, b, numel(r1), ok}; %#ok<SAGROW>
    fprintf('%-13s %-8s %-6s n=%d  identical=%d\n', d, a, b, numel(r1), ok);
end, end, end
T = cell2table(rows, 'VariableNames', {'dataset','arm','basis','nTrials','identical'});
writetable(T, fullfile(fileparts(here), 'outputs_K1', 'loadV028_v002_vs_v001_check_K1_20261121_v002.csv'));
% K=1 sanity: loads, md5 checked, trial count equals K=4
for d = ds
    [rk, mk] = loadV028_v002(d, NHarm=1);  [r4] = loadV028_v002(d);
    fprintf('K=1 %-13s n=%d (K4 n=%d) md5Checked=%d NHARMONICS=%g\n', d, numel(rk), numel(r4), mk.md5Checked, mk.config.NHARMONICS);
    assert(mk.md5Checked && numel(rk) == numel(r4) && mk.config.NHARMONICS == 1);
end
% fail-loud checks
try, loadV028_v002("Cook ASD", NHarm=1, Arm="budget1"); error('NOT RAISED'); catch e, fprintf('expected error: %s\n', e.identifier); assert(e.identifier == "loadV028:nharmArm"); end
try, loadV028_v002("Cook ASD", NHarm=1, Dir=fullfile(getenv('HOME'),'nonexistent')); error('NOT RAISED'); catch e, fprintf('expected error: %s\n', e.identifier); assert(e.identifier == "loadV028:missing"); end
fprintf('TOTAL non-identical: %d of %d\n', nFail, height(T));
assert(nFail == 0);
