function verdictTable = compare42CellVerdict_v001()
% COMPARE42CELLVERDICT_V001  Does Finding #160's own PASS/CONDITIONAL/FAIL
% table change if computed from N_REPS=200 data instead of N_REPS=20?
%
% Replicates checkInvertibilityForEmpirical_v2_002.m's own methodology
% exactly (rise-coverage only, PASS>=95%, CONDITIONAL>=90%, FAIL<90%,
% no_beta_obs-equivalent trials excluded from the denominator) on BOTH
% the existing v012 (N=20) corpus and the new v013/v014 (N=200) corpus,
% for direct comparison against Finding #160's own citable table.
%
% Fraser, D.S. (2026)

srcDir = fileparts(mfilename('fullpath'));
datasets = {'Zarandi','Cook_CTRL','Cook_ASD','Dhieb','Hickman_PLAC','Hickman_HALO','Fraser'};
newVersion = {'v013','v013','v014','v014','v013','v013','v013'};
pipelineLabels = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];

rows = cell(numel(datasets)*6, 1);
rowI = 0;

for i = 1:numel(datasets)
    ds = datasets{i};
    Sold = load(fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v012.mat', ds)), 'results');
    Snew = load(fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', ds, newVersion{i})), 'results');

    for pp = 1:6
        [verdictOld, covOld, nValidOld] = verdictOneCell_local(Sold.results, pp);
        [verdictNew, covNew, nValidNew] = verdictOneCell_local(Snew.results, pp);

        rowI = rowI + 1;
        rows{rowI} = struct('dataset', string(strrep(ds,'_',' ')), 'pipeline', pipelineLabels(pp), ...
            'verdictN20', verdictOld, 'covN20', covOld, 'nValidN20', nValidOld, ...
            'verdictN200', verdictNew, 'covN200', covNew, 'nValidN200', nValidNew, ...
            'changed', verdictOld ~= verdictNew);
    end
end

verdictTable = struct2table([rows{:}]);
fprintf('=== 42-cell PASS/CONDITIONAL/FAIL: N=20 (v012) vs N=200 (v013/v014) ===\n\n');
disp(verdictTable);

nChanged = sum(verdictTable.changed);
fprintf('\nCells with a DIFFERENT verdict at N=200 vs N=20: %d/42\n', nChanged);
if nChanged > 0
    fprintf('\nChanged cells:\n');
    disp(verdictTable(verdictTable.changed, :));
end

end

function [verdict, coverage, nValid] = verdictOneCell_local(results, pp)
    statuses = arrayfun(@(r) r.invertStatus(pp), results);
    valid = statuses ~= "no_beta_obs";
    nValid = sum(valid);
    if nValid == 0
        verdict = "no-data"; coverage = NaN; return
    end
    nRise = sum(statuses(valid) == "rise");
    coverage = nRise / nValid;
    if coverage >= 0.95
        verdict = "PASS";
    elseif coverage >= 0.90
        verdict = "CONDITIONAL";
    else
        verdict = "FAIL";
    end
end
