function verdictTable = compare42CellVerdict_v002(opts)
% COMPARE42CELLVERDICT_V002  42-cell PASS/CONDITIONAL/FAIL table at N_REPS=20
% (v012, Finding #160's citable table) against an N_REPS=200 corpus.
%
% Fork of compare42CellVerdict_v001.m. Only change: the N_REPS=200 corpus is
% an argument, defaulting to v015 (the paper's N_REPS=200 basis since
% Finding #216), instead of the hardcoded v013/v014 mix. Verdict rule is
% identical (rise-coverage over trials with a valid invertStatus; PASS >=95%,
% CONDITIONAL >=90%, FAIL below).
%
%   T = compare42CellVerdict_v002();                    % v015
%   T = compare42CellVerdict_v002(N200Version="v013v014");  % v001 behaviour
%
% Fraser, D.S. (2026)

arguments
    opts.N200Version (1,1) string {mustBeMember(opts.N200Version, ["v015","v013v014"])} = "v015"
end

srcDir   = fileparts(mfilename('fullpath'));
datasets = {'Zarandi','Cook_CTRL','Cook_ASD','Dhieb','Hickman_PLAC','Hickman_HALO','Fraser'};
if opts.N200Version == "v015"
    newVersion = repmat({'v015'}, 1, numel(datasets));
else
    newVersion = {'v013','v013','v014','v014','v013','v013','v013'};
end
pipelineLabels = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];

rows = cell(numel(datasets)*6, 1);
rowI = 0;
for i = 1:numel(datasets)
    ds    = datasets{i};
    fOld  = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v012.mat', ds));
    fNew  = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', ds, newVersion{i}));
    Sold  = loadResults_local(fOld);
    Snew  = loadResults_local(fNew);
    for pp = 1:6
        [vOld, cOld, nOld] = verdictOneCell_local(Sold, pp, fOld);
        [vNew, cNew, nNew] = verdictOneCell_local(Snew, pp, fNew);
        rowI = rowI + 1;
        rows{rowI} = struct('dataset', string(strrep(ds,'_',' ')), 'pipeline', pipelineLabels(pp), ...
            'verdictN20', vOld, 'covN20', cOld, 'nValidN20', nOld, ...
            'verdictN200', vNew, 'covN200', cNew, 'nValidN200', nNew, ...
            'changed', vOld ~= vNew);
    end
end

verdictTable = struct2table([rows{:}]);
verdictTable.Properties.Description = "N=20: v012; N=200: " + opts.N200Version;

fprintf('=== 42-cell verdicts: N=20 (v012) vs N=200 (%s) ===\n\n', opts.N200Version);
disp(verdictTable);
for col = ["verdictN20","verdictN200"]
    v = verdictTable.(col);
    fprintf('%-12s PASS %2d | CONDITIONAL %2d | FAIL %2d | no-data %d\n', col, ...
        sum(v=="PASS"), sum(v=="CONDITIONAL"), sum(v=="FAIL"), sum(v=="no-data"));
end
fprintf('Cells with a different verdict: %d/42\n', sum(verdictTable.changed));
end

%% ------------------------------------------------------------------------
function results = loadResults_local(f)
    if ~isfile(f)
        error('compare42CellVerdict:MissingCorpus', 'FAILED PATH: %s', f);
    end
    S = load(f, 'results');
    if ~isfield(S, 'results')
        error('compare42CellVerdict:NoResults', 'No ''results'' variable in %s', f);
    end
    results = S.results;
end

function [verdict, coverage, nValid] = verdictOneCell_local(results, pp, f)
    if any(arrayfun(@(r) numel(r.invertStatus), results) ~= 6)
        error('compare42CellVerdict:PipelineCount', ...
            'invertStatus is not 6 pipelines wide for every trial in %s', f);
    end
    statuses = arrayfun(@(r) r.invertStatus(pp), results);
    valid  = statuses ~= "no_beta_obs";
    nValid = sum(valid);
    if nValid == 0
        verdict = "no-data"; coverage = NaN; return
    end
    coverage = sum(statuses(valid) == "rise") / nValid;
    if coverage >= 0.95
        verdict = "PASS";
    elseif coverage >= 0.90
        verdict = "CONDITIONAL";
    else
        verdict = "FAIL";
    end
end
