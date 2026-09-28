function results = exploreTier3PipelineChoice_v001(opts)
% EXPLORETIER3PIPELINECHOICE_V001  Full six-pipeline (+median) spread for
% tier3DesignPrecision_v001's Mode="single" precision report, across all
% four Hickman/Cook arms.
%
% MOTIVATION (2026-09-20): SG-LMLS was checkTier3AgainstMainProject_v001's
% default only because it matches compareHickmanPlacHalo_v002.m's/
% compareCookCtrlAsd_v002.m's own PrimaryPP=4 "sensitivity" pipeline --
% chosen historically for the main project's own reasons, not shown to be
% the right default for a general-purpose tool a researcher would run on
% arbitrary new data. Dagmar's own question: is SG-LMLS defensible here,
% or does picking any one "primary" pipeline as a tool default paper over
% real between-pipeline spread a user should see. This script answers
% that empirically, not by assumption -- for all four Hickman/Cook arms,
% run Mode="single" once per pipeline (all six named, plus "median") and
% report seStage2, semAdequate and the variance decomposition side by
% side, unlike checkTier3AgainstMainProject_v001.m which only ever tests
% one pipeline at a time for cross-validation purposes.
%
% This is read-only against already-frozen data (same six .mat files
% checkTier3AgainstMainProject_v001.m already reads) -- no new analysis,
% a different lens on numbers that already exist.
%
% Fraser, D.S. (2026)

arguments
    opts.SrcDir (1,1) string = "."
end

ppNames = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
choices = [ppNames, "median"];
datasets = ["Hickman_PLAC","Hickman_HALO","Cook_CTRL","Cook_ASD"];

results = struct();
for d = 1:numel(datasets)
    tag = datasets(d);
    batch = loadAsFullBatch_local(opts.SrcDir, tag, "v015", ppNames);
    fprintf('\n=== %s ===\n', tag);
    fprintf('%-10s %10s %10s %8s %8s %8s %8s\n', 'pipeline', 'pointEst', 'seStage2', 'adeq?', 'inv%', 'within%', 'betw%');
    for p = 1:numel(choices)
        cl = tier3DesignPrecision_v001(batch, Mode="single", Pipeline=choices(p));
        results.(matlab.lang.makeValidName(tag)).(matlab.lang.makeValidName(choices(p))) = cl;
        fprintf('%-10s %10.4f %10.4f %8d %8.0f %8.0f %8.0f\n', ...
            choices(p), cl.pointEst, cl.seStage2, cl.semAdequate, ...
            100*cl.fracInversion, 100*cl.fracWithinSubject, 100*cl.fracBetweenSubject);
    end
end
end

function batch = loadAsFullBatch_local(srcDir, tag, version, ppNames)
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag, version));
    if ~isfile(f)
        error('exploreTier3PipelineChoice_v001:NotFound', '%s', sprintf('FAILED PATH: %s', f));
    end
    D = load(f, 'results'); r = D.results(:);
    N = numel(r);
    batch = repmat(struct('ok', true, 'subLabel', "", 'sesLabel', "", 'runLabel', "", 'trial', []), N, 1);
    for i = 1:N
        batch(i).subLabel = string(r(i).subjectID);
        batch(i).sesLabel = tag;
        batch(i).runLabel = string(i);
        bg = getRowVec_local(r(i), 'betaGenStar', ppNames);
        ciLo = getRowVec_local(r(i), 'ciLo', ppNames);
        ciHi = getRowVec_local(r(i), 'ciHi', ppNames);
        batch(i).trial = struct('pipelineLabels', ppNames, 'betaGenStar', bg, 'ciLo', ciLo, 'ciHi', ciHi);
    end
end

function row = getRowVec_local(s, fld, ppNames)
    n = numel(ppNames);
    row = nan(1, n);
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && numel(x) == n
            row = double(x(:)');
        end
    end
end
