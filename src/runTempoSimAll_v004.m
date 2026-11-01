function summary = runTempoSimAll_v004(opts)
%RUNTEMPOSIMALL_V004  = v003 calling simulateTempoCorpus_v004; default arms AXcF, AXaF, AXsF, AXnF (spec 10.55, R30) on Hickman PLAC and Fraser (HALO unread, not simulated).
%   Simulate each arm, then run it through the synthetic runner. For an interactive MATLAB session (HPC, 72-worker pool).
%   Restartable: an arm x dataset whose output exists is skipped (ForceRerun to redo). Outputs: results/tempoSim/.
%   A smoke corpus records maxTrials; a real run (SimMaxTrials=Inf) refuses to reuse it.   Fraser, D.S. (2026)  v004
arguments
    opts.Arms (1,:) string = ["AXcF","AXaF","AXsF","AXnF"]
    opts.Datasets (1,:) string = ["Hickman PLAC","Fraser"]   % cheapest first, Fraser last
    opts.N_REPS (1,1) double = 50
    opts.NWorkers (1,1) double = 72
    opts.SimMaxTrials (1,1) double = Inf
    opts.FraserSessions (1,1) double = 10
    opts.ForceRerun (1,1) logical = false
    opts.SummaryFile (1,1) string = "runTempoSimAll_R30_summary.mat"
end
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir);
simDir = fullfile(srcDir,'..','results','tempoSim');
if ~isfile(fullfile(srcDir,'runLoopClosureFftnoise_v015.m')), makeRunLoopClosureV015_v001('SrcDir',srcDir); end
if ~isfile(fullfile(srcDir,'runLoopClosureSynthetic_v001.m')), makeRunLoopClosureSynthetic_v001('SrcDir',srcDir); end
ensurePool_local(opts.NWorkers);
summary = table('Size',[0 5],'VariableTypes',["string","string","double","double","string"], ...
    'VariableNames',["arm","dataset","nTrials","minutes","note"]);
tAll = tic;
for arm = opts.Arms
    need = false;
    for ds = opts.Datasets
        f = fullfile(simDir, sprintf('tempoSim_%s_%s.mat', arm, strrep(ds,' ','_')));
        if ~isfile(f), need = true; continue; end
        mt = load(f,'maxTrials').maxTrials;
        if ~isequaln(mt, opts.SimMaxTrials)
            error('runTempoSimAll:corpusMismatch','%s was simulated with MaxTrials=%g but this run asks for %g. Delete %s.', f, mt, opts.SimMaxTrials, simDir);
        end
    end
    if need
        fprintf('\n--- simulating arm %s ---\n', arm);
        simulateTempoCorpus_v004(arm, Datasets=opts.Datasets, FraserSessions=opts.FraserSessions, MaxTrials=opts.SimMaxTrials);
    end
    for ds = opts.Datasets
        fprintf('\n=== arm %s | %s ===\n', arm, ds);
        ensurePool_local(opts.NWorkers);
        t0 = tic;
        r = runLoopClosureSynthetic_v001(char("SYN-"+arm+"-"+ds), 'TrialSelection','all', ...
            'N_REPS',opts.N_REPS, 'UseParfor',true, 'ForceRerun',opts.ForceRerun);
        summary(end+1,:) = {arm, ds, numel(r), toc(t0)/60, ""}; %#ok<AGROW>
        save(fullfile(simDir,char(opts.SummaryFile)),'summary');
    end
end
fprintf('\n=== done in %.2f h ===\n', toc(tAll)/3600); disp(summary);
end

function p = ensurePool_local(nWorkers)
p = gcp('nocreate');
if ~isempty(p)
    try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
    if ~alive, try, delete(p); catch, end, p = []; end
end
if isempty(p), p = parpool('local', nWorkers, 'IdleTimeout', Inf); end
if isempty(p) || ~p.Connected, error('runTempoSimAll:noPool','Could not establish a parallel pool.'); end
end
