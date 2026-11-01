function summary = runTempoSimAll_v002(opts)
%RUNTEMPOSIMALL_V002  = v001 calling simulateTempoCorpus_v002; default arms AXpF and AXbF (spec 10.44 and 10.45, R23) on Fraser, Hickman PLAC and HALO.
%RUNTEMPOSIMALL_V001  T2/T2b driver: simulate each arm, then run it through the runner.
%   For an interactive MATLAB session (HPC: 72-worker pool). Restartable: the runner
%   skips any arm x dataset whose output already exists (pass ForceRerun to redo).
%   Spec: docs/SPEC_T2_TempoSlopeRecovery_v001.md. Outputs: results/tempoSim/.
%
%   SMOKE TEST (few cores, minutes):
%     runTempoSimAll_v001(Arms="A", Datasets="Zarandi", N_REPS=3, NWorkers=4, ...
%                         SimMaxTrials=8, ForceRerun=true)
%   A smoke corpus records maxTrials in its file; a real run (SimMaxTrials=Inf) refuses
%   to reuse it, so a smoke file can never masquerade as a real arm. Delete
%   results/tempoSim/ after a smoke test.
arguments
    opts.Arms (1,:) string = ["AXpF","AXbF"]   % R23 (spec 10.44, 10.45): decoupled arm, then the seam-matched coupled control; the v001 tier lists are not repeated here
    % cheapest first, Fraser last (restart loses least)
    opts.Datasets (1,:) string = ["Hickman PLAC","Hickman HALO","Fraser"]
    opts.N_REPS (1,1) double = 50
    opts.NWorkers (1,1) double = 72
    opts.SimMaxTrials (1,1) double = Inf
    opts.FraserSessions (1,1) double = 10
    opts.ForceRerun (1,1) logical = false
    opts.SummaryFile (1,1) string = "runTempoSimAll_R23_summary.mat"   % pass a different name for a side run, so the study summary is not overwritten (it was, 17 Oct 2026, by a smoke test)
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
    % ---- simulate whatever is missing; refuse smoke/real mix-ups -------------
    need = false;
    for ds = opts.Datasets
        f = fullfile(simDir, sprintf('tempoSim_%s_%s.mat', arm, strrep(ds,' ','_')));
        if ~isfile(f), need = true; continue; end
        mt = load(f,'maxTrials').maxTrials;
        if ~isequaln(mt, opts.SimMaxTrials)
            error('runTempoSimAll:corpusMismatch', ...
                '%s was simulated with MaxTrials=%g but this run asks for %g. Delete %s.', f, mt, opts.SimMaxTrials, simDir);
        end
    end
    if need
        fprintf('\n--- simulating arm %s ---\n', arm);
        simulateTempoCorpus_v002(arm, Datasets=opts.Datasets, FraserSessions=opts.FraserSessions, MaxTrials=opts.SimMaxTrials);
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
% as runLoopClosureAll_v015: IdleTimeout disabled, fail loud if no pool
p = gcp('nocreate');
if ~isempty(p)
    try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
    if ~alive, try, delete(p); catch, end, p = []; end
end
if isempty(p), p = parpool('local', nWorkers, 'IdleTimeout', Inf); end
if isempty(p) || ~p.Connected, error('runTempoSimAll:noPool','Could not establish a parallel pool.'); end
end
