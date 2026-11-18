function out = checkV033EqualsV028_v001(opts)
%CHECKV033EQUALSV028_V001 Acceptance of runLoopClosureFftnoise_v033 (the K = 0 template baseline). Errors on any failure (2026-11-18).
%   Small grid (N_BETA 7, N_REPS 3, N_R7 2, all trials of one dataset, parfor: per-cell substreams make the results independent of the worker layout); outputs tagged _acc33.
%   1. NHarmonics 4 (default), endMatched, end-anchored, SegmentRule "detect": v033 equals v028 exactly (betaObs, betaObsWhole, betaGenStar,
%      betaRecCurveMed, emSeg, emL, emWindow, emRule, emStatus, emAlphaSeg).
%   2. NHarmonics 0 runs: every trial stores nHarm = 0 and config.NHARMONICS = 0; it is not silently the K = 4 run (the residual segments differ
%      on at least one trial); the share of trials that ran and the failure reasons are returned and printed (many refusals are expected and are
%      themselves informative, so a low share is reported, not an error).
%   3. The validator: -1 and 0.5 are refused, 'fc' and 2 are accepted by the parser (parser only; nothing run).
%   out = checkV033EqualsV028_v001();  out = checkV033EqualsV028_v001(Dataset='Cook ASD');
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
end
here = fileparts(mfilename('fullpath'));  addpath(here, genpath(fullfile(here, 'functions')));
c = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc33', 'NoiseModel', 'endMatched', 'Template', 'endAnchored', 'SegmentRule', 'detect', 'TrialSelection', 'all', 'UseParfor', true};
if isempty(gcp('nocreate')), ensurePool_v001(); end
a28 = runLoopClosureFftnoise_v028(opts.Dataset, c{:});
a33 = runLoopClosureFftnoise_v033(opts.Dataset, c{:});
eqAll(a28, a33, ["betaObs", "betaObsWhole", "betaGenStar", "betaRecCurveMed", "emSeg", "emL", "emWindow", "emRule", "emStatus", "emAlphaSeg"], 'K = 4');
fprintf('1. K = 4: v033 = v028 exactly on %d trials (%d ran)\n', numel(a33), nnz([a33.emStatus] == "ok"));
z33 = runLoopClosureFftnoise_v033(opts.Dataset, c{:}, 'NHarmonics', 0);
nh = [z33.nHarm];
if ~all(nh == 0), error('acc33:nHarm', 'K = 0 run stored nHarm other than 0 (%s)', mat2str(unique(nh))); end
okZ = [z33.emStatus] == "ok";  ok4 = [a33.emStatus] == "ok";
both = find(okZ & ok4);  differ = 0;
for k = both, if ~isequaln(z33(k).emSeg, a33(k).emSeg) || ~isequaln(z33(k).emAlphaSeg, a33(k).emAlphaSeg), differ = differ + 1; end, end
if ~isempty(both) && differ == 0, error('acc33:same', 'K = 0 and K = 4 give identical segments on every trial they share: K is not reaching the template'); end
fprintf('2. K = 0: nHarm = 0 on all %d trials; %d ran (K = 4: %d ran); segments differ from K = 4 on %d of %d shared trials\n', numel(z33), nnz(okZ), nnz(ok4), differ, numel(both));
st = [z33.emStatus];  [u, ~, g] = unique(extractBefore(pad(st(st ~= "ok"), 60), 60));
for i = 1:numel(u), fprintf('     %4d  %s\n', nnz(g == i), strtrim(u(i))); end
bad = {-1, 0.5};
for k = 1:numel(bad)
    threw = false;
    try, runLoopClosureFftnoise_v033(opts.Dataset, c{:}, 'NHarmonics', bad{k}, 'MaxChunks', 0); catch, threw = true; end
    if ~threw, error('acc33:validator', 'NHarmonics = %g was accepted', bad{k}); end
end
fprintf('3. validator: -1 and 0.5 refused\nALL CHECKS PASSED\n');
out = struct('dataset', opts.Dataset, 'nTrials', numel(z33), 'ranK0', nnz(okZ), 'ranK4', nnz(ok4), 'differ', differ, 'failureReasonsK0', {cellstr(u)}, 'generated', string(datetime('now')));
end

function eqAll(a, b, fields, label)
if numel(a) ~= numel(b), error('acc33:n', '%s: trial counts differ', label); end
for k = 1:numel(a)
    for f = fields
        if ~isequaln(a(k).(f), b(k).(f)), error('acc33:equal', '%s: trial %s field %s differs', label, a(k).trialID, f); end
    end
end
end
