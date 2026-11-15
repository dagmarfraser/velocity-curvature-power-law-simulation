%% checkSimexJoinSpec_v001  Deep Review GD17 item 7a: what the SIMEX-to-monotonicity join counts mean.
% Reads results/crossReferenceSimexMonotonicity_v002.mat (V2 PRIMARY = places each trial at its own constellation-median beta_gen*; a trial with none is
% "unplaced" for all six pipelines and left out; per-pipeline invRise then classifies the placed cells, so non-invertible cells remain). For each corpus
% it records the SIMEX cells, the cells dropped as unplaced, the cells joined, and the number of trials behind the dropped cells (V1 rows minus V2 rows).
% Fails loud if dropped cells are not 6 per unplaced trial or the counts do not add up.
% Run from src/. Writes results/checkSimexJoinSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
f = fullfile(ROOT, "results", "crossReferenceSimexMonotonicity_v002.mat"); if ~isfile(f), error("simexSpec:file", "%s", "FAILED PATH: " + f); end
S = load(f);  res = S.res;  AJ = S.allJoined;  NPIPE = 6;
rows = res(res.variant == "V2 PRIMARY", :);  out = table();
for c = 1:height(rows)
    J1 = AJ{c, 2};  J2 = AJ{c, 3};
    k1 = J1.trialID + "|" + J1.pipeline;  k2 = J2.trialID + "|" + J2.pipeline;
    dropped = setdiff(k1, k2);  tr = unique(extractBefore(dropped, "|"));
    assert(numel(dropped) == rows.nUnplaced(c), "simexSpec:cells", "%s: dropped cells differ from nUnplaced", rows.group(c));
    assert(rows.nSimex(c) == height(J2) + rows.nUnplaced(c), "simexSpec:sum", "%s: nSimex ~= joined + unplaced", rows.group(c));
    assert(isempty(dropped) || numel(dropped) == NPIPE*numel(tr), "simexSpec:six", "%s: dropped cells are not 6 per trial", rows.group(c));
    out = [out; table(rows.group(c), rows.nSimex(c), height(J2), rows.nUnplaced(c), numel(tr), rows.nInv(c), rows.nNon(c), ...
        'VariableNames', ["group" "nSimex" "nJoined" "nUnplaced" "nUnplacedTrials" "nInv" "nNon"])]; %#ok<AGROW>
end
disp(out);  save(fullfile(ROOT, "results", "checkSimexJoinSpec_v001.mat"), "out");
