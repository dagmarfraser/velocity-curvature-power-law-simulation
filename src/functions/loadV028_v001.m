function [results, meta] = loadV028_v001(dataset, opts)
%LOADV028_V001 Load a v028 per-trial result file in the v015 field layout (2026-11-16).
%   [results, meta] = loadV028_v001("Cook CTRL")                       primary arm, window basis
%   loadV028_v001("Fraser", Arm="budget1")                              dataset "Fraser" is read as FraserN48
%   loadV028_v001("Cook CTRL", ObsBasis="whole")                        whole-trial basis (the registered one)
%   loadV028_v001("Cook CTRL", Seed=271828)                             seed rerun (r50)
%   Arms: primary (endMatched, detectability rule), budget1 (endMatched, 1% budget), shapeXu (shapedXuSeg, detectability).
%   The v028 files hold every v015 field (checked 2026-11-16: 33 of 33 present, same trial order as v015) plus window and whole-trial
%   versions. ObsBasis "window" returns the file as saved (betaObs, betaGenStar, betaGenStarMed, ciLo, ciHi, invertStatus, vgfObs are the
%   window basis). ObsBasis "whole" copies betaObsWhole, betaGenStarWhole, betaGenStarMedWhole, ciLoWhole, ciHiWhole, invertStatusWhole,
%   vgfObsWhole over those names, so a script written for v015 runs unchanged on the whole-trial basis. betaObs on the whole trial is
%   bit-identical to v015 betaObs (checked on Cook CTRL, 94 of 94). betaRecCurveMed, betaRecSlice and the forward map are always those
%   of the shortened synthetic. meta carries config, betaGenVec, pipelineLabels, loopCCC, file and the basis chosen.
%   Looks first in src/ beside this file's parent, then in opts.Dir. Fails loud if the file is missing.
arguments
    dataset (1,1) string
    opts.Arm (1,1) string {mustBeMember(opts.Arm, ["primary", "budget1", "shapeXu", "harmonicsFc"])} = "primary"
    opts.ObsBasis (1,1) string {mustBeMember(opts.ObsBasis, ["window", "whole"])} = "window"
    opts.Seed (1,1) double = 1729
    opts.Dir (1,1) string = ""
end
srcDir = fileparts(fileparts(mfilename('fullpath')));
tag = strrep(dataset, " ", "_");  if tag == "Fraser", tag = "FraserN48"; end
switch opts.Arm
    case "primary",     arm = "endMatched_segdetect";
    case "budget1",     arm = "endMatched_bud1";
    case "shapeXu",     arm = "shapedXuSeg_segdetect";
    case "harmonicsFc", arm = "endMatched_segdetect_nHfc";
end
extra = "";  if opts.Seed ~= 1729, extra = "_seed" + opts.Seed + "_r50"; end
name = "loopClosureResults_" + tag + "_all_" + arm + extra + "_v028.mat";
dirs = [srcDir, opts.Dir];  f = "";
for d = dirs(strlength(dirs) > 0), if isfile(fullfile(d, name)), f = fullfile(d, name);  break, end, end
if f == "", error('loadV028:missing', 'v028 file %s not found in %s', name, strjoin(dirs, ' or ')); end
S = load(f);  results = S.results;
if opts.ObsBasis == "whole"
    map = ["betaObs" "betaObsWhole"; "betaGenStar" "betaGenStarWhole"; "betaGenStarMed" "betaGenStarMedWhole"; "ciLo" "ciLoWhole"; ...
           "ciHi" "ciHiWhole"; "invertStatus" "invertStatusWhole"; "vgfObs" "vgfObsWhole"];
    for k = 1:size(map, 1)
        assert(isfield(results, map(k, 2)), 'loadV028:field', '%s lacks %s', name, map(k, 2));
        [results.(map(k, 1))] = results.(map(k, 2));
    end
end
meta = struct('file', f, 'config', S.config, 'betaGenVec', S.betaGenVec, 'pipelineLabels', {S.pipelineLabels}, 'loopCCC', S.loopCCC, ...
    'arm', opts.Arm, 'obsBasis', opts.ObsBasis, 'dataset', tag);
end
