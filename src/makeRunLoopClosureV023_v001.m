function outPath = makeRunLoopClosureV023_v001(opts)
%MAKERUNLOOPCLOSUREV023_V001  Derive runLoopClosureFftnoise_v023 from v022: the modernised shapeXu sensitivity arm (2026-11-15).
%
%   WHY: Dagmar (2026-11-15 09:10) asked for the shapeXu sensitivity arm to be modernised so that it, too, is examined from
%   beta_obs on the same-length window (primary) and on the whole trial (registered basis). The existing v015/v020 shaped_xu
%   runs use full-length synthetics and whole-trial beta_obs only, so they cannot be compared like with like.
%   WHAT v023 CHANGES IN v022, AND ONLY THIS
%     1. NoiseModel 'shapedXuSeg' (new). Everything is as 'endMatched' (same per-axis segments and common length L from
%        endToEndSegment_v004, synthetic shortened to L, beta_obs on the central window and on the whole trial, two
%        inversions against the same forward map, the same per-cell substreams) EXCEPT the noise: each replicate is the
%        production shaped_xu construction (generateLoopClosureNoise_v003('shaped_xu'): the segment's |FFT| on Xu fractional
%        phase) applied to each axis's segment, with the Xu carrier at the segment's alpha (alphaPWLW_v001, as endMatched
%        selects on). No alpha selection (production shaped_xu had none) and no R7 gate. So the two arms differ in the noise
%        construction only.
%     2. The checks that switched endMatched behaviour on (segment block, smokeA, whole-trial inversion) now accept either
%        model. The default stays 'shaped_xu', and 'endMatched' is unchanged, so v023 reproduces v022 for both.
%     3. Names: function, error ids, output files and messages carry v023.
%   ACCEPTANCE: checkV023Acceptance_v001: defaults and endMatched equal v022 exactly (same machine, serial); shapedXuSeg
%   runs with all fields and differs from endMatched. Fail Loud: each replacement asserts its match count.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v022.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v023.m");
if ~isfile(inPath), error("makeV023:missingSource", "v022 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV023:exists", "v023 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"NoiseModel option", "@(x) ismember(x, {'auto','fftnoise','xu','shaped_xu','bootstrap','endMatched'}));   % v022: endMatched", ...
    "@(x) ismember(x, {'auto','fftnoise','xu','shaped_xu','bootstrap','endMatched','shapedXuSeg'}));   % v022: endMatched; v023: shapedXuSeg", 1};
E{end+1} = {"segment block", "    if strcmp(noiseModel, 'endMatched')   % v022", "    if ismember(noiseModel, {'endMatched', 'shapedXuSeg'})   % v022; v023: shapedXuSeg too", 1};
E{end+1} = {"smokeA", "    if strcmp(noiseModel, 'endMatched'), smokeA(:) = NaN; end", "    if ismember(noiseModel, {'endMatched', 'shapedXuSeg'}), smokeA(:) = NaN; end   % v023: shapedXuSeg too;", 1};
E{end+1} = {"whole-trial inversion", "    if ~strcmp(noiseModel, 'endMatched'), return, end", "    if ~ismember(noiseModel, {'endMatched', 'shapedXuSeg'}), return, end   % v023: shapedXuSeg too", 1};
caseOld = "                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end" + nl + "                case 'bootstrap'";
caseNew = "                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end" + nl + ...
    "                case 'shapedXuSeg'   % v023: production shaped_xu on each axis's segment, carrier at the segment's PWLW alpha, same substreams" + nl + ...
    "                    rsC = RandStream(""threefry4x64_20"", ""Seed"", rngSeed);  rsC.Substream = tr.cellKey * 1e5 + bi * 1e3 + rep;" + nl + ...
    "                    rsOld = RandStream.setGlobalStream(rsC);" + nl + ...
    "                    [nMaj, nMin] = generateLoopClosureNoise_v003('shaped_xu', segMaj, segMin, FS, aSeg(1), aSeg(2));" + nl + ...
    "                    RandStream.setGlobalStream(rsOld);" + nl + ...
    "                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end" + nl + ...
    "                case 'bootstrap'";
E{end+1} = {"shapedXuSeg noise case", caseOld, caseNew, 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV023:anchor", "%s: expected %d match(es), found %d. v022 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v022");  nMat = count(s, "_v022.mat");
if nName < 5 || nMat < 2, error("makeV023:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v022", "runLoopClosureFftnoise_v023");  s = replace(s, "_v022.mat", "_v023.mat");
s = replace(s, "fprintf('[v022] Output exists", "fprintf('[v023] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v023(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v023(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v023  v022 + NoiseModel 'shapedXuSeg' (the shapeXu sensitivity arm on the same segments, length, windows and seeds as 'endMatched'). GENERATED: see makeRunLoopClosureV023_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV023:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
