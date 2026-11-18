function outPath = makeRunLoopClosureV027_v001(opts)
%MAKERUNLOOPCLOSUREV027_V001  Derive runLoopClosureFftnoise_v027 from v026: a NoiseScale option for the end-matched noise.
%
%   WHY (Dagmar, 2026-11-15 10:44: "can we ground the 1% vs 2% in something more principled?"): the segment budget is a
%   bound on how much the wrap can add to the derivative noise energy. Its consequence is measured by inflating the
%   end-matched noise's variance by that share and reading the beta_gen* shift (noiseBudgetSensitivity_v001).
%   WHAT v027 CHANGES IN v026, AND ONLY THIS
%     1. Option 'NoiseScale' (default 1): each 'endMatched' noise draw is multiplied by sqrt(NoiseScale), so its variance
%        (and its velocity and acceleration energy) rises by NoiseScale - 1. The draws themselves are unchanged (same
%        substreams), so runs at different scales differ only by the scaling. Non-default values add "_ns<value>" to the
%        output name; stored in config (NOISE_SCALE).
%     2. Names carry v027. At the default v027 is v026.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v026.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v027.m");
if ~isfile(inPath), error("makeV027:missingSource", "v026 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV027:exists", "v027 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"option", "addParameter(p, 'SegmentBudget',", ...
    "addParameter(p, 'NoiseScale',     1, @(x) isnumeric(x) && isscalar(x) && x > 0);   % v027" + nl + "addParameter(p, 'SegmentBudget',", 1};
E{end+1} = {"tag", "NH_OPT = opt.NHarmonics;", ...
    "NOISE_SCALE = opt.NoiseScale;  if NOISE_SCALE ~= 1, OPT_TAG = OPT_TAG + ""_ns"" + string(NOISE_SCALE); end   % v027" + nl + "NH_OPT = opt.NHarmonics;", 1};
E{end+1} = {"calls", "RNG_SEED, NH_OPT, SEG_BUDGET);", "RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE);", 4};
E{end+1} = {"signature", "rngSeed, nHopt, segBudget)", "rngSeed, nHopt, segBudget, noiseScale)", 1};
E{end+1} = {"scale", "                    [nMaj, nMin] = endMatchedNoise_v001(segMaj, segMin, FS, aSeg(1), aSeg(2));", ...
    "                    [nMaj, nMin] = endMatchedNoise_v001(segMaj, segMin, FS, aSeg(1), aSeg(2));" + nl + ...
    "                    nMaj = nMaj * sqrt(noiseScale);  nMin = nMin * sqrt(noiseScale);   % v027: variance x noiseScale (1 = unchanged)", 1};
E{end+1} = {"config", "'SEGMENT_BUDGET', SEG_BUDGET);", "'SEGMENT_BUDGET', SEG_BUDGET, 'NOISE_SCALE', NOISE_SCALE);   % v027", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV027:anchor", "%s: expected %d match(es), found %d. v026 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v026");  nMat = count(s, "_v026.mat");
if nName < 5 || nMat < 2, error("makeV027:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v026", "runLoopClosureFftnoise_v027");  s = replace(s, "_v026.mat", "_v027.mat");
s = replace(s, "fprintf('[v026] Output exists", "fprintf('[v027] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v027(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v027(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v027  v026 + NoiseScale option (end-matched noise variance multiplier, for the budget sensitivity). GENERATED: see makeRunLoopClosureV027_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV027:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
