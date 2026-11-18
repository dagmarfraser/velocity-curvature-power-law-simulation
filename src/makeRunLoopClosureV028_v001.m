function outPath = makeRunLoopClosureV028_v001(opts)
%MAKERUNLOOPCLOSUREV028_V001  Derive runLoopClosureFftnoise_v028 from v027: the segment rule as an option (detectability).
%
%   WHY (Dagmar, 2026-11-15 13:41, "spin up the big iron"): the rerun uses the end-anchored template with the parameter-free
%   detectability rule as primary and the fixed 1% budget as a sensitivity arm (checkSegmentRuleFailures_v002: trials without
%   a segment, detect / 1%: Hickman PLAC 5% / 13%, HALO 4% / 10%, Cook ASD 0% / 4%, others equal; noiseBudgetSensitivity_v001:
%   +1% noise energy moves the median beta_gen* by < 0.001).
%   WHAT v028 CHANGES IN v027, AND ONLY THIS
%     1. Option 'SegmentRule': "auto" (default: v027's behaviour, Tol 1 unless SegmentBudget is given), "tol", "budget" or
%        "detect" (endToEndSegment_v006). The search calls endToEndSegment_v006, which with "auto" equals v005 exactly.
%        Non-default rules add "_seg<rule>" to the output name (the budget keeps its "_bud" tag); stored in config
%        (SEGMENT_RULE) and per trial (emRule, as before).
%     2. Names carry v028. At its defaults v028 is v027.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v027.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v028.m");
if ~isfile(inPath), error("makeV028:missingSource", "v027 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV028:exists", "v028 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"option", "addParameter(p, 'NoiseScale',", ...
    "addParameter(p, 'SegmentRule',    'auto', @(x) ismember(string(x), [""auto"",""tol"",""budget"",""detect""]));   % v028" + nl + "addParameter(p, 'NoiseScale',", 1};
E{end+1} = {"tag", "NOISE_SCALE = opt.NoiseScale;", ...
    "SEG_RULE = string(opt.SegmentRule);  if ~ismember(SEG_RULE, [""auto"", ""budget""]), OPT_TAG = OPT_TAG + ""_seg"" + SEG_RULE; end   % v028" + nl + "NOISE_SCALE = opt.NoiseScale;", 1};
E{end+1} = {"calls", "RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE);", "RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE, SEG_RULE);", 4};
E{end+1} = {"signature", "rngSeed, nHopt, segBudget, noiseScale)", "rngSeed, nHopt, segBudget, noiseScale, segRule)", 1};
E{end+1} = {"segment call", "endToEndSegment_v005([resid_maj(1:Mc), resid_min(1:Mc)], Tol=1, Budget=segBudget);", ...
    "endToEndSegment_v006([resid_maj(1:Mc), resid_min(1:Mc)], Rule=segRule, Tol=1, Budget=segBudget);   % v028: v006 (auto = v005)", 1};
E{end+1} = {"config", "'NOISE_SCALE', NOISE_SCALE);", "'NOISE_SCALE', NOISE_SCALE, 'SEGMENT_RULE', SEG_RULE);   % v028", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV028:anchor", "%s: expected %d match(es), found %d. v027 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v027");  nMat = count(s, "_v027.mat");
if nName < 5 || nMat < 2, error("makeV028:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v027", "runLoopClosureFftnoise_v028");  s = replace(s, "_v027.mat", "_v028.mat");
s = replace(s, "fprintf('[v027] Output exists", "fprintf('[v028] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v028(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v028(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v028  v027 + SegmentRule option (auto | tol | budget | detect; endToEndSegment_v006). GENERATED: see makeRunLoopClosureV028_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV028:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
