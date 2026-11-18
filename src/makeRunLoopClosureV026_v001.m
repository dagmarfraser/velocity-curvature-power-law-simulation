function outPath = makeRunLoopClosureV026_v001(opts)
%MAKERUNLOOPCLOSUREV026_V001  Derive runLoopClosureFftnoise_v026 from v025: fix v024's window bug; segment budget option.
%
%   BUG FIXED (found 2026-11-15 on BlueBEAR, jobs 61831718-61831725): makeRunLoopClosureV024_v001 replaced the start of the
%   line "c0 = ...;  wA = ...;  wB = ...;" with a version ending in a comment, which commented out the definitions of wA and
%   wB. Every endMatched/shapedXuSeg trial that reached that line stopped with "Unrecognized function or variable 'wA'", so
%   v024 and v025 produced no end-matched inversions. checkV024Acceptance_v001 passed because it skipped failed trials
%   (a vacuous pass); checkV026Acceptance_v001 requires real trials to run.
%   WHAT v026 CHANGES IN v025, AND ONLY THIS
%     1. The window line restored: c0 centred on the covered residual (v024's intent), wA = c0 + tsClip - 1, wB = wA + L - 1.
%     2. Option 'SegmentBudget' (default []: the fixed Tol 1 rule of v022-v025). A positive number switches the segment search
%        to endToEndSegment_v005's budget rule (the wrap adds at most that share of each axis's first- and second-difference
%        energy; 0.01 = 1%). Non-default values add "_bud<value*100>" to the output name; stored in config (SEGMENT_BUDGET)
%        and per trial (emRule). The search uses endToEndSegment_v005 in both cases (v005 without a budget = v004 exactly,
%        checked on 7,986 trial-template rows by checkSegmentRuleFailures_v001).
%     3. Names carry v026.
%   ACCEPTANCE: checkV026Acceptance_v001.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v025.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v026.m");
if ~isfile(inPath), error("makeV026:missingSource", "v025 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV026:exists", "v026 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"window line (bug)", "            c0 = floor((Mc - Lsyn) / 2) + 1;   % v024: centred on the covered span  wA = c0 + tsClip - 1;  wB = wA + Lsyn - 1;   % central window of the real trial, length Lsyn", ...
    "            c0 = floor((Mc - Lsyn) / 2) + 1;  wA = c0 + tsClip - 1;  wB = wA + Lsyn - 1;   % v026: restored (v024 had commented out wA, wB); central window of the covered span, length Lsyn", 1};
E{end+1} = {"option", "addParameter(p, 'NHarmonics',", ...
    "addParameter(p, 'SegmentBudget',  [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0 && x < 1));   % v026" + nl + "addParameter(p, 'NHarmonics',", 1};
E{end+1} = {"tag", "NH_OPT = opt.NHarmonics;", ...
    "SEG_BUDGET = opt.SegmentBudget;  if ~isempty(SEG_BUDGET), OPT_TAG = OPT_TAG + ""_bud"" + string(round(100 * SEG_BUDGET, 4)); end   % v026" + nl + "NH_OPT = opt.NHarmonics;", 1};
E{end+1} = {"calls", "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT);", ...
    "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT, SEG_BUDGET);", 4};
E{end+1} = {"signature", "noiseModel, endAnchored, carrierAlpha, rngSeed, nHopt)", "noiseModel, endAnchored, carrierAlpha, rngSeed, nHopt, segBudget)", 1};
E{end+1} = {"segment call", "[sI0, sI1, emInfo] = endToEndSegment_v004([resid_maj(1:Mc), resid_min(1:Mc)], Tol=1, Mode=""commonLength"");", ...
    "[sI0, sI1, emInfo] = endToEndSegment_v005([resid_maj(1:Mc), resid_min(1:Mc)], Tol=1, Budget=segBudget);   % v026: v005 (= v004 when segBudget is empty)" + nl + ...
    "            res.emRule = emInfo.rule;", 1};
E{end+1} = {"emRule field", "'emWindow', [NaN NaN]; 'emTail', NaN;", "'emWindow', [NaN NaN]; 'emTail', NaN; 'emRule', """";", 1};
E{end+1} = {"config", "'NHARMONICS', NH_OPT);", "'NHARMONICS', NH_OPT, 'SEGMENT_BUDGET', SEG_BUDGET);   % v026", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV026:anchor", "%s: expected %d match(es), found %d. v025 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v025");  nMat = count(s, "_v025.mat");
if nName < 5 || nMat < 2, error("makeV026:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v025", "runLoopClosureFftnoise_v026");  s = replace(s, "_v025.mat", "_v026.mat");
s = replace(s, "fprintf('[v025] Output exists", "fprintf('[v026] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v026(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v026(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v026  v025 + v024's window bug fixed + SegmentBudget option (endToEndSegment_v005). GENERATED: see makeRunLoopClosureV026_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV026:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
