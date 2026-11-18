function outPath = makeRunLoopClosureV033_v001(opts)
%MAKERUNLOOPCLOSUREV033_V001  Derive runLoopClosureFftnoise_v033 from v028: the template harmonic count may be 0 (2026-11-18).
%
%   WHY (Dagmar, 2026-11-18, "we probably need a 0 harmonic run"): the harmonic-count ladder (1, 2, 3, 4, 6, 10, fc) has no baseline. K = 0 strips no
%   harmonic: the windowed template is only a constant and a linear trend per Hann window, so the residual is the trend-removed signal. It shows
%   how the end-matched segment search, the segment exponent and everything downstream behave when the movement's own harmonics are left in.
%   WHAT v033 CHANGES IN v028, AND ONLY THIS
%     1. Option 'NHarmonics' accepts 0 (was: integers from 1, or 'fc'). The validator's lower bound goes from 1 to 0. The template code
%        (templateSubtract_local) already builds a constant and a trend and loops h = 1:nH, so nH = 0 needs no other change.
%     2. Names: the function, its output files and its messages say v033 ("_nH0_v033.mat" for K = 0; other counts write "_nH<K>_v033.mat").
%   ACCEPTANCE: at NHarmonics = 4 (the default) v033 equals v028 (checked by checkV033EqualsV028_v001 before any K = 0 job is submitted).
%   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v028.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v033.m");
if ~isfile(inPath), error("makeV033:missingSource", "v028 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV033:exists", "v033 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
from = "addParameter(p, 'NHarmonics',     4, @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x == round(x)) || strcmp(string(x), ""fc""));   % v025";
to   = "addParameter(p, 'NHarmonics',     4, @(x) (isnumeric(x) && isscalar(x) && x >= 0 && x == round(x)) || strcmp(string(x), ""fc""));   % v025; v033: 0 allowed (constant and trend only)";
n = count(s, from);
if n ~= 1, error("makeV033:anchor", "NHarmonics validator: expected 1 match, found %d. v028 has drifted; nothing written.", n); end
s = replace(s, from, to);
nName = count(s, "runLoopClosureFftnoise_v028");  nMat = count(s, "_v028.mat");
if nName < 5 || nMat < 2, error("makeV033:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v028", "runLoopClosureFftnoise_v033");  s = replace(s, "_v028.mat", "_v033.mat");
s = replace(s, "fprintf('[v028] Output exists", "fprintf('[v033] Output exists");
hdr = "function results = runLoopClosureFftnoise_v033(datasetName, varargin)" + nl;
if count(s, hdr) ~= 1, error("makeV033:header", "function header not found exactly once; nothing written."); end
s = replace(s, hdr, hdr + "% runLoopClosureFftnoise_v033  v028 + NHarmonics may be 0 (baseline: constant and trend only). GENERATED: see makeRunLoopClosureV033_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV033:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
