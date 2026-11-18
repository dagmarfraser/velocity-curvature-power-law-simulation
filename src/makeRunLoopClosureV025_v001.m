function outPath = makeRunLoopClosureV025_v001(opts)
%MAKERUNLOOPCLOSUREV025_V001  Derive runLoopClosureFftnoise_v025 from v024: the template's harmonic count as an option.
%
%   WHY (template comparison, 2026-11-15, templateComparison_v002): with 4 harmonics the residual keeps movement power at k*f0
%   (3-9% of acceleration power by arm); with every harmonic up to the BWFD cutoff it loses 14-19% of acceleration noise in
%   the notches. Stage 2 compares beta_gen* across templates, which needs the runner to take the harmonic count.
%   WHAT v025 CHANGES IN v024, AND ONLY THIS
%     1. Option 'NHarmonics': a positive integer (default 4, v024's literal) or 'fc' = max(1, floor(BW cutoff / f0)) per trial
%        (the BW cutoff is derivConfigs(1).filterParams(2), 10 Hz by default). Non-default values add "_nH<value>" to the output
%        name. The count used is stored per trial (nHarm) and in config (NHARMONICS).
%     2. The template call uses that count instead of the literal 4. Nothing else changes, so v025 at its defaults is v024.
%   ACCEPTANCE: v025 defaults equal v024 exactly (checked by the stage-2 script templateComparisonBeta_v001 before it runs).
%   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v024.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v025.m");
if ~isfile(inPath), error("makeV025:missingSource", "v024 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV025:exists", "v025 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"option", "addParameter(p, 'Template',       'endAnchored', @(x) ismember(x, {'endAnchored','production'}));   % v020", ...
    "addParameter(p, 'Template',       'endAnchored', @(x) ismember(x, {'endAnchored','production'}));   % v020" + nl + ...
    "addParameter(p, 'NHarmonics',     4, @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x == round(x)) || strcmp(string(x), ""fc""));   % v025", 1};
E{end+1} = {"tag", "if ~END_ANCHORED, OPT_TAG = OPT_TAG + ""_tplProd""; end   % v020: v015's production windows", ...
    "if ~END_ANCHORED, OPT_TAG = OPT_TAG + ""_tplProd""; end   % v020: v015's production windows" + nl + ...
    "NH_OPT = opt.NHarmonics;  if ~(isnumeric(NH_OPT) && NH_OPT == 4), OPT_TAG = OPT_TAG + ""_nH"" + string(NH_OPT); end   % v025", 1};
E{end+1} = {"calls", "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED);", ...
    "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT);", 4};
E{end+1} = {"signature", "noiseModel, endAnchored, carrierAlpha, rngSeed)", "noiseModel, endAnchored, carrierAlpha, rngSeed, nHopt)", 1};
E{end+1} = {"count", "    res.f0 = estimateF0_local(double(tr.x), double(tr.y), FS);", ...
    "    res.f0 = estimateF0_local(double(tr.x), double(tr.y), FS);" + nl + ...
    "    if isnumeric(nHopt), nHarm = nHopt; else, nHarm = max(1, floor(derivConfigs(1).filterParams(2) / res.f0)); end   % v025" + nl + ...
    "    res.nHarm = nHarm;", 1};
E{end+1} = {"template call", "        double(tr.x), double(tr.y), FS, res.f0, 4, 4, endAnchored);", "        double(tr.x), double(tr.y), FS, res.f0, nHarm, 4, endAnchored);   % v025: nHarm (was 4)", 1};
E{end+1} = {"pack", "    r.M                 = res.M;", "    r.M                 = res.M;" + nl + "    if isfield(res, 'nHarm'), r.nHarm = res.nHarm; else, r.nHarm = NaN; end   % v025", 1};
E{end+1} = {"config", "substream cellKey*1e5 + bi*1e3 + rep"");   % v022", "substream cellKey*1e5 + bi*1e3 + rep"", 'NHARMONICS', NH_OPT);   % v022, v025", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV025:anchor", "%s: expected %d match(es), found %d. v024 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v024");  nMat = count(s, "_v024.mat");
if nName < 5 || nMat < 2, error("makeV025:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v024", "runLoopClosureFftnoise_v025");  s = replace(s, "_v024.mat", "_v025.mat");
s = replace(s, "fprintf('[v024] Output exists", "fprintf('[v025] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v025(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v025(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v025  v024 + NHarmonics option (template harmonic count; integer or 'fc'). GENERATED: see makeRunLoopClosureV025_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV025:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
