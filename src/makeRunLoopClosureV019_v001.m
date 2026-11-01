function outPath = makeRunLoopClosureV019_v001(opts)
%MAKERUNLOOPCLOSUREV019_V001  Derive runLoopClosureFftnoise_v019 from v017: the END-ANCHORED template as the production template (Sol's final review, SESSION_LOG Addendum 69).
%
%   WHY: the template zero-residual tail (Addenda 39, 47, 50) is an edge effect: the production template's regular windows can stop short of the last
%   sample, leaving a stretch whose residual is exactly 0. Omitting tail trials changes the population (Cook CTRL cccMed .776 -> .698), so the headline
%   quantities are recomputed with the tail REMOVED BY CONSTRUCTION: one extra final window ending at the last sample, s = N - win + 1
%   (templateEndAnchorTest_v001's variant, which bounded beta_gen* only; this run produces the full per-trial results, so cccMed and the tempo slope can be recomputed).
%   WHAT v019 CHANGES IN v017, AND ONLY THIS
%       template   templateSubtract_local: the window starts are 1:hop:(N-win+1) as before (identical to the while loop), plus one final window
%                  starting at N-win+1 when those stop short of the last sample.
%       names      function, error ids and output files carry v019 (loopClosureResults_<dataset>_all_shaped_xu_v019.mat); config gains END_ANCHORED_TEMPLATE = true.
%   At every defaults-setting v019 is v017 (and so v015's code path) except the template. Chunked, resumable, deadline-aware, per-trial seeded: as runH1HPC_v002.
%   Fail Loud: each replacement asserts its match count; any drift in v017 aborts before anything is written.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v017.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v019.m");
if ~isfile(inPath), error("makeV019:missingSource", "v017 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV019:exists", "v019 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = fileread(inPath);  nl = newline;
E = {};
E{end+1} = {"function declaration", "function results = runLoopClosureFftnoise_v017(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v019(datasetName, varargin)" + nl + "% runLoopClosureFftnoise_v019  v017 + END-ANCHORED template (one final window ending at the last sample). GENERATED: see makeRunLoopClosureV019_v001.m. Edit the generator, not this file." + nl, 1};
E{end+1} = {"error ids", "'runLoopClosureFftnoise_v017:", "'runLoopClosureFftnoise_v019:", 8};
E{end+1} = {"output filename literal", "'loopClosureResults_%s_%s_%s_v017.mat'", "'loopClosureResults_%s_%s_%s_v019.mat'", 2};
E{end+1} = {"option tag rename", """_v017.mat"", OPT_TAG + ""_v017.mat""", """_v019.mat"", OPT_TAG + ""_v019.mat""", 1};
E{end+1} = {"skip message", "fprintf('[v017] Output exists", "fprintf('[v019] Output exists", 1};
E{end+1} = {"banner", "'=== runLoopClosureFftnoise_v017: ", "'=== runLoopClosureFftnoise_v019: ", 1};
E{end+1} = {"log line", "'runLoopClosureFftnoise_v017  |", "'runLoopClosureFftnoise_v019  |", 1};
E{end+1} = {"template window starts", ...
    "    s = 1;" + nl + "    while s + win - 1 <= N" + nl + "        e   = s + win - 1;" + nl, ...
    "    starts = 1:hop:(N - win + 1);   % the production windows (identical to: s = 1; while s + win - 1 <= N, s = s + hop)" + nl + ...
    "    if win > N, starts = []; end" + nl + ...
    "    if ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window" + nl + ...
    "    for s = starts" + nl + "        e   = s + win - 1;" + nl, 1};
E{end+1} = {"template loop end", "        ww(s:e)   = ww(s:e) + han;" + nl + "        s         = s + hop;" + nl + "    end" + nl, ...
    "        ww(s:e)   = ww(s:e) + han;" + nl + "    end" + nl, 1};
E{end+1} = {"config flag", "'OPT_TAG', OPT_TAG, 'CHUNK_SIZE', CHUNK);   % v016, v017", "'OPT_TAG', OPT_TAG, 'CHUNK_SIZE', CHUNK, 'END_ANCHORED_TEMPLATE', true);   % v016, v017, v019", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV019:anchor", "%s: expected %d match(es), found %d. v017 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
fid = fopen(outPath, "w");  fwrite(fid, s);  fclose(fid);
fprintf("Wrote %s (%d replacements)\n", outPath, numel(E));
end
