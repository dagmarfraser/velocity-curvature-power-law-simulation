function outPath = makeRunLoopClosureV020_v001(opts)
%MAKERUNLOOPCLOSUREV020_V001  Derive runLoopClosureFftnoise_v020 from v019 for the N=48 Fraser corpus (2026-11-03).
%
%   WHY: the N=48 Fraser corpus (importDB_fraser_v002) has to be run under BOTH templates on BlueBEAR, the same
%   hardware as the other datasets' v015 and v019 runs (Dagmar, 2026-11-03): v015's production template stays the
%   reporting basis (Z12) and v019's end-anchored template the sensitivity analysis (D30).
%   WHAT v020 CHANGES IN v019, AND ONLY THIS
%     1. Option 'Template': 'endAnchored' (default, v019) or 'production' (v015's windows: the final end-anchored
%        window is NOT added). 'production' adds "_tplProd" to the output name; config gains TEMPLATE and
%        END_ANCHORED_TEMPLATE carries the actual setting.
%     2. Dataset 'FraserN48': importDB_fraser_v002 (gaps split, direct touches dropped, ID101/102/214 out, ID401-403
%        in), noise file noiseCharacterisation_fraserN48.mat (rerunNoiseCharacterisation_v006), and sigToMM =
%        1/5.209 (UIKit points per mm on the M4 13" iPad; the 'Fraser' case's 1/10.41793 is about 2x off and is
%        left as it was so v015/v019 stay reproducible). sigToMM reaches a_mm, b_mm and sigmaMM only, not beta.
%     3. Fail loud on the EBR fits: v019's three empty catch blocks around regressDataEBR (observed fit,
%        forward map, smoke check) now COUNT the error per trial and keep the last message
%        (r.nEbrErr, r.ebrErrMsg). The value written is unchanged (NaN), so at every v019 setting the numerics are
%        v019's. Non-convergence NaNs from regressDataEBR (its disp) are not counted here: that is the
%        registered limitBreak = 0 behaviour (auditFraserN48ObservedFits_v001 counts them for the observed fits).
%     4. Names: function, error ids and output files carry v020.
%   ACCEPTANCE: at Template='endAnchored' on dataset 'Fraser', v020's betaObs, betaGenStar and invertStatus must
%   equal v019's on the same machine and seed (checkV020Acceptance_v001).
%   Fail Loud: each replacement asserts its match count; any drift in v019 aborts before anything is written.
%   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v019.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v020.m");
if ~isfile(inPath), error("makeV020:missingSource", "v019 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV020:exists", "v020 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = fileread(inPath);  nl = newline;
E = {};
% 4. names
E{end+1} = {"function declaration", "function results = runLoopClosureFftnoise_v019(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v020(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v020  v019 + Template option (endAnchored | production), dataset FraserN48 (importDB_fraser_v002, 5.209 pt/mm), counted EBR failures. GENERATED: see makeRunLoopClosureV020_v001.m. Edit the generator, not this file." + nl, 1};
E{end+1} = {"error ids", "'runLoopClosureFftnoise_v019:", "'runLoopClosureFftnoise_v020:", 8};
E{end+1} = {"output filename literal", "'loopClosureResults_%s_%s_%s_v019.mat'", "'loopClosureResults_%s_%s_%s_v020.mat'", 2};
E{end+1} = {"option tag rename", """_v019.mat"", OPT_TAG + ""_v019.mat""", """_v020.mat"", OPT_TAG + ""_v020.mat""", 1};
E{end+1} = {"skip message", "fprintf('[v019] Output exists", "fprintf('[v020] Output exists", 1};
E{end+1} = {"banner", "'=== runLoopClosureFftnoise_v019: ", "'=== runLoopClosureFftnoise_v020: ", 1};
E{end+1} = {"log line", "'runLoopClosureFftnoise_v019  |", "'runLoopClosureFftnoise_v020  |", 1};
% 1. Template option
E{end+1} = {"Template option", "addParameter(p, 'OutTag'", ...
    "addParameter(p, 'Template',       'endAnchored', @(x) ismember(x, {'endAnchored','production'}));   % v020" + nl + "addParameter(p, 'OutTag'", 1};
E{end+1} = {"Template tag", "if EDGE_CLIP ~= 20, OPT_TAG", ...
    "END_ANCHORED = strcmp(opt.Template, 'endAnchored');   % v020" + nl + ...
    "if ~END_ANCHORED, OPT_TAG = OPT_TAG + ""_tplProd""; end   % v020: v015's production windows" + nl + ...
    "if EDGE_CLIP ~= 20, OPT_TAG", 1};
E{end+1} = {"config flag", "'END_ANCHORED_TEMPLATE', true);   % v016, v017, v019", ...
    "'END_ANCHORED_TEMPLATE', END_ANCHORED, 'TEMPLATE', opt.Template, 'SIG_TO_MM', sigToMM);   % v016, v017, v019, v020", 1};
E{end+1} = {"processOneTrial calls", "nDeriv, nRegress, nPipelines, NOISE_MODEL);", "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED);", 4};
E{end+1} = {"processOneTrial signature", "nRegress, nPipelines, noiseModel)" + nl, "nRegress, nPipelines, noiseModel, endAnchored)" + nl, 1};
E{end+1} = {"template call", "res.f0, 4, 4);", "res.f0, 4, 4, endAnchored);", 1};
E{end+1} = {"template signature", "templateSubtract_local(x, y, fs, f0, nH, nCW)", "templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)", 1};
E{end+1} = {"end-anchor switch", "    if ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window", ...
    "    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window; v020: only when Template = 'endAnchored'", 1};
% 2. FraserN48
E{end+1} = {"FraserN48 dataset", "    otherwise" + nl, ...
    "    case 'FraserN48'   % v020" + nl + ...
    "        matFile  = fullfile(srcDir, 'noiseCharacterisation_fraserN48.mat');" + nl + ...
    "        sigToMM  = 1 / 5.209;   % UIKit points per mm, M4 13"" iPad (Step 1 audit, 2026-11-03)" + nl + ...
    "        FS_DS    = 240;" + nl + ...
    "        importFn = @() importDB_fraser_v002(Verbose=false);" + nl + ...
    "    otherwise" + nl, 1};
% 3. counted EBR failures
E{end+1} = {"EBR counter init", "    doR7    = strcmp(noiseModel, 'shaped_xu');" + nl, ...
    "    doR7    = strcmp(noiseModel, 'shaped_xu');" + nl + ...
    "    res.nEbrErr = 0;  res.ebrErrMsg = """";   % v020: regressDataEBR errors, counted (were silent)" + nl, 1};
E{end+1} = {"EBR counter pack", "    r.smokeA            = res.smokeA;" + nl, ...
    "    r.smokeA            = res.smokeA;" + nl + "    r.nEbrErr           = res.nEbrErr;     % v020" + nl + "    r.ebrErrMsg         = res.ebrErrMsg;   % v020" + nl, 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV020:anchor", "%s: expected %d match(es), found %d. v019 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
% 3. the three empty catches (observed fit, forward map, smoke check): bare 'catch' line followed by 'end'
pat = "(\n)( +)catch *\n( +)end\n";
nCatch = numel(regexp(s, pat, "match"));
if nCatch ~= 3, error("makeV020:anchor", "empty catch blocks: expected 3, found %d. v019 has drifted; nothing written.", nCatch); end
s = regexprep(s, pat, "$1$2catch ME   % v020: counted, not silent" + nl + ...
    "$2    res.nEbrErr = res.nEbrErr + 1;  res.ebrErrMsg = string(ME.identifier) + "": "" + string(ME.message);" + nl + "$3end" + nl);
% Fail loud: nothing of v019's identity may survive outside comments
if contains(s, "_v019.mat") || contains(s, "'runLoopClosureFftnoise_v019")
    error("makeV020:residue", "%s", "v019 names survive in the generated code; nothing written.");
end
fid = fopen(outPath, "w");  fwrite(fid, s);  fclose(fid);
fprintf("Wrote %s (%d replacements + %d catch blocks)\n", outPath, numel(E), nCatch);
end
