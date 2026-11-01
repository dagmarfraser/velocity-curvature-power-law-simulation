function outPath = makeRunLoopClosureV016_v001(opts)
%MAKERUNLOOPCLOSUREV016_V001  Derive runLoopClosureFftnoise_v016 from v015 (D-D16 stage S2).
%
%   WHY A GENERATOR: as makeRunLoopClosureV015_v001 (anchored, asserted replacements cannot drift; a
%   retyped 1,080-line runner could). Edit this generator, never the generated file.
%
%   WHAT v016 ADDS (docs/SPEC_D16_SamplingRateCutoff_v001.md, section 3, S2), AND ONLY THIS
%   Three options, with v015's values as defaults:
%       'Subsample' 1   integer stride; x and y of every selected trial are subsampled cleanly, x(1:Subsample:end),
%                       no filter (Dagmar, 2026-10-24 13:04: "we subsample cleanly"), before anything else, and the
%                       runner then works at FS_RUN = FS_DS / Subsample (beta_obs, f0, template, IRASA, forward map,
%                       inversion). (Drafted first with decimate(); changed before this generator was ever run.)
%       'BWCutoff'  10  Butterworth cutoff in Hz (filterType 2, [2 BWCutoff 1])
%       'SGRef'     17  Savitzky-Golay reference frame at 100 Hz (filterType 6, [4 SGRef])
%   'EdgeClip' already exists in v015 (default 20 samples). Any non-default option is added to the output
%   and log names (e.g. loopClosureResults_Fraser_all_shaped_xu_sub2_v016.mat); all four are saved in config.
%   At the defaults v016 runs v015's code path exactly: no subsampling, the same filter literals, the same name
%   pattern with v016 in place of v015.
%
%   NOT CHANGED, AND THEREFORE NOT RATE-MATCHED: r.alphaIRA and r.sigma come from the 240 Hz noise file
%   (bioResults), as in v015; they only label the output. The forward map's own noise (alphaMaj/alphaMin,
%   the residual) is re-estimated at FS_RUN inside processOneTrial_local, as in v015.
%
%   ACCEPTANCE TEST (runD16HPC_v001 runs it before anything else): v016 at the defaults on Zarandi must
%   reproduce v015's betaObs exactly and its stochastic fields within runLoopClosureAll_v015's tolerances
%   (Finding #204: beta_gen* is not run-reproducible at a fixed seed across parallel configurations).
%
%   Fail Loud: every replacement asserts its match count; any drift in v015 aborts before anything is written.
%   Provenance: written 2026-10-24 (Session 146), Dagmar's decisions of that day (both clip arms; the sweep).

arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath  = fullfile(srcDir, "runLoopClosureFftnoise_v015.m");
outPath = fullfile(srcDir, "runLoopClosureFftnoise_v016.m");
if ~isfile(inPath), error("makeV016:missingSource", "v015 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite
    error("makeV016:exists", "v016 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath);
end
s = fileread(inPath);  nl = newline;
J = @(c) strjoin(c, nl) + nl;          % join lines (not sprintf: "%" in MATLAB source must stay literal)

% ---- anchored replacements: {description, old, new, expectedCount} -------
E = {};
E{end+1} = {"function declaration", ...
    "function results = runLoopClosureFftnoise_v015(datasetName, varargin)", ...
    "function results = runLoopClosureFftnoise_v016(datasetName, varargin)", 1};
E{end+1} = {"error ids", "'runLoopClosureFftnoise_v015:", "'runLoopClosureFftnoise_v016:", 4};   % leading quote: the banner also contains the name and a colon
E{end+1} = {"output filename literal", "'loopClosureResults_%s_%s_%s_v015.mat'", "'loopClosureResults_%s_%s_%s_v016.mat'", 2};
E{end+1} = {"skip message", "fprintf('[v015] Output exists", "fprintf('[v016] Output exists", 1};
E{end+1} = {"banner", ...
    "'=== runLoopClosureFftnoise_v015: %s | %s | model=%s | parfor=%d ===\n'", ...
    "'=== runLoopClosureFftnoise_v016: %s | %s | model=%s | parfor=%d ===\n'", 1};
E{end+1} = {"log line", ...
    "'runLoopClosureFftnoise_v015  |  %s  |  %s  |  model=%s  |  parfor=%d\n'", ...
    "'runLoopClosureFftnoise_v016  |  %s  |  %s  |  model=%s  |  parfor=%d\n'", 1};
E{end+1} = {"new parameters", ...
    "addParameter(p, 'EdgeClip',       20,    @isscalar);" + nl, J([
    "addParameter(p, 'EdgeClip',       20,    @isscalar);"
    "addParameter(p, 'Subsample',      1,     @(x) isscalar(x) && x >= 1 && x == round(x));   % v016: clean stride, no filter"
    "addParameter(p, 'BWCutoff',       10,    @(x) isscalar(x) && x > 0);   % v016: Butterworth cutoff, Hz"
    "addParameter(p, 'SGRef',          17,    @(x) isscalar(x) && x > 0);   % v016: SG reference frame at 100 Hz"]), 1};
E{end+1} = {"new option variables", ...
    "EDGE_CLIP       = opt.EdgeClip;" + nl, J([
    "EDGE_CLIP       = opt.EdgeClip;"
    "SUBSAMPLE       = opt.Subsample;  % v016"
    "BW_CUTOFF       = opt.BWCutoff;   % v016"
    "SG_REF          = opt.SGRef;      % v016"]), 1};
E{end+1} = {"filter parameters", "'filterParams', {[2 10 1],[4 17]});", "'filterParams', {[2 BW_CUTOFF 1],[4 SG_REF]});   % v016: v015 literals were [2 10 1],[4 17]", 1};
E{end+1} = {"run rate and option tag", ...
    "%% ---------- Load noise characterisation", J([
    "%% ---------- v016: run rate and option tag ---------------------------------"
    "FS_RUN  = FS_DS / SUBSAMPLE;"
    "OPT_TAG = """";"
    "if SUBSAMPLE ~= 1,  OPT_TAG = OPT_TAG + ""_sub"" + SUBSAMPLE; end"
    "if BW_CUTOFF ~= 10, OPT_TAG = OPT_TAG + ""_bw"" + strrep(string(BW_CUTOFF), ""."", ""p""); end"
    "if SG_REF ~= 17,    OPT_TAG = OPT_TAG + ""_sg"" + strrep(string(SG_REF), ""."", ""p""); end"
    "if EDGE_CLIP ~= 20, OPT_TAG = OPT_TAG + ""_clip"" + EDGE_CLIP; end"
    ""]) + "%% ---------- Load noise characterisation", 1};
E{end+1} = {"tag the output name", ...
    "if isfile(outFile) && ~opt.ForceRerun", J([
    "outFile = char(replace(string(outFile), ""_v016.mat"", OPT_TAG + ""_v016.mat""));   % v016: no-op at the defaults"
    ""]) + "if isfile(outFile) && ~opt.ForceRerun", 1};
E{end+1} = {"report options", ...
    "fprintf('Importing trials...\n');", J([
    "fprintf('  v016 options: subsample stride %d (run at %g Hz), BW cutoff %g Hz, SG reference frame %g, edge clip %d, tag ""%s""\n', ..."
    "    SUBSAMPLE, FS_RUN, BW_CUTOFF, SG_REF, EDGE_CLIP, OPT_TAG);"]) + "fprintf('Importing trials...\n');", 1};
E{end+1} = {"subsample the selected trials", ...
    "trialsSel = trials(matchIdx);" + nl + "clear trials" + nl, J([
    "trialsSel = trials(matchIdx);"
    "clear trials"
    "if SUBSAMPLE > 1   % v016: clean subsampling, every SUBSAMPLE-th sample from the first, no filter"
    "    for k = 1:numel(trialsSel)"
    "        trialsSel(k).x = trialsSel(k).x(1:SUBSAMPLE:end);"
    "        trialsSel(k).y = trialsSel(k).y(1:SUBSAMPLE:end);"
    "        if isfield(trialsSel, 't'),  trialsSel(k).t  = (0:numel(trialsSel(k).x) - 1)' / FS_RUN; end"
    "        if isfield(trialsSel, 'fs'), trialsSel(k).fs = FS_RUN; end"
    "    end"
    "    fprintf('  v016: %d trials subsampled (stride %d, no filter) to %g Hz.\n', numel(trialsSel), SUBSAMPLE, FS_RUN);"
    "end"]), 1};
E{end+1} = {"run rate into the trial processing", "FS_DS, sigToMM, betaGenVec", "FS_RUN, sigToMM, betaGenVec", 3};
E{end+1} = {"tag the log name", ...
    "    tag, TRIAL_SELECTION, datestr(now, 'yyyymmdd_HHMMSS')));", ...
    "    [tag char(OPT_TAG)], TRIAL_SELECTION, datestr(now, 'yyyymmdd_HHMMSS')));", 1};
E{end+1} = {"options into config", "'monoParams', monoParams);", ...
    "'monoParams', monoParams, 'SUBSAMPLE', SUBSAMPLE, 'FS_DS', FS_DS, 'FS_RUN', FS_RUN, ..." + nl + ...
    "    'BW_CUTOFF', BW_CUTOFF, 'SG_REF', SG_REF, 'EDGE_CLIP', EDGE_CLIP, 'OPT_TAG', OPT_TAG);   % v016", 1};
E{end+1} = {"header note", ...
    "function results = runLoopClosureFftnoise_v016(datasetName, varargin)" + nl, J([
    "function results = runLoopClosureFftnoise_v016(datasetName, varargin)"
    "% runLoopClosureFftnoise_v016  v015 + options for sampling rate and filters (D-D16 stage S2)."
    "%"
    "% GENERATED, not hand-written: see makeRunLoopClosureV016_v001.m. Edit the generator, not this file."
    "%"
    "% THE ONLY CHANGES FROM v015: options 'Subsample' (1), 'BWCutoff' (10 Hz) and 'SGRef' (17), beside v015's"
    "% 'EdgeClip' (20). Clean subsampling (x(1:Subsample:end), no filter) is applied to the selected trials and"
    "% everything then runs at FS_RUN = FS_DS / Subsample. Non-default options are added to the output and log names and all"
    "% are saved in config. r.alphaIRA and r.sigma still come from the 240 Hz noise file (labels only); the"
    "% forward map re-estimates its noise at FS_RUN, as v015 does at FS_DS."
    "% ACCEPTANCE TEST: at the defaults v016 is v015 (runD16HPC_v001 checks it on Zarandi first)."
    "% See docs/SPEC_D16_SamplingRateCutoff_v001.md. v015's own header follows unchanged."]), 1};

for i = 1:numel(E)
    [desc, old, new, want] = E{i}{:};
    got = count(s, old);
    if got ~= want
        error("makeV016:anchor", "Anchor '%s': expected %d match(es), found %d. Nothing written: " + ...
            "v015 has drifted from what this generator expects.", desc, want, got);
    end
    s = replace(s, old, new);
end
if contains(s, "FS_DS, sigToMM") || ~contains(s, "FS_RUN, sigToMM")
    error("makeV016:check", "%s", "Trial processing still receives FS_DS; nothing written.");
end

fid = fopen(outPath, "w");  if fid < 0, error("makeV016:write", "Could not open for writing:\n  %s", outPath); end
fwrite(fid, s);  fclose(fid);
fprintf("Wrote %s\n  (%d anchored edits, %d bytes)\n", outPath, numel(E), strlength(s));
end
