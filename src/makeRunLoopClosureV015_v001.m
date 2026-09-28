function outPath = makeRunLoopClosureV015_v001(opts)
%MAKERUNLOOPCLOSUREV015_V001  Derive runLoopClosureFftnoise_v015 from v013.
%
%   WHY A GENERATOR RATHER THAN A HAND-WRITTEN v015
%   v013 is 1048 lines, most of it a carefully-kept historical header. The
%   change v015 needs is two lines. Retyping the file risks silent
%   transcription error in code nobody would re-read; a generator with
%   anchored, asserted replacements cannot drift. This also keeps the
%   project's duplication-over-refactor convention (v015 is a standalone
%   file, v013's outputs are untouched) without the transcription risk.
%
%   WHAT v015 CHANGES, AND ONLY THIS
%   v013 computes the per-trial forward map -- the beta_gen -> beta_rec
%   curve the inversion is performed against -- and then DISCARDS it,
%   persisting only the inversion's outputs plus betaRecSlice, which is a
%   [N_REPS x 6] replicate cloud at one grid point, NOT the curve.
%   Without the curve there is no way to predict a held-out pipeline's
%   beta_obs, so the strong (out-of-sample) form of leave-one-pipeline-out
%   is impossible against existing corpora. See Finding #199 (scoped form),
%   docs/CritiqueResponse_ResearchTrajectory_v001.md rows A2/A3.
%
%   v015 persists exactly the two reductions the inversion itself consumes
%   (v013 lines 582-583):
%       betaRecCurveMed = squeeze(median(res.betaRec, 1, 'omitnan'))  [6 x N_BETA]
%       betaRecCurveSD  = squeeze(std(res.betaRec, 0, 1, 'omitnan'))  [6 x N_BETA]
%   Using the same statistics means the out-of-sample prediction is made
%   against the identical curve the inversion used, not a reconstruction
%   of it. Cost: ~2.4 kB/trial (~7 MB for Fraser).
%
%   NOTHING ELSE CHANGES. betaGenStar, invertStatus, the gate
%   (findBothBranches_v008), CI handling and all numerics are v013's,
%   untouched. v015 must therefore reproduce v013's betaGenStar exactly;
%   that is the acceptance test, and runLoopClosureAll_v015.m checks it.
%
%   Fail Loud, Never Fake: every replacement asserts its expected match
%   count. Any drift in v013 aborts before anything is written.
%
%   Provenance: written 2026-09-16. Dagmar's decision: go straight to HPC
%   rather than pilot locally (local estimate ~59 h on 10 workers vs ~8 h
%   on 72; see SESSION_LOG Session 110).

arguments
    opts.SrcDir   (1,1) string = ""
    opts.Overwrite (1,1) logical = false
end

if opts.SrcDir == ""
    srcDir = fileparts(mfilename("fullpath"));
else
    srcDir = opts.SrcDir;
end

inPath  = fullfile(srcDir, "runLoopClosureFftnoise_v013.m");
outPath = fullfile(srcDir, "runLoopClosureFftnoise_v015.m");

if ~isfile(inPath)
    error("makeV015:missingSource", "v013 not found:\n  %s", inPath);
end
if isfile(outPath) && ~opts.Overwrite
    error("makeV015:exists", ...
        "v015 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath);
end

s = fileread(inPath);

% ---- anchored replacements: {description, old, new, expectedCount} -------
E = {};

E{end+1} = {"function declaration", ...
    "function results = runLoopClosureFftnoise_v013(datasetName, varargin)", ...
    "function results = runLoopClosureFftnoise_v015(datasetName, varargin)", 1};

E{end+1} = {"output filename literal", ...
    "'loopClosureResults_%s_%s_%s_v013.mat'", ...
    "'loopClosureResults_%s_%s_%s_v015.mat'", 2};

E{end+1} = {"skip message", ...
    "fprintf('[v013] Output exists", ...
    "fprintf('[v015] Output exists", 1};

E{end+1} = {"banner", ...
    "'=== runLoopClosureFftnoise_v013: %s | %s | model=%s | parfor=%d ===\n'", ...
    "'=== runLoopClosureFftnoise_v015: %s | %s | model=%s | parfor=%d ===\n'", 1};

E{end+1} = {"log line", ...
    "'runLoopClosureFftnoise_v013  |  %s  |  %s  |  model=%s  |  parfor=%d\n'", ...
    "'runLoopClosureFftnoise_v015  |  %s  |  %s  |  model=%s  |  parfor=%d\n'", 1};

for id = ["UnknownDataset","NoIRASA","TrialNotFound","NoMatch"]
    E{end+1} = {"error id " + id, ...
        "runLoopClosureFftnoise_v013:" + id, ...
        "runLoopClosureFftnoise_v015:" + id, 1}; %#ok<AGROW>
end

% the one substantive change.
% NOTE: built with strjoin(...), NOT sprintf. sprintf would reinterpret the
% "%" comment markers in the generated MATLAB source as format specifiers, and
% ["a" "b"] builds a 1x2 string ARRAY rather than concatenating. Both bit on
% first run, 2026-09-16.
nl = newline;
E{end+1} = {"persist forward map", ...
    "    r.betaRecSlice      = betaRecSlice;" + nl, ...
    strjoin([ ...
    "    r.betaRecSlice      = betaRecSlice;"
    "    % v015: persist the forward map the inversion is performed against."
    "    % These are the SAME two reductions invertBeta_local consumes (median"
    "    % and std over reps), so any out-of-sample prediction made from them"
    "    % uses the identical curve the inversion used, not a reconstruction."
    "    % [nPipelines x N_BETA]. v013 computed and discarded these."
    "    r.betaRecCurveMed   = squeeze(median(res.betaRec, 1, 'omitnan'));"
    "    r.betaRecCurveSD    = squeeze(std(res.betaRec, 0, 1, 'omitnan'));"], nl) + nl, 1};

% ---- header note --------------------------------------------------------
E{end+1} = {"header note", ...
    "function results = runLoopClosureFftnoise_v015(datasetName, varargin)" + nl, ...
    strjoin([ ...
    "function results = runLoopClosureFftnoise_v015(datasetName, varargin)"
    "% runLoopClosureFftnoise_v015  v013 + persist the per-trial forward map."
    "%"
    "% GENERATED, not hand-written: see makeRunLoopClosureV015_v001.m. Edit the"
    "% generator, not this file, or the next regeneration will overwrite you."
    "%"
    "% THE ONLY CHANGE FROM v013 is in packResult_local: two extra fields,"
    "%   betaRecCurveMed  [nPipelines x N_BETA]  median over reps"
    "%   betaRecCurveSD   [nPipelines x N_BETA]  std over reps"
    "% These are the exact two reductions invertBeta_local already computes"
    "% (v013 lines 582-583) and then throws away. v013 persisted only"
    "% betaRecSlice, a [N_REPS x 6] replicate cloud at ONE grid point, which"
    "% is not a forward map and cannot predict a held-out pipeline's betaObs."
    "%"
    "% WHY: the strong, out-of-sample form of leave-one-pipeline-out -- infer"
    "% beta_gen* from five pipelines, predict the sixth's betaObs through its"
    "% own forward map -- is the residue of GPT SOL's circularity critique"
    "% (section 2) and the manuscript's own most substantive open item. The"
    "% scoped form (Finding #199) compares inverted values only and cannot"
    "% answer it."
    "%"
    "% ACCEPTANCE TEST: v015 changes no numerics, so betaGenStar, invertStatus,"
    "% ciLo/ciHi and betaObs must match v013's EXACTLY on any dataset run with"
    "% the same seed and config. runLoopClosureAll_v015.m checks this and"
    "% refuses to proceed if it fails."
    "%"
    "% v013's own header follows unchanged."], nl) + nl, 1};

for i = 1:numel(E)
    [desc, old, new, want] = E{i}{:};
    got = count(s, old);
    if got ~= want
        error("makeV015:anchor", ...
            "Anchor '%s': expected %d match(es), found %d. Nothing written -- " + ...
            "v013 has drifted from what this generator expects.", desc, want, got);
    end
    s = replace(s, old, new);
end

fid = fopen(outPath, "w");
if fid < 0
    error("makeV015:write", "Could not open for writing:\n  %s", outPath);
end
fwrite(fid, s);
fclose(fid);

fprintf("Wrote %s\n  (%d anchored edits, %d bytes)\n", outPath, numel(E), numel(s));
end
