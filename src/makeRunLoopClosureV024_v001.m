function outPath = makeRunLoopClosureV024_v001(opts)
%MAKERUNLOOPCLOSUREV024_V001  Derive runLoopClosureFftnoise_v024 from v023: keep the template's uncovered tail out of the noise.
%
%   WHY (coherence check 2026-11-15, checkTemplateTailSegment_v001): the production template lays its overlap-add windows from
%   the start of the trial, so the kept residual can end in a stretch no window reached: exactly zero, after a step (K-MF5).
%   That stretch occurs in 22% of Fraser, 30-33% of Hickman and 42-45% of Cook trials (median 4-9% of the residual), and v023's
%   end-matched segment includes it in 95-98% of them, so the step's broadband energy enters the noise of 'endMatched' and
%   'shapedXuSeg'. The stretch is not residual (no template was subtracted there); it is excluded.
%   WHAT v024 CHANGES IN v023, AND ONLY THIS
%     1. templateSubtract_local also returns tailN, the number of uncovered samples at the end of the kept residual
%        (0 whenever the last window reaches the end clip, e.g. always with Template 'endAnchored'). No numeric change.
%     2. For 'endMatched' and 'shapedXuSeg' only: the segment search runs on the covered residual, resid(1:M - tailN), and the
%        central beta_obs window is centred on that covered span. New per-trial field emTail (= tailN). Trials without a
%        tail get exactly v023's segments, windows and noise.
%     3. Names carry v024. Defaults ('shaped_xu', full residual) are untouched, so v024 reproduces v023 there.
%   ACCEPTANCE: checkV024Acceptance_v001.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v023.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v024.m");
if ~isfile(inPath), error("makeV024:missingSource", "v023 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV024:exists", "v024 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;
E = {};
E{end+1} = {"template def", "function [resX, resY, fitX, fitY, cs] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)", ...
    "function [resX, resY, fitX, fitY, cs, tailN] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)   % v024: + tailN", 1};
E{end+1} = {"template tail", "    ce   = min(max(N - cl, cs + 100), N);", ...
    "    ce   = min(max(N - cl, cs + 100), N);" + nl + ...
    "    lastEnd = 0;  if ~isempty(starts), lastEnd = starts(end) + win - 1; end   % v024: end of the last overlap-add window" + nl + ...
    "    tailN = max(0, ce - lastEnd);   % v024: samples at the end of the kept residual that no window reached (exactly zero)", 1};
E{end+1} = {"template call", "[xResid, yResid, xFit, yFit, tsClip] = templateSubtract_local( ...", "[xResid, yResid, xFit, yFit, tsClip, tsTail] = templateSubtract_local( ...   % v024: + tail", 1};
E{end+1} = {"segment search", "[sI0, sI1, emInfo] = endToEndSegment_v004([resid_maj(:), resid_min(:)], Tol=1, Mode=""commonLength"");", ...
    "Mc = M - tsTail;  res.emTail = tsTail;   % v024: search only the covered residual" + nl + ...
    "            [sI0, sI1, emInfo] = endToEndSegment_v004([resid_maj(1:Mc), resid_min(1:Mc)], Tol=1, Mode=""commonLength"");", 1};
E{end+1} = {"window centre", "            c0 = floor((M - Lsyn) / 2) + 1;", "            c0 = floor((Mc - Lsyn) / 2) + 1;   % v024: centred on the covered span", 1};
E{end+1} = {"emTail field", "'emWindow', [NaN NaN];", "'emWindow', [NaN NaN]; 'emTail', NaN;", 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV024:anchor", "%s: expected %d match(es), found %d. v023 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v023");  nMat = count(s, "_v023.mat");
if nName < 5 || nMat < 2, error("makeV024:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v023", "runLoopClosureFftnoise_v024");  s = replace(s, "_v023.mat", "_v024.mat");
s = replace(s, "fprintf('[v023] Output exists", "fprintf('[v024] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v024(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v024(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v024  v023 + the production template's uncovered (zero) tail excluded from the end-matched segment search. GENERATED: see makeRunLoopClosureV024_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV024:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
