function outPath = makeRunLoopClosureV022_v001(opts)
%MAKERUNLOOPCLOSUREV022_V001  Derive runLoopClosureFftnoise_v022 from v021: the end-matched noise model (decisions of 2026-11-14/15).
%
%   WHY: Dagmar chose (2026-11-14, option 3b) the end-matched random-phase noise as the primary noise model, with the
%   segment length from the data (fairShare rule), the synthetic trajectory shortened to that length ("do it", 22:57), and
%   beta_obs computed both on a central window of the real trial with the synthetic's length (primary, like with like) and on
%   the whole trial (sensitivity, the registered basis) (2026-11-15, 00:46). Evidence:
%   docs/EndMismatch_FairShareRule_FullLength_20261114_v001.md, docs/ShortenedSynthetic_SegmentAndBetaObs_20261114_v001.md.
%   WHAT v022 CHANGES IN v021, AND ONLY THIS
%     1. NoiseModel 'endMatched' (new; the default stays 'shaped_xu', so v022 at its defaults is v021). Per trial:
%        endToEndSegment_v004([resid_maj, resid_min], Tol 1, Mode "commonLength") gives one length L and a segment per axis;
%        L below the runner minimum 4*(EDGE_CLIP + 50) or a non-finite segment alpha (alphaPWLW_v001, [1, min(20, fs/2 - 1)])
%        makes the trial NaN with the reason in emStatus (Fail Loud). Templates are built at length L. Noise per replicate is
%        endMatchedNoise_v001 (random phases on each axis's own segment, alpha selection Tol 0.3, 20 draws); a replicate with
%        no accepted draw is skipped by the existing NaN check and counted in nNoiseRejected.
%     2. Seeds for 'endMatched' follow docs/README_SeedPractice_v001.md: threefry4x64_20, master RngSeed, substream
%        cellKey*1e5 + betaIndex*1e3 + rep, where cellKey is the trial's bioResults row (a fixed identity, not its position in
%        the run). Results therefore do not depend on parfor, chunking or the trial list. Asserts N_REPS < 1000, N_BETA < 100.
%     3. beta_obs for 'endMatched': betaObs (and vgfObs) are on the central window of the real trial with length L (trial
%        index = residual index + template clip - 1); betaObsWhole/vgfObsWhole keep the whole-trial values, and the trial is
%        inverted against the same forward map a second time from them (betaGenStarWhole, betaGenStarMedWhole, ciLoWhole,
%        ciHiWhole, invertStatusWhole). smokeA (the fftnoise determinism check) is NaN for 'endMatched'.
%     4. templateSubtract_local also returns its clip cs (no numeric change). trialsSel gains cellKey. processOneTrial_local
%        takes the master seed. Names: function, error ids, output files and messages carry v022.
%   ACCEPTANCE: at the defaults v022 must equal v021 exactly (isequaln of betaObs, betaGenStar, invertStatus) on the same
%   machine and seed (checkV022Acceptance_v001). Fail Loud: each replacement asserts its match count.   Fraser, D.S. (2026)  v001
arguments
    opts.SrcDir    (1,1) string  = ""
    opts.Overwrite (1,1) logical = false
end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fileparts(mfilename("fullpath")); end
inPath = fullfile(srcDir, "runLoopClosureFftnoise_v021.m");  outPath = fullfile(srcDir, "runLoopClosureFftnoise_v022.m");
if ~isfile(inPath), error("makeV022:missingSource", "v021 not found:\n  %s", inPath); end
if isfile(outPath) && ~opts.Overwrite, error("makeV022:exists", "v022 already exists:\n  %s\nPass Overwrite=true to regenerate.", outPath); end
s = string(fileread(inPath));  nl = newline;  J = @(L) strjoin(L, nl);
E = {};
E{end+1} = {"NoiseModel option", "@(x) ismember(x, {'auto','fftnoise','xu','shaped_xu','bootstrap'}));", ...
    "@(x) ismember(x, {'auto','fftnoise','xu','shaped_xu','bootstrap','endMatched'}));   % v022: endMatched", 1};
E{end+1} = {"cellKey", "trialsSel = trials(matchIdx);" + nl, ...
    "trialsSel = trials(matchIdx);" + nl + "for k = 1:numel(trialsSel), trialsSel(k).cellKey = selIdx(k); end   % v022: fixed identity for the endMatched substreams" + nl, 1};
E{end+1} = {"processOneTrial calls", "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA);", ...
    "nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED);", 4};
E{end+1} = {"processOneTrial signature", "nRegress, nPipelines, noiseModel, endAnchored, carrierAlpha)" + nl, ...
    "nRegress, nPipelines, noiseModel, endAnchored, carrierAlpha, rngSeed)" + nl, 1};
E{end+1} = {"template clip call", "    [xResid, yResid, xFit, yFit] = templateSubtract_local( ...", "    [xResid, yResid, xFit, yFit, tsClip] = templateSubtract_local( ...   % v022: + clip", 1};
E{end+1} = {"template clip def", "function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)", ...
    "function [resX, resY, fitX, fitY, cs] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)   % v022: + cs", 1};
emBlock = J([ ...
    "    MTPL = M;   % v022: template length (the whole residual unless endMatched)"
    "    if strcmp(noiseModel, 'endMatched')   % v022"
    "        res.emL = NaN;  res.emLfrac = NaN;  res.emSeg = NaN(2, 2);  res.emAlphaSeg = [NaN NaN];  res.emStatus = ""ok"";  res.nNoiseRejected = 0;"
    "        res.betaObsWhole_vec = res.betaObs_vec;  res.vgfObsWhole_vec = res.vgfObs_vec;"
    "        try"
    "            if N_REPS >= 1000 || N_BETA >= 100, error('runLoopClosureFftnoise_v022:emSubstream', 'substream layout needs N_REPS < 1000 and N_BETA < 100'); end"
    "            [sI0, sI1, emInfo] = endToEndSegment_v004([resid_maj(:), resid_min(:)], Tol=1, Mode=""commonLength"");"
    "            Lsyn = emInfo.L;"
    "            if Lsyn < 4*(EDGE_CLIP + 50), error('runLoopClosureFftnoise_v022:emShort', 'end-matched length %d is below the runner minimum %d', Lsyn, 4*(EDGE_CLIP + 50)); end"
    "            segMaj = resid_maj(sI0(1):sI1(1));  segMin = resid_min(sI0(2):sI1(2));  segMaj = segMaj(:);  segMin = segMin(:);"
    "            aSeg = [alphaPWLW_v001(segMaj, FS, [1.0, min(20.0, FS/2 - 1)]), alphaPWLW_v001(segMin, FS, [1.0, min(20.0, FS/2 - 1)])];"
    "            if ~all(isfinite(aSeg)), error('runLoopClosureFftnoise_v022:emAlpha', 'segment alpha not finite (%g, %g)', aSeg(1), aSeg(2)); end"
    "            c0 = floor((M - Lsyn) / 2) + 1;  wA = c0 + tsClip - 1;  wB = wA + Lsyn - 1;   % central window of the real trial, length Lsyn"
    "            [res.betaObs_vec, res.vgfObs_vec, nErrW, msgW] = betaObsOn_local(double(tr.x(wA:wB)), double(tr.y(wA:wB)), FS, EDGE_CLIP, derivConfigs, regressConfigs, nDeriv, nRegress);"
    "            res.nEbrErr = res.nEbrErr + nErrW;  if nErrW > 0, res.ebrErrMsg = msgW; end"
    "            res.emL = Lsyn;  res.emLfrac = emInfo.Lfrac;  res.emSeg = [sI0(:), sI1(:)];  res.emAlphaSeg = aSeg;  res.emWindow = [wA wB];"
    "            MTPL = Lsyn;"
    "        catch ME"
    "            res.emStatus = ""FAIL: "" + string(ME.message);"
    "            res.betaObs_vec(:) = NaN;  res.vgfObs_vec(:) = NaN;"
    "            res.betaRec = nanTemplate;  res.smokeA = nanSlimVec;"
    "            return"
    "        end"
    "    end"
    ""
    "    nPhi  = 10000;"]);
E{end+1} = {"endMatched block", "    nPhi  = 10000;", emBlock, 1};
E{end+1} = {"template time", "    tVec  = (0:M-1)' / FS;", "    tVec  = (0:MTPL-1)' / FS;   % v022: MTPL = M unless endMatched", 1};
E{end+1} = {"template size", "    templates = zeros(M, 2, N_BETA);", "    templates = zeros(MTPL, 2, N_BETA);   % v022", 1};
caseBlock = J([ ...
    "                case 'endMatched'   % v022: per-cell substream (README_SeedPractice), each axis's own end-matched segment"
    "                    rsC = RandStream(""threefry4x64_20"", ""Seed"", rngSeed);  rsC.Substream = tr.cellKey * 1e5 + bi * 1e3 + rep;"
    "                    rsOld = RandStream.setGlobalStream(rsC);"
    "                    [nMaj, nMin] = endMatchedNoise_v001(segMaj, segMin, FS, aSeg(1), aSeg(2));"
    "                    RandStream.setGlobalStream(rsOld);"
    "                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end"
    "                case 'bootstrap'"]);
E{end+1} = {"noise case", "                case 'bootstrap'", caseBlock, 1};
E{end+1} = {"smokeA", "    res.smokeA = smokeA;", "    if strcmp(noiseModel, 'endMatched'), smokeA(:) = NaN; end   % v022: the determinism check is whole-trial fftnoise; not comparable" + nl + "    res.smokeA = smokeA;", 1};
E{end+1} = {"whole-trial inversion (chunk, parfor)", "invertInfo); %#ok<PFBNS>", ...
    "invertInfo); %#ok<PFBNS>" + nl + "            r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, NOISE_MODEL);   % v022", 2};
E{end+1} = {"whole-trial inversion (serial)", "sigmaSelected(s), NOISE_MODEL, invertInfo);" + nl, ...
    "sigmaSelected(s), NOISE_MODEL, invertInfo);" + nl + "        r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, NOISE_MODEL);   % v022" + nl, 1};
E{end+1} = {"config", "'CARRIER_ALPHA', CARRIER_ALPHA);", ...
    "'CARRIER_ALPHA', CARRIER_ALPHA, 'V022_ENDMATCHED', ""endToEndSegment_v004 commonLength Tol 1; endMatchedNoise_v001 Tol 0.3 MaxDraw 20; betaObs = central window of length L, betaObsWhole = whole trial; substream cellKey*1e5 + bi*1e3 + rep"");   % v022", 1};
localFns = J([ ...
    "% --------------------------------------------------------------------------"
    "function [bObs, vObs, nErr, msg] = betaObsOn_local(x, y, FS, EDGE_CLIP, derivConfigs, regressConfigs, nDeriv, nRegress)"
    "% v022: the runner's beta_obs block (processOneTrial_local, top) on a given stretch of trial; ppCM = (rIdx-1)*nDeriv + dIdx."
    "    betaObs = NaN(nDeriv, nRegress);  vgfObs = NaN(nDeriv, nRegress);  nErr = 0;  msg = """";"
    "    for dIdx = 1:nDeriv"
    "        [dx, dy] = differentiateKinematicsEBR(x, y, derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);"
    "        vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2);  vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);"
    "        ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3);  ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);"
    "        sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  lm = [1, -1/3];"
    "        for rIdx = 1:nRegress"
    "            try"
    "                [b, v] = regressDataEBR(sp, kp, regressConfigs(rIdx).type, lm, 0, 0);"
    "                betaObs(dIdx, rIdx) = b;  vgfObs(dIdx, rIdx) = v;"
    "                if rIdx == 1, lm = [v, b]; end"
    "            catch ME"
    "                nErr = nErr + 1;  msg = string(ME.identifier) + "": "" + string(ME.message);"
    "            end"
    "        end"
    "    end"
    "    bObs = betaObs(:)';  vObs = vgfObs(:)';"
    "end"
    ""
    "% --------------------------------------------------------------------------"
    "function r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, noiseModel)"
    "% v022: for endMatched, the record of the end-matched noise and a second inversion from the whole-trial beta_obs against the"
    "% same forward map. Every endMatched result gets the same fields (NaN where the trial stopped early), so the struct array joins."
    "    if ~strcmp(noiseModel, 'endMatched'), return, end"
    "    nanP = NaN(nPipelines, 1);"
    "    r.betaObsBasis = ""window"";"
    "    if ~isfield(res, 'emStatus'), r.betaObs(:) = NaN;  r.vgfObs(:) = NaN; end   % never reached the end-matched step: no window beta_obs"
    "    fl = {'emL', NaN; 'emLfrac', NaN; 'emSeg', NaN(2, 2); 'emAlphaSeg', [NaN NaN]; 'emWindow', [NaN NaN]; 'emStatus', ""not reached (residual below the runner minimum or alpha failure)""; 'nNoiseRejected', NaN; 'betaObsWhole_vec', res.betaObs_vec; 'vgfObsWhole_vec', res.vgfObs_vec};"
    "    for k = 1:size(fl, 1)"
    "        if isfield(res, fl{k, 1}), v = res.(fl{k, 1}); else, v = fl{k, 2}; end"
    "        r.(strrep(fl{k, 1}, '_vec', '')) = v;"
    "    end"
    "    r.betaGenStarWhole = nanP;  r.betaGenStarMedWhole = NaN;  r.ciLoWhole = nanP;  r.ciHiWhole = nanP;  r.invertStatusWhole = repmat(""no_beta_obs"", 1, nPipelines);"
    "    if isfield(res, 'betaObsWhole_vec') && isfield(res, 'betaRec') && any(isfinite(res.betaRec(:)))"
    "        resW = res;  resW.betaObs_vec = res.betaObsWhole_vec;"
    "        [bW, loW, hiW, bWmed, ~, infoW] = invertBeta_local(resW, betaGenVec, nPipelines, monoParams);"
    "        r.betaGenStarWhole = bW;  r.betaGenStarMedWhole = bWmed;  r.ciLoWhole = loW;  r.ciHiWhole = hiW;  r.invertStatusWhole = [infoW.status];"
    "    end"
    "end"
    ""
    "% --------------------------------------------------------------------------"
    "function f0 = estimateF0_local(x, y, fs)"]);
E{end+1} = {"local functions", "% --------------------------------------------------------------------------" + nl + "function f0 = estimateF0_local(x, y, fs)", localFns, 1};
for k = 1:numel(E)
    n = count(s, E{k}{2});
    if n ~= E{k}{4}, error("makeV022:anchor", "%s: expected %d match(es), found %d. v021 has drifted; nothing written.", E{k}{1}, E{k}{4}, n); end
    s = replace(s, E{k}{2}, E{k}{3});
end
nName = count(s, "runLoopClosureFftnoise_v021");  nMat = count(s, "_v021.mat");
if nName < 5 || nMat < 2, error("makeV022:names", "%s", sprintf("unexpected name counts (%d, %d)", nName, nMat)); end
s = replace(s, "runLoopClosureFftnoise_v021", "runLoopClosureFftnoise_v022");  s = replace(s, "_v021.mat", "_v022.mat");
s = replace(s, "fprintf('[v021] Output exists", "fprintf('[v022] Output exists");
s = replace(s, "function results = runLoopClosureFftnoise_v022(datasetName, varargin)" + nl, ...
    "function results = runLoopClosureFftnoise_v022(datasetName, varargin)" + nl + ...
    "% runLoopClosureFftnoise_v022  v021 + NoiseModel 'endMatched' (fairShare common-length segments, shortened synthetic, beta_obs on a same-length central window with the whole trial as sensitivity). GENERATED: see makeRunLoopClosureV022_v001.m. Edit the generator, not this file." + nl);
fid = fopen(outPath, "w", "n", "UTF-8");  if fid < 0, error("makeV022:write", "Cannot write %s", outPath); end
fwrite(fid, unicode2native(char(s), "UTF-8"));  fclose(fid);
fprintf("Wrote %s (%d names, %d .mat literals renamed)\n", outPath, nName, nMat);
end
