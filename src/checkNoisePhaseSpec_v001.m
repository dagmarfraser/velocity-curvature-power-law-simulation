%% checkNoisePhaseSpec_v001  Deep Review GD18: fftnoise and shaped_xu share the trial's Fourier amplitude spectrum and differ in phase.
% Parses src/functions/generateLoopClosureNoise_v003.m (the production generator named in the v015 runner) and stops if a pattern is missing:
%   shaped_xu: Eemp = abs(fft(resid)); Fout = Eemp .* exp(1i*angle(Fb)); then a single scale factor restores the injection sigma.
%   fftnoise : the case calls fftnoise() on the FFT of the residual (phase-scramble, amplitude preserved).
% Run from src/. Writes results/checkNoisePhaseSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
f = fullfile(ROOT, "src", "functions", "generateLoopClosureNoise_v003.m");  if ~isfile(f), error("nph:file", "%s", "FAILED PATH: " + f); end
r = fileread(f);
need = @(pat, what) assert(~isempty(regexp(r, pat, "once")), "nph:parse", "Pattern not found: %s", what);
need("Eemp\s*=\s*abs\(fft\(resid\)\)", "empirical magnitude Eemp = abs(fft(resid))");
need("Fout\s*=\s*Eemp\s*\.\*\s*exp\(1i\*angle\(Fb\)\)", "Fout = Eemp .* exp(1i*angle(Fb))");
need("Fb\s*=\s*fft\(base\)", "Fb = fft(base)");
need("n\s*=\s*n\s*\*\s*\(sig\s*/\s*sN\)", "single scale factor sig/sN");
need("nMaj\s*=\s*fftnoise\(srcMaj\)", "fftnoise case");
out = struct("shapedXuFormula", "|FFT(resid)| .* exp(1i*angle(FFT(Xu_base)))", "scaleFactorOnly", true, "generator", "generateLoopClosureNoise_v003.m"); %#ok<NASGU>
save(fullfile(ROOT, "results", "checkNoisePhaseSpec_v001.mat"), "-struct", "out");
fprintf("checkNoisePhaseSpec_v001: patterns found in %s\n", f);
