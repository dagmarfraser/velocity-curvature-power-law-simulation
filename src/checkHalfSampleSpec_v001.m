%% checkHalfSampleSpec_v001  Deep Review GD18: numbers in Part 1B's half-sample-offset sentence, from predictHarmonicGain_v002.mat (T_bw).
% BW-OLS gain-0.9 onsets in f0/f_c at f_c 5, 10, 20 Hz: measured pipeline, emulation with the BWFD half-sample delay (emul), emulation with it removed (emulZP), theory.
% Spread = max - min over the three cutoffs. Delay effect at each cutoff = emulZPLoc/emulLoc - 1 (local gain window). Fails loud if rows are missing.
% Run from src/. Writes results/checkHalfSampleSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
f = fullfile(ROOT, "results", "predictHarmonicGain_v002.mat");  if ~isfile(f), error("hs:file", "%s", "FAILED PATH: " + f); end
T = load(f, "T_bw").T_bw;  T = T(string(T.regression) == "OLS", :);  T = sortrows(T, "fc");
if ~isequal(T.fc(:)', [5 10 20]), error("hs:rows", "expected BW-OLS rows at fc 5, 10, 20"); end
sp = @(x) max(x) - min(x);
out = struct("fc", T.fc', "spreadMeasured", sp(T.measured), "spreadEmul", sp(T.emul), "spreadEmulZP", sp(T.emulZP), "spreadTheory", sp(T.theory), ...
    "delayEffectPct", 100*(T.emulZPLoc ./ T.emulLoc - 1)');
fprintf("spread measured %.4f emul %.4f emulZP %.4f theory %.4f; delay effect %% at fc 5/10/20: %s\n", out.spreadMeasured, out.spreadEmul, out.spreadEmulZP, out.spreadTheory, mat2str(round(out.delayEffectPct,2)));
save(fullfile(ROOT, "results", "checkHalfSampleSpec_v001.mat"), "-struct", "out");
