%% checkTempoBinSpec_v001  Deep Review GD17 item 7b: the f0 bin edges and the per-bin SG-IRLS target medians (including Cook CTRL) behind the "sliding target" paragraph.
% Parses extractGainTarget_v002.m for EDGES / BIN_NAMES (a missing pattern is an error) and reads B3 from results/gainTarget_v002.mat
% (dataset x bin: median gain, median target, n reliable targets). Saves the in-domain table for the ledger.
% Run from src/. Writes results/checkTempoBinSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
sc = fileread(fullfile(ROOT, "src", "extractGainTarget_v002.m"));
assert(~isempty(regexp(sc, "EDGES\s+= \[0 0\.35 0\.7 1\.4 2\.8 Inf\];\s+BIN_NAMES = \[""<0\.35"" ""0\.35-0\.7"" ""0\.7-1\.4"" ""1\.4-2\.8"" "">2\.8""\];", "once")), "binSpec:parse", "bin edges not found");
f = fullfile(ROOT, "results", "gainTarget_v002.mat"); if ~isfile(f), error("binSpec:file", "%s", "FAILED PATH: " + f); end
B = load(f, "B3").B3;  B.Properties.VariableNames = ["dataset" "bin" "nTrials" "medGain" "medGainCount" "medTarget" "nReliable"];
T = B(:, ["dataset" "bin" "medTarget" "nReliable"]);
q = @(d, b) T.medTarget(T.dataset == d & T.bin == b);  n = @(d, b) T.nReliable(T.dataset == d & T.bin == b);
bins4 = ["<0.35" "0.35-0.7" "0.7-1.4" "1.4-2.8"];
assert(isempty(n("Cook_CTRL", "1.4-2.8")) && isempty(n("Cook_CTRL", ">2.8")), "binSpec:ctrl", "Cook CTRL unexpectedly has trials in a bin above 1.4 Hz");
CTRL = arrayfun(@(b) q("Cook_CTRL", b), bins4(1:3));  CTRLn = arrayfun(@(b) n("Cook_CTRL", b), bins4(1:3));
HI = struct("hickmanPLAC", [q("Hickman_PLAC", ">2.8") n("Hickman_PLAC", ">2.8")], "hickmanHALO", [q("Hickman_HALO", ">2.8") n("Hickman_HALO", ">2.8")], "dhieb", [q("Dhieb", ">2.8") n("Dhieb", ">2.8")]);
S = struct("T", T, "CTRL", CTRL, "CTRLn", CTRLn, "HI", HI, "edges", [0 0.35 0.7 1.4 2.8 Inf]);
disp(S.CTRL); disp(S.CTRLn); disp(HI);  save(fullfile(ROOT, "results", "checkTempoBinSpec_v001.mat"), "-struct", "S");
