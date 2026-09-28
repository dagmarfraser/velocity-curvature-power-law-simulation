%% extractL9ModelSlim_v001.m
% Save the small pieces of the fitted L9 Stage 1 LMM needed by inspectL9Model_v002, from
% an lme ALREADY IN THE BASE WORKSPACE (left there by interrupting inspectL9Model_v001
% after its ~10 h single-core load of src/extractL9_checkpoint_v004.mat). It never
% reloads: if lme is absent it stops, because reloading costs hours.
% Each step is timed and appended to the output immediately, so a stall or interruption
% keeps everything already extracted. If a step stalls, Ctrl+C: the stack names the
% property, and lme stays in memory.
% Run on BlueBEAR from src/, in the SAME session:  cd src; extractL9ModelSlim_v001
% Paths resolve from this file's location (path fix 2026-09-28, in place; the first run
% was from the repo root and wrote the same file).
% Writes results/L9ModelSlim_v001.mat (cn, fx, se, tn, Tm, vn, V; a few kB).
% Fraser, D.S. (2026)  v001

%% CONFIG
SRC_DIR  = fileparts(mfilename("fullpath"));
OUT_SLIM = fullfile(fileparts(SRC_DIR), "results", "L9ModelSlim_v001.mat");

%% Guard: lme must already be in memory
if ~exist("lme", "var")
    error("slimL9:noLme", "%s", "No variable 'lme' in this workspace. Run in the session where " + ...
        "inspectL9Model_v001 loaded it (Ctrl+C that script first). Not reloading: the load takes hours.");
end
if ~isa(lme, "LinearMixedModel")
    error("slimL9:class", "%s", "'lme' is a " + class(lme) + ", not a LinearMixedModel");
end
if ~isfolder(fileparts(OUT_SLIM))
    error("slimL9:outDir", "%s", "Missing folder: " + fileparts(OUT_SLIM));
end
fprintf("lme found: %s, %d observations. Writing %s\n", class(lme), lme.NumObservations, OUT_SLIM);
tSlim = tic;
source = "src/extractL9_checkpoint_v004.mat via inspectL9Model_v001 (BlueBEAR)";
runDate = string(datetime("now"));
save(OUT_SLIM, "source", "runDate");                         % create the file first

%% Steps: each timed and appended
fprintf("[1/5] coefficient names and fixed effects ... ");
cn = string(lme.CoefficientNames(:));  fx = fixedEffects(lme);
save(OUT_SLIM, "cn", "fx", "-append");  fprintf("%d coefficients, %.1f s\n", numel(cn), toc(tSlim));

fprintf("[2/5] coefficient SEs ... ");
se = sqrt(diag(lme.CoefficientCovariance));
save(OUT_SLIM, "se", "-append");  fprintf("%.1f s\n", toc(tSlim));

fprintf("[3/5] fixed-effects formula ... ");
fe = lme.Formula.FELinearFormula;  fprintf("%.1f s\n", toc(tSlim));

fprintf("[4/5] term names, term matrix, variable names ... ");
tn = string(fe.TermNames(:));  Tm = fe.Terms;  vn = string(fe.VariableNames(:));
save(OUT_SLIM, "tn", "Tm", "vn", "-append");  fprintf("%d terms, %.1f s\n", numel(tn), toc(tSlim));

fprintf("[5/5] variable info (formula variables only) ... ");
V = lme.VariableInfo;  V = V(cellstr(vn), :);
save(OUT_SLIM, "V", "-append");  fprintf("%.1f s\n", toc(tSlim));

d = dir(OUT_SLIM);
fprintf("Done in %.1f s. %s (%.1f kB). Copy to the iMac's results/ and run inspectL9Model_v002.\n", ...
    toc(tSlim), OUT_SLIM, d.bytes / 1024);
