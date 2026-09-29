%% depositRDSAudit_v001.m
% Zenodo deposit audit on BlueBEAR (RDS): existence, bytes, variable inventory of the
% two big Stage 1 sources, MD5/SHA-256 of every full-fat file, size of the minimal tier,
% AppleDouble (._*) litter, and the DB-name question.
% Run from src/ in sinteractive (never matlab -nodisplay):  depositRDSAudit_v001
% Writes results/depositRDSAudit_v001.csv (one row per file) and .txt (console log).
% Only write besides those: copies results/LMM_coefficients_top_L9_20260525_093556.txt
% into src/ if absent (no overwrite). Nothing is deleted or moved.
% Fraser, D.S. (2026)  v001

%% CONFIG
HASH_MD5    = true;      % Zenodo displays MD5, so this is the post-upload check
HASH_SHA256 = true;      % kept in the manifest for the archive
INVENTORY   = true;      % whos -file on the checkpoint and stage1_results_latest
COPY_LMMTOP = true;
LMMTOP      = "LMM_coefficients_top_L9_20260525_093556.txt";

% Full-fat tier: [path relative to repo root, expected bytes (NaN = unknown, report only)]
FULL = { ...
    "results/powerlaw_debug_v058.db",            6699646976; ...
    "src/stage1_results_latest.mat",            15904865414; ...
    "src/stage2_adequacy_latest.mat",           15974589601; ...
    "src/vgfLMM_v001_latest.mat",               15932561007; ...
    "src/stage2_vgf_adequacy_latest.mat",       16057450175; ...
    "src/extractL9_checkpoint_v004.mat",        NaN};

% Minimal tier: exact names first, then globs (relative to repo root)
MINIMAL = [ ...
    "src/stage3_conditional_L9_20260510_221036.mat"
    "src/stage4_integration_L9_20260523_133842.mat"
    "src/stage3_vgf_nditional_L9_20260523_133822.mat"
    "src/loopClosureResults_*_all_shaped_xu_v012.mat"
    "src/loopClosureResults_*_all_shaped_xu_v013.mat"
    "src/loopClosureResults_*_all_shaped_xu_v014.mat"
    "src/loopClosureResults_*_all_shaped_xu_v015.mat"
    "src/loopClosureResults_*_all_shaped_xu_v007.mat"
    "src/loopClosureResults_Fraser_all_shaped_xu_v008.mat"
    "src/monotonicSegments_v2_002.mat"
    "src/perCoordinateSEM_v001.mat"
    "src/perCoordinateSEM_v2_001.mat"
    "src/perCoordinateSEM_VGF_v001.mat"
    "src/predictionLookup_v2_002.mat"
    "src/constellation*_v001.mat"
    "src/constellationMetrics_v00*.mat"
    "src/constellationRSA_HPC_v007.mat"
    "src/constellationRSA_HPC_v008.mat"
    "src/noiseCharacterisation_*.mat"
    "src/measureRegressionSaturation_v005.mat"
    "src/loopClosureLOPOOutOfSample_v001.mat"
    "src/loopClosureVarDecomp_v01*_all6_*.mat"
    "src/checkD10VGFCoeffCI_v002_*.mat"
    "results/L9ModelSlim_v001.mat"
    "results/L9ModelInspect_v002.mat"
    "results/tempoFoldSweep*_v001.mat"
    "results/gridTempoAtFold_v004.mat"
    "results/simpleEffects_v004_*"];

%% Root from this script's own location (run from src/)
here = fileparts(mfilename("fullpath"));
root = fileparts(here);
if ~isfolder(fullfile(root, "src")) || ~isfolder(fullfile(root, "results"))
    error("depositAudit:root", "%s", "Expected src/ and results/ under " + root + ". Run from src/.");
end
if ~isunix
    error("depositAudit:os", "%s", "Needs md5sum/sha256sum; run on BlueBEAR.");
end
logFile = fullfile(root, "results", "depositRDSAudit_v001.txt");
diary(logFile); diary on;
cleanupObj = onCleanup(@() diary("off")); %#ok<NASGU>
fprintf("Root: %s\nHost: %s\n\n", root, getenv("HOSTNAME"));

%% 1. Full-fat tier: existence, bytes, hashes
rows = table('Size', [0 7], 'VariableTypes', ["string" "string" "double" "double" "string" "string" "string"], ...
    'VariableNames', ["tier" "path" "bytes" "expected" "match" "md5" "sha256"]);
for i = 1:size(FULL, 1)
    p = FULL{i, 1}; f = fullfile(root, p);
    if ~isfile(f)
        error("depositAudit:missing", "%s", "Full-fat file missing: " + f);
    end
    d = dir(f); b = d.bytes; e = FULL{i, 2};
    if isnan(e), m = "n/a (no expected size)"; elseif b == e, m = "OK"; else, m = "MISMATCH"; end
    if m == "MISMATCH"
        warning("depositAudit:size", "%s", p + ": " + b + " bytes, expected " + e);
    end
    h5 = ""; h256 = "";
    if HASH_MD5,    h5   = hashFile(f, "md5sum"); end
    if HASH_SHA256, h256 = hashFile(f, "sha256sum"); end
    fprintf("%-46s %16d  %s\n   md5 %s\n   sha256 %s\n", p, b, m, h5, h256);
    rows(end+1, :) = {"full", string(p), b, e, m, h5, h256}; %#ok<SAGROW>
end

%% 2. Minimal tier: resolve, size, hash
fprintf("\nMinimal tier\n");
minN = 0; minB = 0;
for i = 1:numel(MINIMAL)
    g = fullfile(root, MINIMAL(i));
    d = dir(g); d = d(~[d.isdir] & ~startsWith({d.name}, "._"));
    if isempty(d)
        warning("depositAudit:noMatch", "%s", "No match: " + MINIMAL(i));
        continue
    end
    for k = 1:numel(d)
        f = fullfile(d(k).folder, d(k).name);
        h5 = ""; h256 = "";
        if HASH_MD5,    h5   = hashFile(f, "md5sum"); end
        if HASH_SHA256, h256 = hashFile(f, "sha256sum"); end
        rel = string(strrep(f, root + filesep, ""));
        rows(end+1, :) = {"minimal", rel, d(k).bytes, NaN, "n/a", h5, h256}; %#ok<SAGROW>
        minN = minN + 1; minB = minB + d(k).bytes;
    end
    fprintf("  %-58s %3d file(s)\n", MINIMAL(i), numel(d));
end
fprintf("Minimal tier: %d files, %.3f GB\n", minN, minB / 1e9);
fullB = sum(rows.bytes(rows.tier == "full"));
fprintf("Full-fat tier: %d files, %.2f GB (%d bytes)\nBoth tiers: %.2f GB\n", ...
    sum(rows.tier == "full"), fullB / 1e9, fullB, (fullB + minB) / 1e9);

%% 3. Variable inventory of the two big Stage 1 sources (header read, no data load)
if INVENTORY
    fprintf("\nVariable inventory (matfile, no load)\n");
    for nm = ["extractL9_checkpoint_v004.mat", "stage1_results_latest.mat"]
        f = fullfile(root, "src", nm);
        t0 = tic; w = whos(matfile(f));
        fprintf("%s  (%.1f s, %d variables)\n", nm, toc(t0), numel(w));
        for k = 1:numel(w)
            fprintf("   %-32s %-28s %14d bytes  [%s]\n", w(k).name, mat2str(w(k).size), w(k).bytes, w(k).class);
        end
        hasLme = any(strcmp({w.class}, "LinearMixedModel")) || any(strcmp({w.name}, "lme"));
        fprintf("   -> holds a fitted lme object: %s\n", string(hasLme));
    end
end

%% 4. DB names present on RDS results/
fprintf("\nDatabase files in results/\n");
d = dir(fullfile(root, "results", "*.db")); d = d(~startsWith({d.name}, "._"));
for k = 1:numel(d), fprintf("   %-34s %14d bytes\n", d(k).name, d(k).bytes); end
for nm = ["powerlaw_debug_v058.db", "powerlaw_multiverse_v058.db", "powerlaw_multiverse_v057.db"]
    fprintf("   %-34s %s\n", nm, string(isfile(fullfile(root, "results", nm))));
end

%% 5. AppleDouble litter (exclude from any tar/upload)
fprintf("\nAppleDouble (._*) files\n");
for sub = ["src", "results", "data", "figures"]
    d = dir(fullfile(root, sub, "**", "._*"));
    fprintf("   %-8s %4d files, %d bytes\n", sub, numel(d), sum([d.bytes]));
end

%% 6. LMM_coefficients_top open item
src = fullfile(root, "results", LMMTOP); dst = fullfile(root, "src", LMMTOP);
fprintf("\n%s\n   in results/: %s   in src/: %s\n", LMMTOP, string(isfile(src)), string(isfile(dst)));
if COPY_LMMTOP && isfile(src) && ~isfile(dst)
    [ok, msg] = copyfile(src, dst);
    if ~ok, error("depositAudit:copy", "%s", msg); end
    fprintf("   copied to src/ (no overwrite)\n");
end

%% Output
writetable(rows, fullfile(root, "results", "depositRDSAudit_v001.csv"));
fprintf("\nWrote results/depositRDSAudit_v001.csv (%d rows) and .txt\n", height(rows));

%% Local function
function h = hashFile(f, tool)
% Runs md5sum/sha256sum; errors loudly on any non-zero status or unparsable output.
[st, out] = system(sprintf('%s "%s"', tool, f));
if st ~= 0
    error("depositAudit:hash", "%s", tool + " failed for " + f + ": " + strtrim(out));
end
tk = regexp(out, "^([0-9a-f]{32,64})\s", "tokens", "once");
if isempty(tk)
    error("depositAudit:parse", "%s", "Could not parse " + tool + " output for " + f + ": " + strtrim(out));
end
h = string(tk{1});
end
