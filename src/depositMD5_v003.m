%% depositMD5_v003.m
% MD5 manifest for the Zenodo deposit: explicit file list, one authoritative host per file.
% Zenodo displays and verifies MD5 only. Replaces v002's wildcards (which pulled in Pilot,
% Dagenais, Peters, PDM and Combined files) and adds a HOME column: the host whose copy is
% deposited. Needed because eight files differ between RDS and the Dropbox tree (Session 127
% comparison): the local v007 corpora carry an extra field, vgfObsMM, so local is deposited.
% Run from src/ on each host; on BlueBEAR use sinteractive (never matlab -nodisplay):
%   depositMD5_v003
% On RDS it hashes the home="rds" rows; on any other host (M5, iMac) the home="local" rows.
% Appends one CSV row per file as each hash completes, so a killed session loses nothing.
% A rerun reuses rows already recorded for this host (errors if a file's size has changed).
% Missing home files are collected and reported as an ERROR after everything else is hashed.
% Writes only results/depositMD5_v003.csv (host,tier,home,path,bytes,md5,source).
% Full-fat MD5s below came from one earlier run (depositRDSAudit_v001) and were never
% independently repeated: REVERIFY_FULL = true recomputes them on RDS and errors on any
% mismatch (slow, once: later runs reuse the recorded rows).
% Fraser, D.S. (2026)  v003

%% CONFIG
REVERIFY_FULL = true;
OUT_NAME      = "depositMD5_v003.csv";

DS    = ["Zarandi" "Cook_CTRL" "Cook_ASD" "Dhieb" "Hickman_PLAC" "Hickman_HALO" "Fraser"];
DS13  = ["Zarandi" "Cook_CTRL" "Hickman_PLAC" "Hickman_HALO" "Fraser"];
DS14  = ["Cook_ASD" "Dhieb"];
DS007 = ["Zarandi" "Cook_CTRL" "Cook_ASD" "Dhieb" "Hickman_PLAC" "Hickman_HALO"];
CONST = ["Fraser" "Zarandi" "Cook" "CookASD" "Dhieb" "HickmanPLAC" "HickmanHALO"];
NOISE = ["fraser" "cook" "cookASD" "dhieb" "hickmanPLAC" "hickmanHALO" "zarandi"];
LC    = "src/loopClosureResults_%s_all_shaped_xu_%s.mat";

% Full-fat tier (home RDS): path, bytes, MD5
FULL_PATH  = ["results/powerlaw_debug_v058.db"
              "src/stage1_results_latest.mat"
              "src/stage2_adequacy_latest.mat"
              "src/vgfLMM_v001_latest.mat"
              "src/stage2_vgf_adequacy_latest.mat"
              "src/extractL9_checkpoint_v004.mat"];
FULL_BYTES = [6699646976; 15904865414; 15974589601; 15932561007; 16057450175; 14248827761];
FULL_MD5   = ["a0c0661d4c2425c6c112edc9e565a149"
              "dd5d0a4a7a2ac840c9dacbd9a0bf3d5b"
              "0dbd89a0114fb7c68fbd0eef5635dc7d"
              "3535c2352f0252b89dc3aa69eedcc21b"
              "f367ff3379d5a6881efae22b1dab1a70"
              "90c27945c05015933e93d73a9c9056fb"];

% Minimal tier, home RDS (byte-identical on both trees, or RDS-only)
RDS_MIN = [ ...
    "src/stage3_conditional_L9_20260510_221036.mat"
    "src/stage4_integration_L9_20260523_133842.mat"
    "src/stage3_vgf_nditional_L9_20260523_133822.mat"
    compose(LC, DS', "v012")
    compose(LC, DS13', "v013")
    compose(LC, DS14', "v014")
    compose(LC, DS', "v015")
    "src/loopClosureResults_Fraser_all_shaped_xu_v008.mat"
    "src/monotonicSegments_v2_002.mat"
    "src/perCoordinateSEM_v001.mat"
    "src/perCoordinateSEM_v2_001.mat"
    "src/perCoordinateSEM_VGF_v001.mat"
    "src/constellationMetrics_v003.mat"
    "src/constellationRSA_HPC_v007.mat"
    "src/constellationRSA_HPC_v008.mat"
    compose("src/noiseCharacterisation_%s.mat", NOISE')
    "src/measureRegressionSaturation_v005.mat"
    "src/checkD10VGFCoeffCI_v002_rawCoeffs.mat"
    "src/checkD10VGFCoeffCI_v002_results.mat"
    "results/L9ModelSlim_v001.mat"
    "results/tempoFoldSweep_v001.mat"
    "results/tempoFoldSweepC_v001.mat"
    "results/tempoFoldSweepE_v001.mat"
    "results/gridTempoAtFold_v004.mat"
    "results/simpleEffects_v004_L9_summary.txt"];

% Minimal tier, home local (differs from RDS, or exists only in the Dropbox tree)
LOC_MIN = [ ...
    compose(LC, DS007', "v007")
    "src/constellationMetrics_v002.mat"
    "src/constellationMetrics_v004.mat"
    "src/predictionLookup_v2_002.mat"
    compose("src/constellation%s_v001.mat", CONST')
    "src/loopClosureLOPOOutOfSample_v001.mat"
    "src/loopClosureVarDecomp_v011_all6_N20.mat"
    "src/loopClosureVarDecomp_v011_all6_N200.mat"
    "src/loopClosureVarDecomp_v012_gated_all6_legacy.mat"
    "src/loopClosureVarDecomp_v012_gated_all6_gatedV012.mat"
    "src/loopClosureVarDecomp_v012_gated_all6_confirmatory.mat"
    "src/loopClosureVarDecomp_v012_gated_all6_v015.mat"
    "src/simpleEffects_L9_20260912_101823.mat"
    "results/L9ModelInspect_v002.mat"
    "results/simpleEffects_v004_L9_20260928_151609.mat"
    "results/simpleEffects_v004_log_20260928_143752.txt"];

%% Manifest table and its own consistency checks
M = [blk("full", "rds", FULL_PATH); blk("minimal", "rds", RDS_MIN); blk("minimal", "local", LOC_MIN)];
if numel(unique(M.path)) ~= height(M)
    [~, ia] = unique(M.path); dup = M.path(setdiff(1:height(M), ia));
    error("depositMD5:dup", "%s", "Duplicate manifest paths: " + strjoin(dup', ", "));
end

%% Root, host, tool (root from this script's own location; run from src/)
root = string(fileparts(fileparts(mfilename("fullpath"))));
if ~isfolder(fullfile(root, "src")) || ~isfolder(fullfile(root, "results"))
    error("depositMD5:root", "%s", "Expected src/ and results/ under " + root + ". Run from src/.");
end
if ~isunix
    error("depositMD5:os", "%s", "Needs md5sum (Linux) or md5 (macOS).");
end
if ismac, hashCmd = "md5 -q"; else, hashCmd = "md5sum"; end
[~, hostRaw] = system("hostname");
host = strtrim(string(hostRaw));
thisHome = "local"; if startsWith(root, "/rds/"), thisHome = "rds"; end
outCsv = fullfile(root, "results", OUT_NAME);
fprintf("Root: %s\nHost: %s (%s), hashing home=%s rows\nOutput: %s\n\n", root, host, hashCmd, thisHome, outCsv);

%% Resume: rows already recorded for this host
prev = containers.Map('KeyType', 'char', 'ValueType', 'any');
if isfile(outCsv)
    L = readlines(outCsv); L = L(strlength(L) > 0);
    for k = 2:numel(L)
        c = split(L(k), ",");
        if c(1) == host, prev(char(c(4))) = struct('bytes', str2double(c(5)), 'md5', c(6)); end
    end
else
    writelines("host,tier,home,path,bytes,md5,source", outCsv);
end

%% Hash the rows whose home is this host
mine = M(M.home == thisHome, :);
missing = strings(0, 1); nFile = zeros(1, 2); nByte = zeros(1, 2);
for i = 1:height(mine)
    p = mine.path(i); f = fullfile(root, p); isFull = mine.tier(i) == "full";
    if ~isfile(f), missing(end+1, 1) = p; continue, end %#ok<SAGROW>
    d = dir(f); b = d.bytes; key = char(p);
    if isFull
        j = find(FULL_PATH == p, 1);
        if b ~= FULL_BYTES(j)
            error("depositMD5:size", "%s", p + " is " + b + " bytes, known " + FULL_BYTES(j) + ": changed since hashing.");
        end
    end
    if isKey(prev, key)
        e = prev(key);
        if e.bytes ~= b
            error("depositMD5:changed", "%s", p + ": " + b + " bytes now, " + e.bytes + " when recorded.");
        end
        h = e.md5; how = "recorded";
    else
        if isFull && ~REVERIFY_FULL
            h = FULL_MD5(j); how = "carried-over";
        else
            h = md5Of(f, hashCmd); how = "computed";
            if isFull && h ~= FULL_MD5(j)
                error("depositMD5:md5", "%s", p + ": md5 " + h + " differs from earlier " + FULL_MD5(j));
            end
        end
        writelines(host + "," + mine.tier(i) + "," + mine.home(i) + "," + p + "," + b + "," + h + "," + how, ...
            outCsv, WriteMode="append");
    end
    t = 1 + ~isFull; nFile(t) = nFile(t) + 1; nByte(t) = nByte(t) + b;
    fprintf("  %-7s %-58s %14d  %s  (%s)\n", mine.tier(i), p, b, h, how);
end

%% Summary, then fail loud on anything missing
fprintf("\nThis host (%s): full-fat %d files %.2f GB; minimal %d files %.3f GB\n", ...
    thisHome, nFile(1), nByte(1) / 1e9, nFile(2), nByte(2) / 1e9);
fprintf("Rows homed elsewhere (not hashed here): %d\nCSV: %s\n", sum(M.home ~= thisHome), outCsv);
if ~isempty(missing)
    error("depositMD5:missing", "%s", numel(missing) + " manifest file(s) missing on their home host " + thisHome + ": " + strjoin(missing', ", "));
end

%% Local functions
function t = blk(tier, home, paths)
% One manifest block: a table of (tier, home, path) rows.
paths = paths(:); n = numel(paths);
t = table(repmat(string(tier), n, 1), repmat(string(home), n, 1), paths, 'VariableNames', ["tier" "home" "path"]);
end

function h = md5Of(f, hashCmd)
% Runs the host's md5 tool; errors loudly on a non-zero status or unparsable output.
[st, out] = system(hashCmd + " """ + f + """");
if st ~= 0
    error("depositMD5:hash", "%s", "md5 failed for " + f + ": " + strtrim(out));
end
tk = regexp(out, "^([0-9a-f]{32})", "tokens", "once");
if isempty(tk)
    error("depositMD5:parse", "%s", "Could not parse md5 output for " + f + ": " + strtrim(out));
end
h = string(tk{1});
end
