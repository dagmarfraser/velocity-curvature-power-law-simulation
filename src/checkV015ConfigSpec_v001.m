%% checkV015ConfigSpec_v001  Deep Review GD18: the four settings the v015 result files do not record (BW_CUTOFF, SG_REF, EDGE_CLIP, FS_DS) and where v015 fixes them.
% Loads only `config` from each loopClosureResults_<dataset>_all_shaped_xu_v015.mat (fails loud if a file is missing, if a field IS present, or if the file's
% own settings differ from N_REPS 200 / N_BETA 25 / seed 1729) and parses runLoopClosureFftnoise_v015.m for the code-fixed values. Saves the dataset-to-FS_DS map.
% Run from src/. Writes results/checkV015ConfigSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
DS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"];  FLD = ["BW_CUTOFF" "SG_REF" "EDGE_CLIP" "FS_DS"];
for d = DS
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");  if ~isfile(f), error("cfg15:file", "%s", "FAILED PATH: " + f); end
    c = load(f, "config").config;
    if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error("cfg15:cfg", "%s", d + ": N_REPS/N_BETA/seed differ from 200/25/1729"); end
    if any(isfield(c, FLD)), error("cfg15:present", "%s", d + ": a field is present; the 'not recorded' statement needs revisiting"); end
end
r = fileread(fullfile(ROOT, "src", "runLoopClosureFftnoise_v015.m"));
need = @(pat, what) assert(~isempty(regexp(r, pat, "once")), "cfg15:parse", "Pattern not found: %s", what);
need("'filterParams',\s*\{\[2 10 1\],\[4 17\]\}", "filterParams {[2 10 1],[4 17]}");
need("addParameter\(p, 'EdgeClip',\s+20,", "EdgeClip default 20");
tok = regexp(r, "case '([^']+)'(?:(?!case ')[\s\S])*?FS_DS\s*=\s*(\d+);", "tokens");
M = containers.Map();  for i = 1:numel(tok), M(strrep(tok{i}{1}, " ", "_")) = str2double(tok{i}{2}); end
S = struct("fsDS", M, "bwCutoff", 10, "sgRef", 17, "edgeClip", 20);   %#ok<NASGU>
disp(keys(M)); disp(cell2mat(values(M)));
out = struct("datasets", string(keys(M)), "fsDS", cell2mat(values(M)), "bwCutoff", 10, "sgRef", 17, "edgeClip", 20);
save(fullfile(ROOT, "results", "checkV015ConfigSpec_v001.mat"), "-struct", "out");
