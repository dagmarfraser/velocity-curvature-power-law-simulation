function [results, meta] = loadV028_v002(dataset, opts)
%LOADV028_V002 v028 per-trial loader with a selectable number of harmonics K (nH arm). New file, 2026-11-21; loadV028_v001 is untouched.
%   [results, meta] = loadV028_v002("Cook ASD")                          K=4 primary arm: returns exactly what loadV028_v001 returns
%   loadV028_v002("Cook ASD", NHarm=1)                                   K=1 (primary arm, window basis) from BRM/scratch_K, md5 checked
%   loadV028_v002("Cook ASD", NHarm=0)                                   K=0 from BRM/outputs_K0 (v033 files; no recorded md5 list there)
%   NHarm in {0,1,2,3,4,6,10}. NHarm=4 (default) is the old arms exactly as in v001 (Arm primary/budget1/shapeXu/harmonicsFc, Seed, ObsBasis).
%   NHarm ~= 4 exists only for Arm="primary", Seed=1729; any other combination is an error, not a fallback.
%   Directories: K=4 looks in SrcDir (old-tree src) first, then Dir, exactly as v001. NHarm ~= 4 looks ONLY in Dir (default BRM/scratch_K for
%   K=1,2,3,6,10; BRM/outputs_K0 for K=0). A missing file is an error. If Dir holds md5_local_K1plus.txt, the file's md5 must be listed there and
%   must match, else error; where no list exists (K=0 folder, K=4 src) no check is made and meta.md5Checked is false.
%   Everything after the file is chosen (ObsBasis mapping, meta) is the v001 code, unchanged.
arguments
    dataset (1,1) string
    opts.Arm (1,1) string {mustBeMember(opts.Arm, ["primary", "budget1", "shapeXu", "harmonicsFc"])} = "primary"
    opts.NHarm (1,1) double {mustBeMember(opts.NHarm, [0 1 2 3 4 6 10])} = 4
    opts.ObsBasis (1,1) string {mustBeMember(opts.ObsBasis, ["window", "whole"])} = "window"
    opts.Seed (1,1) double = 1729
    opts.Dir (1,1) string = ""
    opts.SrcDir (1,1) string = ""
end
home = string(getenv('HOME'));
b2b = fullfile(home, 'Library', 'CloudStorage', 'Dropbox', 'Brain2Bee');
srcDir = opts.SrcDir;  if srcDir == "", srcDir = string(fullfile(b2b, 'PowerLawSimulationPreReg', 'src')); end
if ~isfolder(srcDir), error('loadV028:srcdir', 'SrcDir %s does not exist', srcDir); end
tag = strrep(dataset, " ", "_");  if tag == "Fraser", tag = "FraserN48"; end
switch opts.Arm
    case "primary",     arm = "endMatched_segdetect";
    case "budget1",     arm = "endMatched_bud1";
    case "shapeXu",     arm = "shapedXuSeg_segdetect";
    case "harmonicsFc", arm = "endMatched_segdetect_nHfc";
end
extra = "";  if opts.Seed ~= 1729, extra = "_seed" + opts.Seed + "_r50"; end
ver = "v028";  md5Checked = false;  md5Value = "";
if opts.NHarm == 4
    dirs = [srcDir, opts.Dir];
else
    if opts.Arm ~= "primary" || opts.Seed ~= 1729
        error('loadV028:nharmArm', 'NHarm=%d exists only for Arm="primary" and Seed=1729 (got Arm=%s, Seed=%g)', opts.NHarm, opts.Arm, opts.Seed);
    end
    arm = arm + "_nH" + opts.NHarm;
    D = opts.Dir;
    if D == ""
        if opts.NHarm == 0, D = string(fullfile(b2b, 'PowerLaw_Repro_2026', 'BRM', 'outputs_K0'));
        else,               D = string(fullfile(b2b, 'PowerLaw_Repro_2026', 'BRM', 'scratch_K')); end
    end
    if opts.NHarm == 0, ver = "v033"; end
    dirs = D;
end
name = "loopClosureResults_" + tag + "_all_" + arm + extra + "_" + ver + ".mat";
f = "";
for d = dirs(strlength(dirs) > 0), if isfile(fullfile(d, name)), f = fullfile(d, name);  break, end, end
if f == "", error('loadV028:missing', 'v028 file %s not found in %s', name, strjoin(dirs, ' or ')); end
if opts.NHarm ~= 4
    lst = fullfile(fileparts(f), 'md5_local_K1plus.txt');
    if isfile(lst)
        L = splitlines(string(fileread(lst)));  L = L(strlength(strtrim(L)) > 0);
        parts = split(strtrim(L));  if isvector(parts), parts = parts(:)'; end
        j = find(parts(:, end) == name);
        if isempty(j), error('loadV028:md5unlisted', '%s is not in %s', name, lst); end
        md5Value = md5File_local(f);
        if lower(parts(j(1), 1)) ~= md5Value, error('loadV028:md5', 'md5 of %s is %s, list says %s', f, md5Value, parts(j(1), 1)); end
        md5Checked = true;
    end
end
S = load(f);  results = S.results;
if opts.ObsBasis == "whole"
    map = ["betaObs" "betaObsWhole"; "betaGenStar" "betaGenStarWhole"; "betaGenStarMed" "betaGenStarMedWhole"; "ciLo" "ciLoWhole"; ...
           "ciHi" "ciHiWhole"; "invertStatus" "invertStatusWhole"; "vgfObs" "vgfObsWhole"];
    for k = 1:size(map, 1)
        assert(isfield(results, map(k, 2)), 'loadV028:field', '%s lacks %s', name, map(k, 2));
        [results.(map(k, 1))] = results.(map(k, 2));
    end
end
meta = struct('file', f, 'config', S.config, 'betaGenVec', S.betaGenVec, 'pipelineLabels', {S.pipelineLabels}, 'loopCCC', S.loopCCC, ...
    'arm', opts.Arm, 'obsBasis', opts.ObsBasis, 'dataset', tag);
meta.nHarm = opts.NHarm;  meta.md5Checked = md5Checked;  meta.md5 = md5Value;   % v002 additions; the nine fields above are v001's
end

function h = md5File_local(f)
fid = fopen(f, 'r');  if fid < 0, error('loadV028:open', 'cannot open %s', f); end
c = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('MD5');
while true
    b = fread(fid, 16*2^20, '*uint8');
    if isempty(b), break, end
    md.update(b);
end
h = lower(string(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
