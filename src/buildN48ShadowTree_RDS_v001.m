function buildN48ShadowTree_RDS_v001()
%BUILDN48SHADOWTREE_RDS_V001  Build the N=48 copy of the project on RDS, from MATLAB (2026-11-04). Open in the GUI on BlueBEAR
%   (sinteractive) and press Run. No arguments. Replaces the .sh version of the same name (no shell is available).
%   WHY: the HPC launchers read the Fraser N=45 corpus and outputs by fixed names. In a copy of the project where those Fraser
%   inputs ARE the N=48 ones, the launchers run UNCHANGED. Nothing in the real tree is modified: everything is copied, and
%   every output of the launchers lands in the copy, <project>_N48, next to the project folder.
%   SWAPS (as buildN48ShadowTree_v001 on the Mac):
%     src/loopClosureResults_Fraser_all_shaped_xu_v015.mat <- loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat
%     src/loopClosureResults_Fraser_all_shaped_xu_v019.mat <- loopClosureResults_FraserN48_all_shaped_xu_v020.mat
%     src/noiseCharacterisation_fraser.mat                 <- noiseCharacterisation_fraserN48.mat
%     src/constellationFraser_v001.mat                     <- constellationFraser_v002.mat
%     src/functions/importDB_fraser_v001.m                 <- shim calling importDB_fraser_v002
%     every src/**/*.m: 10.41793 -> 5.208965 (UIKit points per mm)
%   PLUS (HPC only, Rerun Step 3): every other data file (.mat .txt .csv .log) and *_chunks folder with "fraser" in its name
%   under src/ and results/ is MOVED, inside the copy, to N45_moved/, so the launchers cannot skip Fraser work as already
%   done and no N=45 Fraser intermediate is read silently. Other datasets' outputs stay, so only Fraser is recomputed.
%   Fails loud on any missing input or existing copy. Writes <copy>/N48_SWAPS.txt.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  proj = fileparts(srcDir);  sh = [proj '_N48'];
if isfolder(sh), error('N48rds:exists', '%s', ['Copy exists: ' sh ' (move it away to rebuild).']); end
need = ["loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat" "loopClosureResults_FraserN48_all_shaped_xu_v020.mat" ...
        "noiseCharacterisation_fraserN48.mat" "constellationFraser_v002.mat" "functions/importDB_fraser_v002.m"];
for f = need, if ~isfile(fullfile(srcDir, f)), error('N48rds:missing', '%s', "Missing N=48 input: " + fullfile(srcDir, f)); end, end
t0 = tic;  mkdir(sh);
fprintf('Copying src/ (%s) ...\n', char(datetime('now')));      ok_local(copyfile(srcDir, fullfile(sh, 'src')), 'copy src');
fprintf('Copying results/ (%s) ...\n', char(datetime('now'))); ok_local(copyfile(fullfile(proj, 'results'), fullfile(sh, 'results')), 'copy results');
for d = ["figures" "logs" "docs" "N45_moved"], mkdir(fullfile(sh, d)); end
% raw data (read-only use): the importers resolve ../data from their own location, so the copy points at the real data/
% through a symbolic link made by MATLAB's system() (no shell session needed). Copying it instead would be about 59 GB.
[st, out] = system(sprintf('ln -s "%s" "%s"', fullfile(proj, 'data'), fullfile(sh, 'data')));
if st ~= 0, error('N48rds:link', '%s', "Could not link data/: " + string(out)); end
logL = "N=48 RDS copy built " + string(datetime('now')) + " from " + string(proj);
% 1. move N=45 Fraser data files aside (before the swaps, so the N=48 inputs are not moved)
keepNames = ["noiseCharacterisation_fraserN48.mat" "constellationFraser_v002.mat"];
for top = ["src" "results"]
    L = dir(fullfile(sh, top, '**', '*'));
    for k = 1:numel(L)
        nm = string(L(k).name);  p = fullfile(L(k).folder, L(k).name);
        if ~contains(lower(nm), "fraser") || contains(nm, "FraserN48") || ismember(nm, keepNames) || ~isfolder(L(k).folder), continue, end
        isData = ~L(k).isdir && any(endsWith(lower(nm), [".mat" ".txt" ".csv" ".log"]));
        isChunks = L(k).isdir && endsWith(lower(nm), "_chunks");
        if ~(isData || isChunks) || ~(isfile(p) || isfolder(p)), continue, end
        rel = extractAfter(string(p), string(sh) + "/");  dst = fullfile(sh, 'N45_moved', rel);
        if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
        ok_local(movefile(p, dst), "move " + rel);  logL(end + 1) = "moved " + rel; %#ok<AGROW>
    end
end
% 2. swaps
S = fullfile(sh, 'src');
sw = ["loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat" "loopClosureResults_Fraser_all_shaped_xu_v015.mat";
      "loopClosureResults_FraserN48_all_shaped_xu_v020.mat"        "loopClosureResults_Fraser_all_shaped_xu_v019.mat";
      "noiseCharacterisation_fraserN48.mat"                        "noiseCharacterisation_fraser.mat";
      "constellationFraser_v002.mat"                               "constellationFraser_v001.mat"];
for k = 1:height(sw)
    ok_local(copyfile(fullfile(S, sw(k, 1)), fullfile(S, sw(k, 2))), "swap " + sw(k, 2));
    logL(end + 1) = "swap src/" + sw(k, 2) + " <- " + sw(k, 1) + "  md5 " + md5_local(fullfile(S, sw(k, 2))); %#ok<AGROW>
end
shim = fullfile(S, 'functions', 'importDB_fraser_v001.m');
writelines(["function [trials, rejects, info] = importDB_fraser_v001(varargin)"
    "% N=48 SHADOW SHIM (buildN48ShadowTree_RDS_v001): calls importDB_fraser_v002. 'Shapes' is dropped (v002 is Ellipse-only)."
    "args = varargin;  keep = true(size(args));"
    "for i = 1:numel(args), if (ischar(args{i}) || isstring(args{i})) && strcmpi(args{i}, 'Shapes'), keep([i i+1]) = false; end, end"
    "args = args(keep);"
    "[trials, rejects, info] = importDB_fraser_v002(args{:});"
    "end"], shim);
logL(end + 1) = "shim src/functions/importDB_fraser_v001.m -> importDB_fraser_v002";
M = dir(fullfile(S, '**', '*.m'));  nRep = 0;
for k = 1:numel(M)
    p = fullfile(M(k).folder, M(k).name);  t = fileread(p);
    if contains(t, "10.41793")
        fid = fopen(p, 'w');  fwrite(fid, replace(t, "10.41793", "5.208965"));  fclose(fid);
        nRep = nRep + 1;  logL(end + 1) = "scale literal replaced in " + extractAfter(string(p), string(sh) + "/"); %#ok<AGROW>
    end
end
writelines(logL(:), fullfile(sh, 'N48_SWAPS.txt'));
fprintf('\nN=48 copy ready: %s\n  %d N=45 Fraser files moved to N45_moved/, %d swaps, scale literal replaced in %d .m files, %.1f min\n', ...
    sh, sum(startsWith(logL, "moved")), height(sw) + 1, nRep, toc(t0) / 60);
fprintf('Next: cd to %s and open the launchers there (sessions A-E).\n=== COMPLETED OK ===\n', fullfile(sh, 'src'));
end

function ok_local(st, what)
    if ~st, error('N48rds:io', '%s', "Failed: " + what); end
end
function m = md5_local(f)
    [st, out] = system(sprintf('md5sum "%s"', f));
    if st == 0, m = string(strtok(out)); else, m = "(md5sum unavailable)"; end
end
