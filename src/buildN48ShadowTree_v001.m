function shadowRoot = buildN48ShadowTree_v001(opts)
%BUILDN48SHADOWTREE_V001  A copy of the project in which every Fraser input is the N=48 one, for the Part E reruns (2026-11-04).
%   WHY: about 25 Part E scripts read the Fraser N=45 corpus by name (loopClosureResults_Fraser_all_shaped_xu_v015.mat,
%   noiseCharacterisation_fraser.mat, constellationFraser_v001.mat, importDB_fraser_v001) and write results under fixed names.
%   Rather than hand-edit 25 new versions, each is run UNCHANGED inside a shadow copy where those Fraser inputs ARE the N=48
%   ones. Nothing in the real project is read-modified or overwritten: the copy is an APFS clone (cp -c, no extra disk until
%   written) outside Dropbox, and every output lands in the shadow. Final outputs are copied back under results/N48/.
%   SWAPS (and only these), each recorded in <shadow>/N48_SWAPS.txt with md5s:
%     src/loopClosureResults_Fraser_all_shaped_xu_v015.mat  <- loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat (production)
%     src/loopClosureResults_Fraser_all_shaped_xu_v019.mat  <- loopClosureResults_FraserN48_all_shaped_xu_v020.mat (end-anchored)
%     src/noiseCharacterisation_fraser.mat                  <- noiseCharacterisation_fraserN48.mat
%     src/constellationFraser_v001.mat                      <- constellationFraser_v002.mat
%     src/functions/importDB_fraser_v001.m                  <- shim calling importDB_fraser_v002 (Shapes ignored: v002 is Ellipse-only)
%     every src/**/*.m: the literal 10.41793 (pixel density) -> 5.208965 (ipadScaleM4_v001, UIKit points per mm)
%   data/ is a symlink to the real data (read-only use). Downstream intermediates (varDecomp, gate tables, ...) are NOT swapped:
%   they are regenerated inside the shadow by running their producers first.
%   USAGE: root = buildN48ShadowTree_v001();   then run scripts from <root>/src.   Fraser, D.S. (2026)  v001
arguments
    opts.ShadowRoot (1,1) string = fullfile(getenv('HOME'), 'N48shadow', 'PowerLawSimulationPreReg')
end
srcDir = fileparts(mfilename('fullpath'));  proj = fileparts(srcDir);  shadowRoot = opts.ShadowRoot;
if isfolder(shadowRoot), error('N48shadow:exists', '%s', "Shadow exists: " + shadowRoot + " (move it away to rebuild)."); end
mkdir(fileparts(shadowRoot));  mkdir(shadowRoot);
sh = @(c) run_local(c);
sh(sprintf('cp -c -R "%s" "%s/src"', srcDir, shadowRoot));
sh(sprintf('cp -c -R "%s" "%s/results"', fullfile(proj, 'results'), shadowRoot));
sh(sprintf('ln -s "%s" "%s/data"', fullfile(proj, 'data'), shadowRoot));
S = fullfile(shadowRoot, 'src');
swaps = ["loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat" "loopClosureResults_Fraser_all_shaped_xu_v015.mat";
         "loopClosureResults_FraserN48_all_shaped_xu_v020.mat"        "loopClosureResults_Fraser_all_shaped_xu_v019.mat";
         "noiseCharacterisation_fraserN48.mat"                        "noiseCharacterisation_fraser.mat";
         "constellationFraser_v002.mat"                               "constellationFraser_v001.mat"];
log = strings(0);
for k = 1:height(swaps)
    from = fullfile(srcDir, swaps(k, 1));  to = fullfile(S, swaps(k, 2));
    if ~isfile(from), error('N48shadow:missing', '%s', "Missing N=48 input: " + from); end
    delete(to);  copyfile(from, to);
    log(end + 1) = swaps(k, 2) + " <- " + swaps(k, 1) + "  md5 " + md5_local(to); %#ok<AGROW>
end
shim = fullfile(S, 'functions', 'importDB_fraser_v001.m');  delete(shim);
fid = fopen(shim, 'w');
fprintf(fid, '%s\n', ...
 "function [trials, rejects, info] = importDB_fraser_v001(varargin)", ...
 "% N=48 SHADOW SHIM (buildN48ShadowTree_v001): calls importDB_fraser_v002. 'Shapes' is dropped (v002 is Ellipse-only).", ...
 "args = varargin;  keep = true(size(args));", ...
 "for i = 1:numel(args), if (ischar(args{i}) || isstring(args{i})) && strcmpi(args{i}, 'Shapes'), keep([i i+1]) = false; end, end", ...
 "args = args(keep);", ...
 "[trials, rejects, info] = importDB_fraser_v002(args{:});", ...
 "end");
fclose(fid);
log(end + 1) = "functions/importDB_fraser_v001.m <- shim calling importDB_fraser_v002";
files = dir(fullfile(S, '**', '*.m'));  nRep = 0;  repFiles = strings(0);
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);  t = fileread(p);
    if contains(t, "10.41793")
        t = replace(t, "10.41793", "5.208965");  fid = fopen(p, 'w');  fwrite(fid, t);  fclose(fid);
        nRep = nRep + 1;  repFiles(end + 1) = erase(string(p), string(S) + filesep); %#ok<AGROW>
    end
end
log(end + 1) = sprintf("10.41793 -> 5.208965 in %d .m files: %s", nRep, strjoin(repFiles, ", "));
writelines(["N=48 shadow built " + string(datetime('now')) + " by buildN48ShadowTree_v001 from " + string(proj); log(:)], fullfile(shadowRoot, 'N48_SWAPS.txt'));
fprintf('Shadow ready: %s (%d swaps, scale literal replaced in %d files)\n', shadowRoot, height(swaps) + 1, nRep);
end

function run_local(c)
    [st, out] = system(c);
    if st ~= 0, error('N48shadow:system', '%s', "Command failed: " + c + newline + out); end
end
function m = md5_local(f)
    [~, m] = system(sprintf('md5 -q "%s"', f));  m = strtrim(m);
end
