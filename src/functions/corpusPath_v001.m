function f = corpusPath_v001(name, opts)
%CORPUSPATH_V001  The loop-closure corpus file for a dataset, with Fraser on the N=48 corpus (2026-11-04).
%   f = corpusPath_v001("Cook_CTRL")    -> src/loopClosureResults_Cook_CTRL_all_shaped_xu_v015.mat
%   f = corpusPath_v001("Fraser")       -> src/loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat
%   f = corpusPath_v001("Fraser", Template="endAnchored")  -> ..._FraserN48_all_shaped_xu_v020.mat (D30 sensitivity)
%   Why: the N=48 Fraser corpus (importDB_fraser_v002, 5.209 pt/mm, BlueBEAR 2026-11-04, runLoopClosureFftnoise_v020) replaces
%   the N=45 v015 file for Fraser. Its 'production' template is v015's numerics (Z12 reporting basis; checkV020Acceptance_v001).
%   The other datasets keep their v015 files. With Template="endAnchored" every dataset gets its v019 file and Fraser
%   the end-anchored v020 file. Fails loud if the file is missing. Scripts versioned from 2026-11-04 call this instead of
%   building the name, so the Fraser corpus is chosen in one place.   Fraser, D.S. (2026)  v001
arguments
    name (1,1) string {mustBeMember(name, ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"])}
    opts.Template (1,1) string {mustBeMember(opts.Template, ["production","endAnchored"])} = "production"
end
srcDir = fileparts(fileparts(mfilename("fullpath")));
if name == "Fraser"
    if opts.Template == "production", f = fullfile(srcDir, "loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat");
    else,                             f = fullfile(srcDir, "loopClosureResults_FraserN48_all_shaped_xu_v020.mat"); end
else
    if opts.Template == "production", f = fullfile(srcDir, "loopClosureResults_" + name + "_all_shaped_xu_v015.mat");
    else,                             f = fullfile(srcDir, "loopClosureResults_" + name + "_all_shaped_xu_v019.mat"); end
end
if ~isfile(f), error("corpusPath:missing", "%s", "FAILED PATH: " + f); end
end
