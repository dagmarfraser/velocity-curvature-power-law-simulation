function T = buildTable6bV028_v002(opts)
%BUILDTABLE6BV028_V002 Table 6b on the v028 corpus with a selectable K (2026-11-21). Fork of buildTable6bV028_v001 (untouched).
%   Changes, and only these: opts.NHarm (K; 4 = primary arm as in v001), opts.VdFile (the varDecomp .mat to read, REQUIRED; the v015 varDecomp output
%   of the same NHarm/Arm/ObsBasis), opts.OutDir (REQUIRED), opts.OutTag (REQUIRED, in the output file names), opts.SrcDir, loader loadV028_v001 -> loadV028_v002,
%   the "v015" basis branch (regression anchor to the v015 corpus) is dropped (v028 only), and the output is never overwritten.
%   Rows, band rule (1/3 +/- MDC 0.03), columns and statistics are v001's. flat/cluster/ciT come from the varDecomp file; mean and SD are over trials with a
%   finite betaGenStarMed. The varDecomp file must have been built with the same NHarm/Arm/ObsBasis (meta checked; error otherwise).
%   Saves OutDir/buildTable6bV028_v002_<OutTag>.mat and .csv.
arguments
    opts.Arm (1,1) string {mustBeMember(opts.Arm, ["primary", "budget1", "shapeXu", "harmonicsFc"])} = "primary"
    opts.NHarm (1,1) double {mustBeMember(opts.NHarm, [0 1 2 3 4 6 10])} = 4
    opts.ObsBasis (1,1) string {mustBeMember(opts.ObsBasis, ["window", "whole"])} = "window"
    opts.VdFile (1,1) string = ""
    opts.OutDir (1,1) string = ""
    opts.OutTag (1,1) string = ""
    opts.SrcDir (1,1) string = ""
end
if opts.VdFile == "" || opts.OutDir == "" || opts.OutTag == "", error('tbl6bV028v2:args', 'VdFile, OutDir and OutTag are required'); end
if ~isfolder(opts.OutDir), error('tbl6bV028v2:outdir', 'OutDir %s does not exist', opts.OutDir); end
srcDir = opts.SrcDir;  if srcDir == "", srcDir = fullfile(getenv('HOME'), 'Library', 'CloudStorage', 'Dropbox', 'Brain2Bee', 'PowerLawSimulationPreReg', 'src'); end
if ~isfolder(srcDir), error('tbl6bV028v2:srcdir', 'SrcDir %s does not exist', srcDir); end
addpath(genpath(fullfile(srcDir, 'functions')));  addpath(fileparts(mfilename('fullpath')));
outMat = fullfile(opts.OutDir, "buildTable6bV028_v002_" + opts.OutTag + ".mat");  outCsv = fullfile(opts.OutDir, "buildTable6bV028_v002_" + opts.OutTag + ".csv");
if isfile(outMat) || isfile(outCsv), error('tbl6bV028v2:exists', 'output %s or %s already exists; not overwritten', outMat, outCsv); end
ds  = ["Fraser" "Cook ASD" "Cook CTRL" "Hickman HALO" "Hickman PLAC"];
MDC = 0.03;  third = 1/3;  lo0 = third - MDC;  hi0 = third + MDC;
vdFile = opts.VdFile;
assert(isfile(vdFile), 'tbl6bV028v2:src', 'varDecomp file missing: %s', vdFile);
V = load(vdFile);  assert(V.meta.dataSource == "v028", 'tbl6bV028v2:src', 'varDecomp file is not the v028 source');
if isfield(V.meta, 'nHarm')   % v015 files carry nHarm/arm/obsBasis; a v014 file (K=4 primary window only) does not
    assert(V.meta.nHarm == opts.NHarm && V.meta.arm == opts.Arm && V.meta.obsBasis == opts.ObsBasis, 'tbl6bV028v2:meta', ...
        'varDecomp file was built with NHarm=%g Arm=%s ObsBasis=%s, not the requested NHarm=%g Arm=%s ObsBasis=%s', V.meta.nHarm, V.meta.arm, V.meta.obsBasis, opts.NHarm, opts.Arm, opts.ObsBasis);
else
    assert(opts.NHarm == 4 && opts.Arm == "primary" && opts.ObsBasis == "window", 'tbl6bV028v2:meta', 'a v014 varDecomp file is the K=4 primary window source only');
end
rows = cell(0, 13);
for k = 1:numel(ds)
    i = find(string({V.out.label}) == ds(k));  assert(isscalar(i), 'tbl6bV028v2:label', 'dataset %s not in %s', ds(k), vdFile);
    cm = V.out(i).cluster_median;  fl = V.out(i).flat_median;
    R = loadV028_v002(ds(k), Arm=opts.Arm, ObsBasis=opts.ObsBasis, NHarm=opts.NHarm, SrcDir=srcDir);
    x = arrayfun(@(r) double(r.betaGenStarMed), R(:));  x = x(isfinite(x));
    ci = cm.ciT;  inBand = ci(1) >= lo0 && ci(2) <= hi0;  outBand = ci(2) < lo0 || ci(1) > hi0;
    status = "inconclusive";  if inBand, status = "indistinguishable"; elseif outBand, status = "different"; end
    rows(end+1, :) = {ds(k), fl.median, cm.pointEst, ci(1), ci(2), mean(x), std(x), cm.N_trials, numel(x), ci(1) <= third && third <= ci(2), status, ...
        cm.pointEst - third, "v028"}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'dataset', 'flat', 'cluster', 'ciLo', 'ciHi', 'mean', 'sd', 'Nfin', 'nFiniteMed', 'ciIncludesOneThird', 'bandStatus', 'clusterMinusOneThird', 'basis'});
T.dataset = string(T.dataset);  T.bandStatus = string(T.bandStatus);  T.basis = string(T.basis);
T.nHarm = repmat(opts.NHarm, height(T), 1);
disp(T);
save(outMat, 'T');  writetable(T, outCsv);
end
