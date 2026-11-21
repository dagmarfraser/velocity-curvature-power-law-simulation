% runK1Rebuild_v001 (2026-11-21): variance decomposition and Table 6b at K=1 (nH1 files, BRM/scratch_K, md5 checked), outputs in BRM/outputs_K1 with _K1 in the names.
B2B = fullfile(getenv('HOME'),'Library','CloudStorage','Dropbox','Brain2Bee');
od  = fullfile(B2B,'PowerLaw_Repro_2026','BRM','outputs_K1');
addpath(fileparts(mfilename('fullpath')));
tic; loopClosureVarDecomp_v015(DataSource="v028", NHarm=1, OutDir=od, OutTag="K1"); fprintf('varDecomp K=1 run time %.1f s\n', toc);
T = buildTable6bV028_v002(NHarm=1, VdFile=fullfile(od,'loopClosureVarDecomp_v015_all6_v028_nH1_K1.mat'), OutDir=od, OutTag="K1");
