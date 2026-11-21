% runK4Regression_v001 (2026-11-21): v002 builders at K=4 must reproduce the existing K=4 outputs exactly. Stops (error) on any difference.
B2B = fullfile(getenv('HOME'),'Library','CloudStorage','Dropbox','Brain2Bee');
old = fullfile(B2B,'PowerLawSimulationPreReg');
od  = fullfile(B2B,'PowerLaw_Repro_2026','BRM','outputs_K1','regression_K4');
if ~isfolder(od), mkdir(od); end
addpath(fileparts(mfilename('fullpath')));
tic; loopClosureVarDecomp_v015(DataSource="v028", NHarm=4, OutDir=od, OutTag="K4regress"); fprintf('varDecomp K=4 run time %.1f s\n', toc);
newVd = fullfile(od,'loopClosureVarDecomp_v015_all6_v028_K4regress.mat');
A = load(fullfile(old,'src','loopClosureVarDecomp_v014_gated_all6_v028.mat'));  N = load(newVd);
fn = setdiff(fieldnames(A), 'meta');  bad = {};
for i = 1:numel(fn), if ~isequaln(A.(fn{i}), N.(fn{i})), bad{end+1} = fn{i}; end, end %#ok<SAGROW>
fm = fieldnames(A.meta);  badm = {};
for i = 1:numel(fm), if ~strcmp(fm{i},'createdAt') && ~isequaln(A.meta.(fm{i}), N.meta.(fm{i})), badm{end+1} = fm{i}; end, end %#ok<SAGROW>
fprintf('varDecomp: fields compared %s; differing: [%s]; meta differing (createdAt ignored): [%s]\n', strjoin(fn,','), strjoin(bad,','), strjoin(badm,','));
assert(isempty(bad) && isempty(badm), 'K4 regression: varDecomp differs');
T = buildTable6bV028_v002(NHarm=4, VdFile=newVd, OutDir=od, OutTag="K4regress");
P = load(fullfile(old,'results','buildTable6bV028_v001_v028.mat')).T;
assert(isequaln(T(:,1:13), P), 'K4 regression: Table 6b differs');
a = readlines(fullfile(old,'results','buildTable6bV028_v001_v028.csv'));  b = readlines(fullfile(od,'buildTable6bV028_v002_K4regress.csv'));
b = regexprep(b, ',4$', '');  b(1) = regexprep(b(1), ',nHarm$', '');   % nHarm is the only added column
fprintf('Table 6b: table isequaln=1; csv text identical after dropping the added nHarm column: %d\n', isequal(a, b));
assert(isequal(a, b), 'K4 regression: csv text differs');
fprintf('K4 REGRESSION PASSED\n');
