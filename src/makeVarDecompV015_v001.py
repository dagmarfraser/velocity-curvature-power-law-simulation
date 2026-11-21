# Derives BRM/scripts/loopClosureVarDecomp_v015.m from the untouched v014_gated by exact replacements, each asserted to match once (2026-11-21).
import os
B2B = os.path.expanduser('~/Library/CloudStorage/Dropbox/Brain2Bee')
src = open(B2B + '/PowerLawSimulationPreReg/src/loopClosureVarDecomp_v014_gated.m').read()

def rep(a, b):
    global src
    assert src.count(a) == 1, (a, src.count(a))
    src = src.replace(a, b)

rep('function out = loopClosureVarDecomp_v014_gated(opts)\n',
    'function out = loopClosureVarDecomp_v015(opts)\n'
    '% loopClosureVarDecomp_v015 (2026-11-21): fork of v014_gated; v014_gated is untouched. Changes, and only these:\n'
    '%   opts.NHarm (K; 4 = the primary arm as in v014; 0,1,2,3,6,10 are read through loadV028_v002 from BRM/scratch_K or BRM/outputs_K0, md5 checked),\n'
    '%   opts.OutDir (REQUIRED: output folder; the old src/ is never a default), opts.OutTag (appended to the output file name),\n'
    '%   opts.SrcDir (old-tree src, default fixed path), loadV028_v001 -> loadV028_v002, meta.nHarm, and an error if the output file already exists.\n'
    '%   Cluster-bootstrap/synthesis machinery and the RNG stream (global twister, seed 20260619) are unchanged.\n'
    '% ---- original v014_gated header follows ----\n')

rep('        opts.ObsBasis       (1,1) string {mustBeMember(opts.ObsBasis, ["window","whole"])} = "window"                          % v014: DataSource "v028" only\n    end\n',
    '        opts.ObsBasis       (1,1) string {mustBeMember(opts.ObsBasis, ["window","whole"])} = "window"                          % v014: DataSource "v028" only\n'
    '        opts.NHarm          (1,1) double {mustBeMember(opts.NHarm, [0 1 2 3 4 6 10])} = 4                                       % v015: DataSource "v028", Arm "primary" only when ~= 4\n'
    '        opts.OutDir         (1,1) string = ""                                                                                   % v015: required\n'
    '        opts.OutTag         (1,1) string = ""                                                                                   % v015\n'
    '        opts.SrcDir         (1,1) string = ""                                                                                   % v015\n'
    '    end\n'
    '    if opts.OutDir == "", error("loopClosureVarDecomp_v015:outdir", "OutDir is required (v015 never writes into the old src/ by default)"); end\n'
    '    if ~isfolder(opts.OutDir), error("loopClosureVarDecomp_v015:outdir", "OutDir %s does not exist", opts.OutDir); end\n'
    '    if opts.NHarm ~= 4 && (opts.DataSource ~= "v028" || opts.Arm ~= "primary"), error("loopClosureVarDecomp_v015:nharm", "NHarm ~= 4 needs DataSource v028 and Arm primary"); end\n')

rep('    srcDir  = fileparts(mfilename("fullpath"));  addpath(genpath(fullfile(srcDir, "functions")));   % v014: loadV028_v001\n',
    '    srcDir  = opts.SrcDir;  if srcDir == "", srcDir = fullfile(getenv("HOME"), "Library", "CloudStorage", "Dropbox", "Brain2Bee", "PowerLawSimulationPreReg", "src"); end   % v015\n'
    '    if ~isfolder(srcDir), error("loopClosureVarDecomp_v015:srcdir", "SrcDir %s does not exist", srcDir); end\n'
    '    addpath(genpath(fullfile(srcDir, "functions")));  addpath(fileparts(mfilename("fullpath")));   % v015: old-tree helpers, and loadV028_v002 beside this file\n')

rep('    outFile = fullfile(srcDir, sprintf("loopClosureVarDecomp_v014_gated_%s_%s.mat", subsetTag, srcTag));\n',
    '    if opts.NHarm ~= 4, srcTag = srcTag + "_nH" + opts.NHarm; end   % v015\n'
    '    if opts.OutTag ~= "", srcTag = srcTag + "_" + opts.OutTag; end   % v015\n'
    '    outFile = fullfile(opts.OutDir, sprintf("loopClosureVarDecomp_v015_%s_%s.mat", subsetTag, srcTag));\n'
    '    if isfile(outFile), error("loopClosureVarDecomp_v015:exists", "output %s already exists; v015 does not overwrite", outFile); end\n')

rep('        r = loadResults_local(srcDir, tag, opts.DataSource, opts.Arm, opts.ObsBasis);\n',
    '        r = loadResults_local(srcDir, tag, opts.DataSource, opts.Arm, opts.ObsBasis, opts.NHarm);\n')
rep('function r = loadResults_local(srcDir, tag, dataSource, arm, obsBasis)\n',
    'function r = loadResults_local(srcDir, tag, dataSource, arm, obsBasis, nHarm)\n')
rep('        r = loadV028_v001(strrep(tag, "_", " "), Arm=arm, ObsBasis=obsBasis);  r = r(:);  return\n',
    '        r = loadV028_v002(strrep(tag, "_", " "), Arm=arm, ObsBasis=obsBasis, NHarm=nHarm, SrcDir=srcDir);  r = r(:);  return   % v015\n')
rep('        "isDefaultSubset", isDefaultSubset, ...\n',
    '        "isDefaultSubset", isDefaultSubset, "nHarm", opts.NHarm, "arm", opts.Arm, "obsBasis", opts.ObsBasis, ...\n')

open(os.path.dirname(os.path.abspath(__file__)) + '/loopClosureVarDecomp_v015.m', 'w').write(src)
print('written; lines', src.count('\n'))
