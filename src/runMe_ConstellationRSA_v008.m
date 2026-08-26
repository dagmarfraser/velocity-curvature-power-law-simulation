%RUNME_CONSTELLATIONRSA_V008 Path setup + launch for the Pilot-retirement
% RSA/Mantel rerun (runConstellationRSA_HPC_v008.m) on BlueBEAR.
%
% Run this from an interactive MATLAB session started via sinteractive
% (never matlab -nodisplay batch, per this project's own convention) --
% this script does NOT submit a SLURM job itself, it just handles the
% two addpath calls v008 needs (its own src/ folder, plus
% src/mantelToolbox/, which v008 deliberately no longer adds
% automatically -- see v008's own header note) and then runs it.
%
% src/mantelToolbox/mantelTest.m was confirmed 2026-08-26 (direct `diff`
% against stagingForGithub/mantel/toolbox/mantelTest.m, empty output) to
% be byte-identical to the currently published mantel toolbox, including
% the 2026-07-27 tie-handling fix -- genuinely current, not a stale fork.
%
% Usage (from within an sinteractive MATLAB session, cd'd into src/, or
% anywhere -- path resolution below is relative to this file's own
% location, not the current working directory):
%   run('runMe_ConstellationRSA_v008.m')
% or, if already cd'd into src/:
%   runMe_ConstellationRSA_v008
%
% Fraser, D.S. (2026)

srcDir = fileparts(mfilename("fullpath"));
mantelDir = fullfile(srcDir, "mantelToolbox");

fprintf("srcDir:    %s\n", srcDir);
fprintf("mantelDir: %s\n\n", mantelDir);

if ~isfolder(mantelDir)
    error("runMe_ConstellationRSA_v008:mantelDirNotFound", "%s", ...
        sprintf("mantelToolbox folder not found at %s -- confirm this is the correct RDS project path.", mantelDir));
end

addpath(srcDir);
addpath(mantelDir);

if exist("mantelTest", "file") ~= 2
    error("runMe_ConstellationRSA_v008:mantelTestNotFound", "%s", ...
        "mantelTest.m still not found on path after addpath -- check mantelToolbox/ actually contains it.");
end
if exist("buildConstellationRDM_v002", "file") ~= 2
    error("runMe_ConstellationRSA_v008:builderNotFound", "%s", ...
        "buildConstellationRDM_v002.m not found on path after addpath(srcDir).");
end

fprintf("Path confirmed (mantelTest, buildConstellationRDM_v002 both resolve).\n");
fprintf("Launching runConstellationRSA_HPC_v008...\n\n");

runConstellationRSA_HPC_v008
