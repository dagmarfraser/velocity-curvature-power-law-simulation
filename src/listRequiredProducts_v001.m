%% listRequiredProducts_v001.m
% Which MathWorks products do the paper's load-bearing scripts need? (TODO_PreSubmission A11)
% Runs matlab.codetools.requiredFilesAndProducts on each script below (all its dependencies are
% followed) and prints, per script and in union, the products and versions on THIS machine.
% Versions printed are the installed ones here, not necessarily the release a given result was
% run on: BlueBEAR runs used R2024b (module MATLAB/2024b), the Macs R2025b/R2026a.
% Fails loud on a missing script. Run from src/. Writes logs/listRequiredProducts_v001_<date>.log
% and results/listRequiredProducts_v001.mat.

%% CONFIG
SCRIPTS = ["runLoopClosureFftnoise_v015", "Toolchain_caller_v058", "Toolchain_func_v032", ...
    "ModelAdequacy_Stage1_KitchenSink_v2_001", "ModelAdequacy_Stage2_Assessment_v003", ...
    "functions/characteriseBiologicalNoise_v004", "computePerCoordinateSEM_v2_001", ...
    "loopClosureVarDecomp_v012_gated", "compareCookCtrlAsd_v002", "checkAdditiveVsRetimed_v003", ...
    "analyseTrialTempoWindow_v002", "checkPart3PrimaryMetric_v001", "bootstrapTempoReduction_v001", ...
    "fitTempoRandomSlope_v001", "simulateTempoCorpus_v002", "predictHarmonicGain_v003", ...
    "checkFilterCutoffCollapse_v002"];

%% Run
SRC  = fileparts(mfilename("fullpath"));  ROOT = fileparts(SRC);
addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
logFile = fullfile(ROOT, "logs", "listRequiredProducts_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("reqProd:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "listRequiredProducts_v001 run %s on %s, MATLAB %s\n", string(datetime("now")), string(getenv("HOME")), version);
rows = table();
for k = 1:numel(SCRIPTS)
    f = fullfile(SRC, SCRIPTS(k) + ".m");
    if ~isfile(f), fclose(LG); error("reqProd:missing", "%s", "Missing script: " + f); end
    [~, pl] = matlab.codetools.requiredFilesAndProducts(char(f));
    for p = 1:numel(pl)
        rows = [rows; table(SCRIPTS(k), string(pl(p).Name), string(pl(p).Version), ...
            'VariableNames', ["script" "product" "version"])]; %#ok<AGROW>
    end
    say_local(LG, "%-45s %s\n", SCRIPTS(k), strjoin(string({pl.Name}), "; "));
end
U = unique(rows(:, ["product" "version"]), "rows");
say_local(LG, "\nUnion of products (installed versions on this machine):\n%s", evalc("disp(U)"));
save(fullfile(ROOT, "results", "listRequiredProducts_v001.mat"), "rows", "U", "SCRIPTS");
say_local(LG, "Log: %s\n", logFile);  fclose(LG);

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
