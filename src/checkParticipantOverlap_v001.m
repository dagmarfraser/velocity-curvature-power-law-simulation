%% checkParticipantOverlap_v001.m
% D-OV (TODO_PreSubmission): could any participant appear in more than one analysed dataset?
% Reads subjectID from the v015 per-trial corpora and prints, per dataset, the number of
% subjects and the ID scheme, then every ID shared between datasets. A shared ID is only a
% flag: the source studies numbered participants independently, so equal IDs need not be the
% same person, and different IDs cannot rule overlap out. The answer comes from the source
% studies' recruitment records; this script states what the corpora alone can show.
% Run from src/. Writes logs/checkParticipantOverlap_v001_<date>.log.

%% CONFIG
DS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"];
VER = "v015";

%% Run
SRC = fileparts(mfilename("fullpath"));  ROOT = fileparts(SRC);
logFile = fullfile(ROOT, "logs", "checkParticipantOverlap_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("overlap:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkParticipantOverlap_v001 run %s\n", string(datetime("now")));
ids = cell(numel(DS), 1);
for d = 1:numel(DS)
    f = fullfile(SRC, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", DS(d), VER));
    if ~isfile(f), fclose(LG); error("overlap:missing", "%s", "Missing corpus: " + f); end
    S = load(f, "results");
    s = string({S.results.subjectID});
    if any(s == "" | ismissing(s)), fclose(LG); error("overlap:emptyID", "%s", "Empty subjectID in " + DS(d)); end
    ids{d} = unique(s);
    pat = unique(regexprep(ids{d}, "\d", "#"));
    say_local(LG, "%-13s %4d subjects, %5d trials; ID pattern(s): %s; e.g. %s\n", DS(d), numel(ids{d}), ...
        numel(s), strjoin(pat, ", "), strjoin(ids{d}(1:min(3, end)), ", "));
end
say_local(LG, "\nIDs shared between datasets (a flag only; see header):\n");
nShared = 0;
for a = 1:numel(DS)
    for b = a+1:numel(DS)
        c = intersect(ids{a}, ids{b});
        if ~isempty(c)
            nShared = nShared + numel(c);
            say_local(LG, "  %s & %s: %d shared (%s)\n", DS(a), DS(b), numel(c), strjoin(c(1:min(8, end)), ", "));
        end
    end
end
if nShared == 0, say_local(LG, "  none\n"); end
say_local(LG, "Log: %s\n", logFile);  fclose(LG);

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
