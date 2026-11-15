%% auditLimitBreak_v002.m
% v002 (2026-11-03): v001 + the N=48 call sites in the paper-feeding list (runLoopClosureFftnoise_v015, v019, v020 and
% auditFraserN48ObservedFits_v001), and a hard error, not a warning, when a paper-feeding site passes anything but 0.
% Writes results/limitBreakAudit_v002.mat; v001 and its mat are left as they were.
%
% Audit every regressDataEBR call site in src/ and tool/ for its limitBreak argument.
% regressDataEBR (L99-118, L160-183): limitBreak = 0 -> NaN on ANY fitnlm warning (incl. the
% iteration limit), so no unconverged LMLS/IRLS estimate is used; limitBreak = 1 keeps the
% final-iteration estimate. OLS-type calls (regressType 1-3) never read the argument.
% Parsing: continuation lines joined; arguments split on top-level commas only (bracket depth),
% so regressConfigs(rIdx).type and [a, b] seeds are handled. Obsolete/archive folders skipped.
% Writes results/limitBreakAudit_v002.mat (Finding #232).

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
DIRS    = ["src" "tool"];
SKIP    = ["obsolete" "archive"];
OUT_MAT = fullfile(ROOT, "results", "limitBreakAudit_v002.mat");
PAPER   = ["Toolchain_func_v032" "runLoopClosureFftnoise_v012" "runLoopClosureFftnoise_v013" ...
           "runLoopClosureFftnoise_v014" "constellationCCC_v007" "constellationFraser_v001" ...
           "constellationCook_v002" "constellationHickman_v001" "simexBetaRecovery_v001" ...
           "measureRegressionSaturation_v003" "saturationSweepEngine_v001" "tempoFoldEngine_v001" ...
           "processTrialLoopClosure_v001" ...
           "runLoopClosureFftnoise_v015" "runLoopClosureFftnoise_v019" "runLoopClosureFftnoise_v020" ...
           "auditFraserN48ObservedFits_v001"];   % v002: the N=48 runners;   % paper-feeding call sites (Finding #232)

%% Scan
files = [];
for d = DIRS, files = [files; dir(fullfile(ROOT, d, "**", "*.m"))]; end %#ok<AGROW>
files = files(~contains(string({files.folder}), SKIP) & ~startsWith(string({files.name}), "auditLimitBreak"));   % not itself
A = table();
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    L = splitlines(string(fileread(p)));
    hit = find(contains(L, "regressDataEBR(") & ~startsWith(strtrim(L), "%") & ...
               ~contains(L, ["function " "sprintf(" "failed:"]));
    for j = hit(:)'
        s = L(j);  n = j;
        while contains(s, "...") && n < numel(L), n = n + 1; s = extractBefore(s, "...") + " " + strtrim(L(n)); end
        args = splitArgs_local(char(extractAfter(s, "regressDataEBR(")));
        rt = "";  lb = "(none)";
        if numel(args) >= 3, rt = string(args{3}); end
        if numel(args) >= 6, lb = string(args{6}); end
        A = [A; table(erase(string(p), string(ROOT) + filesep), j, numel(args), rt, lb, ...
            'VariableNames', ["file" "line" "nArgs" "regressType" "limitBreak"])]; %#ok<AGROW>
    end
end
[~, stem] = fileparts(A.file);  A.paperFeeding = ismember(stem, PAPER);

%% Report
fprintf("%d call sites in %d files\n", height(A), numel(unique(A.file)));
disp(groupcounts(A, "limitBreak"))
P = A(A.paperFeeding, :);
fprintf("Paper-feeding call sites: %d; limitBreak values: %s\n", height(P), strjoin(unique(P.limitBreak)', ", "));
bad = P(~ismember(P.limitBreak, ["0" "limitBreak"]), :);   % v002: as v001 (variable-valued sites resolved in the caller), but a hard stop
if height(bad) > 0
    disp(bad)
    error("auditLimitBreak_v002:paper", "%s", sprintf("%d paper-feeding call site(s) do not pass limitBreak = 0 (listed above).", height(bad)));
end
missing = setdiff(PAPER, unique(string(stem(A.paperFeeding))));   % v002: every listed file must be found
if ~isempty(missing), error("auditLimitBreak_v002:missing", "%s", "Paper-feeding files with no regressDataEBR call found: " + strjoin(missing, ", ")); end
noArgNonlinear = A(A.limitBreak == "(none)" & ~ismember(A.regressType, ["1" "2" "3"]), :);
fprintf("Calls without limitBreak on a non-OLS type: %d\n", height(noArgNonlinear));
if height(noArgNonlinear) > 0, disp(noArgNonlinear), end
fprintf("\nlimitBreak = 1 (unconverged fits kept):\n");  disp(A(A.limitBreak == "1", ["file" "line" "regressType"]))
fprintf("Variable-valued limitBreak (resolve in the caller):\n");  disp(A(A.limitBreak == "limitBreak", ["file" "line"]))

if ~isfolder(fileparts(OUT_MAT)), error("auditLimitBreak:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "A", "PAPER", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function args = splitArgs_local(c)
    % Split a call's argument text on top-level commas; stop at the closing parenthesis.
    d = 1;  buf = '';  args = {};
    for ch = c
        if any(ch == '([{'), d = d + 1; elseif any(ch == ')]}'), d = d - 1; end
        if d == 0, args{end+1} = strtrim(buf); return, end %#ok<AGROW>
        if ch == ',' && d == 1, args{end+1} = strtrim(buf); buf = ''; else, buf(end+1) = ch; end %#ok<AGROW>
    end
    args{end+1} = strtrim(buf);
end
