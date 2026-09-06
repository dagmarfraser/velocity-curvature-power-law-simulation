% AUDITSEMLOOKUPANDQUADRANTS_V001  Reconcile the two SEM lookups in play, and
% produce the corrected precision-vs-identifiability cross-tabulation.
%
% WHY THIS EXISTS
%   Session 103 (2026-09-02) found that the paper reports SEM from one lookup
%   and Fig 0 plots it from another, and that they disagree. The corrected
%   quadrant counts were computed in a console session with no named producer.
%   This script is that producer, so the numbers in the Fig 0 caption, in
%   claude.md's submission-critical row, and in
%   docs/CoherencePass_Skeleton_v002_v001.md item 1 all have one traceable
%   origin and can be regenerated after any mat is rebuilt.
%
% THE TWO LOOKUPS
%   CANONICAL (what the paper's own published SEM values use, per
%   knockDownFlags_v001.m lines 84-88 and checkFraserSEMCentroid_v001.m lines
%   64-69): snap to the nearest (alpha, sigma, fs) grid node, then take the
%   MEAN of sem over ALL betaGen and ALL VGF rows at that node, per pipeline.
%
%   SCHEMATIC (what plotPaperRoadmapSchematic_v001.m line 86 does): snap the
%   same way, additionally filter to the single betaGen node nearest 1/3, then
%   take row.sem(1). Its REGISTRY has no VGF field, so 14 rows survive that
%   filter and the first is taken blind -- always VGF = 90.017, the grid's
%   lowest node, for every dataset.
%
%   The schematic lookup is reported here for comparison only. It should not
%   be used; this script exists partly to quantify how far wrong it goes.
%
% OUTPUTS (console only, no files written)
%   1. Reproduction check of the paper's seven published SG-IRLS SEM values.
%   2. Canonical SEM, all 7 datasets x 6 pipelines, with adequacy flags.
%   3. Schematic-lookup SEM for the same cells, and the VGF spread that makes
%      the difference.
%   4. Identifiability verdicts recomputed from the v012 corpora, using
%      compare42CellVerdict_v001.m's own rule.
%   5. The precision x identifiability cross-tabulation.
%   6. The corpus-wide adequate-coordinate fraction, with its threshold and
%      denominator stated explicitly.
%
% Fraser, D.S. (2026)  v001

clearvars

%% CONFIG -----------------------------------------------------------------
CFG.srcDir      = fileparts(mfilename('fullpath'));
CFG.semMat      = 'perCoordinateSEM_v2_001.mat';
CFG.loopPattern = 'loopClosureResults_%s_all_shaped_xu_v012.mat';
CFG.pipelines   = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
CFG.semAdequate = 0.0108;      % MDC/2.77, the paper's stated criterion
CFG.semRounded  = 0.011;       % the rounded variant, see the note in step 6
CFG.passCov     = 0.95;        % compare42CellVerdict_v001.m's own bands
CFG.condCov     = 0.90;
CFG.fullGridN   = 18045720/5;  % configurations / replications = coordinates

% Empirical centroids. Regenerate with auditNoiseCoordinates_v001.m rather
% than editing by hand; PUBLISHED is the paper's SG-IRLS SEM list, used here
% only as a reproduction target.
CFG.reg = struct( ...
  'name',     {'Fraser','Zarandi','Cook CTRL','Cook ASD','Hickman PLAC','Hickman HALO','Dhieb'}, ...
  'tag',      {'Fraser','Zarandi','Cook_CTRL','Cook_ASD','Hickman_PLAC','Hickman_HALO','Dhieb'}, ...
  'alpha',    {4.289,   3.184,     4.771,      5.062,     5.337,          5.424,         2.497}, ...
  'sigma',    {2.009,   4.77,      8.15,       7.84,      7.17,           7.42,          7.50}, ...
  'fs',       {240,     100,       133,        133,       133,            133,           100}, ...
  'published',{0.0032,  0.0079,    0.0051,     0.0046,    0.0040,         0.0040,        0.0189});

if isempty(CFG.srcDir), CFG.srcDir = pwd; end
cd(CFG.srcDir);

fprintf('=== AUDIT: SEM lookup reconciliation + quadrant cross-tab ===\n');
fprintf('    run %s\n\n', datetime('now','Format','yyyy-MM-dd HH:mm'));

if ~isfile(CFG.semMat)
    error('auditSEMLookup:noSEMFile', '%s', sprintf( ...
        'FAILED PATH: %s not found in %s.', CFG.semMat, CFG.srcDir));
end
S = load(CFG.semMat, 'coordTable');
T = S.coordTable;

gA = sort(unique(T.alpha));
gS = sort(unique(T.sigma));
gF = sort(unique(T.fs));
gB = sort(unique(T.betaGen));
[~, biNear] = min(abs(gB - 1/3));

nDS = numel(CFG.reg);
nPP = numel(CFG.pipelines);

semCanon  = nan(nDS, nPP);
semSchem  = nan(nDS, nPP);
vgfSpread = nan(nDS, 3);      % [min max ratio] for SG-IRLS at the node
snapped   = nan(nDS, 3);

%% 1-3: both lookups ------------------------------------------------------
for d = 1:nDS
    [~,ai] = min(abs(gA - CFG.reg(d).alpha));
    [~,si] = min(abs(gS - CFG.reg(d).sigma));
    [~,fi] = min(abs(gF - CFG.reg(d).fs));
    snapped(d,:) = [gA(ai) gS(si) gF(fi)];

    atNode = T.alpha==gA(ai) & T.sigma==gS(si) & T.fs==gF(fi);
    if ~any(atNode)
        error('auditSEMLookup:noRows', '%s', sprintf( ...
            'No coordTable rows for %s at snapped (alpha=%.3f, sigma=%.2f, fs=%d).', ...
            CFG.reg(d).name, gA(ai), gS(si), gF(fi)));
    end

    for p = 1:nPP
        m = atNode & T.pipeline == CFG.pipelines(p);
        semCanon(d,p) = mean(T.sem(m), 'omitnan');

        mb = m & T.betaGen == gB(biNear);
        sub = T(mb, :);
        if ~isempty(sub)
            semSchem(d,p) = sub.sem(1);       % reproduces the schematic's blind index
        end
        if CFG.pipelines(p) == "SG-IRLS"
            % Spread across VGF ONLY, i.e. at the fixed betaGen node. That is the
            % quantity the schematic bug ranges over: it pins betaGen first, so VGF
            % is the only axis left varying when row.sem(1) is taken blind. Spread
            % over betaGen AND VGF together is far larger (130x to 600x) but is not
            % what that bug is exposed to.
            v = T.sem(mb);
            vgfSpread(d,:) = [min(v) max(v) max(v)/min(v)];
        end
    end
end

fprintf('--- 1. Reproduction of the paper''s published SG-IRLS SEM list ---\n');
fprintf('%-14s %10s %10s %8s\n','dataset','published','canonical','match');
allMatch = true;
for d = 1:nDS
    p6 = semCanon(d, CFG.pipelines == "SG-IRLS");
    okM = abs(p6 - CFG.reg(d).published) < 5e-4;
    allMatch = allMatch && okM;
    fprintf('%-14s %10.4f %10.5f %8s\n', CFG.reg(d).name, CFG.reg(d).published, p6, string(okM));
end
if ~allMatch
    warning('auditSEMLookup:reproFailed', '%s', ...
        ['At least one published SEM value did not reproduce under the ' ...
         'canonical lookup. Either the mat has changed or the lookup has. ' ...
         'Resolve before trusting anything below.']);
end

fprintf('\n--- 2. Canonical SEM, all cells (adequate = < %.4f) ---\n', CFG.semAdequate);
printGrid(semCanon, CFG, '%9.5f');
adequate = semCanon < CFG.semAdequate;
fprintf('    adequate cells: %d / %d\n', sum(adequate(:)), numel(adequate));

fprintf('\n--- 3. Schematic lookup (row.sem(1)) for the same cells ---\n');
printGrid(semSchem, CFG, '%9.5f');
fprintf('    adequate under schematic lookup: %d / %d\n', ...
    sum(semSchem(:) < CFG.semAdequate), numel(semSchem));
fprintf('    SG-IRLS sem spread across VGF at the fixed betaGen node (min, max, ratio):\n');
for d = 1:nDS
    fprintf('      %-14s %8.5f %8.5f %7.1fx\n', CFG.reg(d).name, vgfSpread(d,:));
end

%% 4: identifiability verdicts -------------------------------------------
fprintf('\n--- 4. Identifiability verdicts recomputed from v012 corpora ---\n');
verdict = strings(nDS, nPP);
covPct  = nan(nDS, nPP);
for d = 1:nDS
    f = fullfile(CFG.srcDir, sprintf(CFG.loopPattern, CFG.reg(d).tag));
    if ~isfile(f)
        error('auditSEMLookup:noLoopMat', '%s', sprintf('FAILED PATH: %s not found.', f));
    end
    L = load(f, 'results');
    for p = 1:nPP
        st = arrayfun(@(r) r.invertStatus(p), L.results);
        ok = st ~= "no_beta_obs";
        if ~any(ok), verdict(d,p) = "no-data"; continue, end
        c = sum(st(ok) == "rise") / sum(ok);
        covPct(d,p) = 100*c;
        if     c >= CFG.passCov, verdict(d,p) = "PASS";
        elseif c >= CFG.condCov, verdict(d,p) = "CONDITIONAL";
        else,                    verdict(d,p) = "FAIL";
        end
    end
end
printGrid(covPct, CFG, '%9.2f');
nP = sum(verdict(:)=="PASS"); nC = sum(verdict(:)=="CONDITIONAL"); nF = sum(verdict(:)=="FAIL");
fprintf('    PASS=%d  CONDITIONAL=%d  FAIL=%d   (Finding #160: 4 / 6 / 32)\n', nP, nC, nF);
if ~isequal([nP nC nF], [4 6 32])
    warning('auditSEMLookup:verdictMismatch', '%s', ...
        ['Recomputed verdicts differ from Finding #160. Check mat versions ' ...
         'and the coverage bands before using any figure built on this.']);
end

%% 5: cross-tabulation ----------------------------------------------------
fprintf('\n--- 5. Precision x identifiability cross-tabulation ---\n');
fprintf('%-14s %12s %12s\n','verdict','adequate','inadequate');
for v = ["PASS","CONDITIONAL","FAIL"]
    fprintf('%-14s %12d %12d\n', v, sum(adequate(:) & verdict(:)==v), ...
                                    sum(~adequate(:) & verdict(:)==v));
end
fprintf('\n    adequate AND PASS      : %d\n', sum(adequate(:) & verdict(:)=="PASS"));
fprintf('    adequate but NOT PASS  : %d\n', sum(adequate(:) & verdict(:)~="PASS"));
fprintf('    (Fig 0''s caption currently claims the first of these is 0.)\n');

%% 6: corpus-wide adequate fraction --------------------------------------
fprintf('\n--- 6. Corpus-wide adequate-coordinate fraction ---\n');
nRows = height(T);
fprintf('    coordTable rows                  : %d\n', nRows);
fprintf('    full grid coordinates            : %d  (18,045,720 configs / 5 reps)\n', CFG.fullGridN);
fprintf('    coordinates with no SEM row      : %d (%.2f%%)\n', ...
    CFG.fullGridN - nRows, 100*(CFG.fullGridN - nRows)/CFG.fullGridN);
for thr = [CFG.semAdequate CFG.semRounded]
    n = sum(T.sem < thr);
    fprintf('    sem < %.4f : %d rows -> %.2f%% of coordTable, %.2f%% of full grid\n', ...
        thr, n, 100*n/nRows, 100*n/CFG.fullGridN);
end
fprintf(['    NOTE: the paper''s "79.8%%" corresponds to the 0.011 threshold over the\n' ...
         '    FULL-GRID denominator, i.e. counting coordinates that produced no SEM as\n' ...
         '    inadequate. That is conservative and defensible, but the threshold used\n' ...
         '    is 0.011 rather than the 0.0108 stated as the criterion elsewhere.\n']);

fprintf('\nDone. No files written.\n');

%% ------------------------------------------------------------------------
function printGrid(M, CFG, fmt)
    fprintf('%-14s', '');
    fprintf('%11s', CFG.pipelines);
    fprintf('\n');
    for d = 1:size(M,1)
        fprintf('%-14s', CFG.reg(d).name);
        for p = 1:size(M,2)
            fprintf(['  ' fmt], M(d,p));
        end
        fprintf('\n');
    end
end
