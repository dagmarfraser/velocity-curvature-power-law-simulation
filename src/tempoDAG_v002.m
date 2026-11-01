%% tempoDAG_v002.m  The paper's causal diagram (TODO C2): two readings of tempo, one measurement layer.
% Simplifies tempoDAG_v001 (17 nodes; full graph kept for the supplement) to the 12 nodes the
% paper needs, and draws tempo's relation to the group or condition both ways:
%   (a) confound: subject traits U drive both group and tempo (a common cause);
%   (b) mediator: the group changes tempo (G -> T), which then acts on beta.
% Variant (a2), tempo setting the condition (T -> G, as in the confound simulation, Finding #258),
% is analysed but not drawn. Edges as in tempoDAG_v001 and SPEC_TempoDAG_v001 §7 (Dagmar,
% 2026-10-22): T -> Bg contested; session and subject layer implicit. Retired-thread guard: the
% noise node is the POSITION-level residual only; no velocity- or acceleration-level node.
% Prints adjustment sets, d-separation checks and the collider check that the paper quotes.
% Run from src/. Needs dagDSeparated_v001.m and dagAdjustmentSets_v002.m (v001 skipped the empty set
% when exactly one candidate was left; found on this script's first run, 2026-10-24).
% Writes figures/tempoDAG_v002.png, results/tempoDAG_v002.mat and a dated log in logs/.

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
OUT_PNG = fullfile(ROOT, "figures", "tempoDAG_v002.png");
OUT_MAT = fullfile(ROOT, "results", "tempoDAG_v002.mat");

%% Nodes: name, label, observed, x, y, label offset dx, dy, alignment
nodeDefs = {
    'U',  'subject traits',            false,  0,    6.0,   0,    0.55, 'center'
    'G',  'group or condition',        true,   0,    2.0,   0,   -0.60, 'center'
    'T',  'tempo f_0',                 false,  3.4,  4.8,  -0.30, 0.15, 'right'
    'Bg', '\beta_{gen}',               false,  3.4,  1.0,  -0.30, 0,    'right'
    'N',  'noise (position residual)', false,  6.8,  6.4,   0,    0.55, 'center'
    'X',  'movement',                  false,  6.8,  3.0,   0,   -0.60, 'center'
    'R',  'recorded trajectory',       true,   9.8,  3.0,   0,   -0.60, 'center'
    'Bo', '\beta_{obs}',               true,  12.8,  5.4,   0,    0.55, 'center'
    'Th', 'estimated f_0',             true,  12.8,  3.0,   0.30, 0.35, 'left'
    'Nh', 'estimated noise',           true,  12.8,  0.6,   0,   -0.60, 'center'
    'Bs', '\beta_{gen}^*',             true,  15.8,  4.3,   0,    0.55, 'center'
    'S',  'passes the gate',           true,  15.8,  1.8,   0,   -0.60, 'center' };
nName = string(nodeDefs(:, 1));  nLab = string(nodeDefs(:, 2));
nObs  = cell2mat(nodeDefs(:, 3));  nXY = cell2mat(nodeDefs(:, 4:5));
nOff  = cell2mat(nodeDefs(:, 6:7));  nAl = string(nodeDefs(:, 8));
obsNodes = nName(nObs);

%% Edges common to both readings: from, to, class
common = {
    'U','T','assumed';  'U','Bg','assumed';  'U','N','assumed'
    'G','Bg','contested';  'T','Bg','contested'
    'T','X','def';  'Bg','X','def';  'T','N','obs'
    'X','R','def';  'N','R','def';  'R','Bo','sim'
    'R','Th','def';  'R','Nh','def';  'Th','Nh','def'
    'Bo','Bs','def';  'Th','Bs','def';  'Nh','Bs','def'
    'Bo','S','def';  'Th','S','def';  'Nh','S','def' };
reading = struct( ...
    "a",  {[common; {'U','G','assumed'}]}, ...     % confound: common cause
    "a2", {[common; {'T','G','assumed'}]}, ...     % confound: tempo sets the condition (not drawn)
    "b",  {[common; {'G','T','obs'}]});            % mediator: group changes tempo
RLAB = struct("a", "(a) tempo as a confound", "a2", "(a2) tempo sets the condition", "b", "(b) tempo as a mediator");

%% Setup
logFile = fullfile(ROOT, "logs", "tempoDAG_v002_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("dag2:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "tempoDAG_v002 run %s\nObserved nodes: %s\n", string(datetime("now")), strjoin(obsNodes, ", "));
res = struct();

for rk = ["a" "a2" "b"]
    E = reading.(rk);  Gr = digraph(string(E(:, 1)), string(E(:, 2)));
    if ~isdag(Gr), error("dag2:cyclic", "%s", "Reading " + rk + " has a cycle"); end
    if ~all(ismember(string(Gr.Nodes.Name), nName)) || numnodes(Gr) ~= numel(nName), error("dag2:nodes", "%s", "Node set mismatch in " + rk); end
    say_local(LG, "\n== %s: %d edges ==\n", RLAB.(rk), numedges(Gr));
    % 1. Backdoor adjustment sets: observed nodes only, then with true tempo, then with subject traits too
    for est = ["Bo" "Bg"]
        for extra = {strings(0, 1), "T", ["T"; "U"]}
            cand = [obsNodes; extra{1}];  cand = cand(cand ~= "G");
            S = dagAdjustmentSets_v002(Gr, "G", est, cand);
            tag = "observed"; if ~isempty(extra{1}), tag = "observed + " + strjoin(extra{1}, " + "); end
            say_local(LG, "G -> %s, candidates %s: %s\n", est, tag, fmtSets_local(S));
            res.(rk).("adj_" + est + "_" + strrep(strrep(tag, " + ", "_"), "observed", "obs")) = fmtSets_local(S);
        end
    end
    % 2. The tempo route alone: does true or estimated tempo block it? (as tempoDAG_v001 section 4)
    if rk == "b"
        Gt = rmnode(Gr, {'U', 'Bg', 'N', 'Bs', 'S'});
        d = [dagDSeparated_v001(Gt, "G", "Bo", "T"), dagDSeparated_v001(Gt, "G", "Bo", "Th"), dagDSeparated_v001(Gt, "G", "Bo", ["Th"; "Nh"])];
        say_local(LG, "Tempo route G -> T -> X -> R -> Bo: blocked by {T} %d, by {estimated f0} %d, by {estimated f0, estimated noise} %d\n", d);
        res.b.routeBlocked = d;
        % 3. Mediator-outcome confounding: with no direct G -> Bg, does conditioning on true tempo open G - Bg through U?
        Gn = rmedge(Gr, "G", "Bg");
        m = [dagDSeparated_v001(Gn, "G", "Bg"), dagDSeparated_v001(Gn, "G", "Bg", "T"), dagDSeparated_v001(Gn, "G", "Bg", ["T"; "U"])];
        say_local(LG, "No direct G -> Bg: G _||_ Bg marginally %d, given {T} %d, given {T, U} %d (0 = dependent)\n", m);
        res.b.mediatorConfounding = m;
    end
    % 4. Collider: drop the contested T -> Bg; does selecting on the gate or conditioning on a record descendant open T - Bg?
    G0 = rmedge(Gr, "T", "Bg");  base = ["U"; "G"];
    z = ["S" "Bo" "Th" "Bs" "R"];  c = arrayfun(@(q) dagDSeparated_v001(G0, "T", "Bg", [base; q]), z);
    c0 = dagDSeparated_v001(G0, "T", "Bg", base);
    say_local(LG, "No T -> Bg edge: T _||_ Bg | {U, G} %d; adding %s gives %s (0 = dependent)\n", c0, strjoin(z, ", "), strjoin(string(c), ", "));
    res.(rk).collider = [c0 c];
end

%% Anchor (fail loud): in the mediator reading the total effect needs no adjustment, with or without U offered
if res.b.adj_Bo_obs_T_U ~= "{} (no adjustment needed)" || res.b.adj_Bg_obs_T_U ~= "{} (no adjustment needed)"
    error("dag2:anchor", "%s", "Mediator reading: the empty set should suffice for G -> Bo and G -> Bg");
end
say_local(LG, "\nAnchor passed: the mediator reading needs no adjustment for the total effect.\n");

%% Figure: readings (a) and (b), same layout
fg = figure("Color", "w", "Position", [60 60 1500 1350]);  tl = tiledlayout(2, 1, "TileSpacing", "compact");
for rk = ["a" "b"]
    E = reading.(rk);  Gr = digraph(string(E(:, 1)), string(E(:, 2)));
    [~, loc] = ismember(string(Gr.Nodes.Name), nName);  xy = nXY(loc, :);
    nexttile;  hold on
    p = plot(Gr, "XData", xy(:, 1), "YData", xy(:, 2), "NodeLabel", {}, "ArrowSize", 10, "MarkerSize", 8, ...
        "NodeColor", [.15 .15 .15], "EdgeColor", [.45 .45 .45]);
    highlight(p, find(~nObs(loc)), "Marker", "o");  highlight(p, find(nObs(loc)), "Marker", "s");
    for e = 1:size(E, 1)
        st = edgeStyle_local(E{e, 3});
        highlight(p, findnode(Gr, E{e, 1}), findnode(Gr, E{e, 2}), "EdgeColor", st.col, "LineStyle", st.ls, "LineWidth", st.lw);
    end
    for k = 1:numnodes(Gr)
        text(xy(k, 1) + nOff(loc(k), 1), xy(k, 2) + nOff(loc(k), 2), nLab(loc(k)), "HorizontalAlignment", nAl(loc(k)), "FontSize", 11, "Interpreter", "tex");
    end
    title(RLAB.(rk), "FontWeight", "normal", "FontSize", 12);  axis off;  xlim([-1.8 17.4]);  ylim([-0.3 7.2]);
end
cls = ["def" "sim" "obs" "assumed" "contested"];
lab = ["by construction" "demonstrated in simulation" "observed association" "assumed" "contested"];
hh = gobjects(1, numel(cls));
for k = 1:numel(cls), st = edgeStyle_local(cls(k)); hh(k) = plot(nan, nan, "Color", st.col, "LineStyle", st.ls, "LineWidth", st.lw); end
lg = legend(hh, lab, "Box", "off", "Orientation", "horizontal");  lg.Layout.Tile = "south";
title(tl, "Tempo in two readings, one measurement layer (circle: latent; square: observed)", "FontWeight", "normal");
exportgraphics(fg, OUT_PNG, "Resolution", 300);
save(OUT_MAT, "res", "reading", "nodeDefs");
say_local(LG, "\nFigure: %s\nSaved: %s\nLog: %s\n", OUT_PNG, OUT_MAT, logFile);  fclose(LG);

%% =========================================================================
function s = fmtSets_local(S)
    if isempty(S), s = "NOT IDENTIFIED"; return, end
    parts = strings(1, numel(S));
    for q = 1:numel(S)
        if isempty(S{q}), parts(q) = "{} (no adjustment needed)"; else, parts(q) = "{" + strjoin(S{q}, ", ") + "}"; end
    end
    s = strjoin(parts, " or ");
end

function st = edgeStyle_local(cls)
    switch string(cls)
        case "def",       st = struct("col", [.45 .45 .45],  "ls", "-",  "lw", 1.0);
        case "sim",       st = struct("col", [0 .447 .698],  "ls", "-",  "lw", 2.4);
        case "obs",       st = struct("col", [0 .62 .451],   "ls", "-",  "lw", 1.6);
        case "assumed",   st = struct("col", [.45 .45 .45],  "ls", "--", "lw", 1.0);
        case "contested", st = struct("col", [.835 .369 0],  "ls", "--", "lw", 2.4);
        otherwise, error("dag2:class", "%s", "Unknown edge class " + string(cls));
    end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
