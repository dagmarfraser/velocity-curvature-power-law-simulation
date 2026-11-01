% tempoDAG_v001.m  Tempo DAG: nodes, evidence-tagged edges, identification tables, figure.
% Run from src/. Needs dagDSeparated_v001.m and dagAdjustmentSets_v001.m (same folder).
% Spec: docs/SPEC_TempoDAG_v001.md. Writes ../figures/tempoDAG_v001.png.
% Scope guard (claude.md, Retired threads): noise is a POSITION-level residual only.
% No node here carries a velocity-level or acceleration-level spectral claim.

%% CONFIG
MEASURED_U = false;    % true: subject traits U treated as measured (matched or covariate)
DRAW_FIG   = true;

here   = fileparts(mfilename('fullpath'));
figDir = fullfile(here, '..', 'figures');
figPng = fullfile(figDir, 'tempoDAG_v001.png');

%% Nodes: name, label (tex), observed, x, y
nodeDefs = {
    'U',   'subject traits',        false,  0,  7
    'G',   'group or drug',         true,   0,  5
    'Pdm', 'paradigm',              true,   0,  2.5
    'Geo', 'path geometry',         true,   0,  0.5
    'Dev', 'device, f_s',           true,   0, -1.5
    'T',   'tempo f_0',             false,  3,  5.5
    'Bg',  '\beta_{gen}',           false,  3,  3.5
    'N',   'noise \alpha, \sigma',  false,  4.5, -0.5
    'X',   'noise-free movement',   false,  5.5, 4.5
    'E',   'noise realisation',     false,  6.6,-0.5
    'R',   'recorded trajectory',   true,   8,  2
    'P',   'pipeline, cutoff',      true,   10.5, 4.5
    'Bo',  '\beta_{obs}',           true,   10.5, 3
    'Th',  'f_0 estimate',          true,   10.5, 1
    'Nh',  '\alpha, \sigma estimate', true, 10.5,-0.5
    'Bs',  '\beta_{gen}^*',         true,   13,  1.5
    'S',   'passes gate',           true,   13, -0.5 };
if MEASURED_U, nodeDefs{1, 3} = true; end
nName = string(nodeDefs(:, 1));
nLab  = string(nodeDefs(:, 2));
nObs  = cell2mat(nodeDefs(:, 3));
nXY   = cell2mat(nodeDefs(:, 4:5));

%% Edges: from, to, class, note.  Classes: def, sim, obs, lit, assumed, contested
edgeDefs = {
    'U',  'G',  'assumed',   'recruitment, matching'
    'U',  'T',  'assumed',   'age, skill, practice'
    'U',  'Bg', 'assumed',   ''
    'U',  'N',  'assumed',   ''
    'G',  'T',  'obs',       'Cook ASD faster: p = .045 on subject means, p = .13 on trials'
    'G',  'Bg', 'contested', 'the scientific question'
    'Pdm','T',  'assumed',   'pacing, instructions'
    'Pdm','N',  'assumed',   'hypothesis (metronome whitening), confounded with device and geometry in Zarandi; position level only'
    'Pdm','Geo','assumed',   'task template'
    'Geo','T',  'assumed',   'f0 = VGF / L_beta'
    'Geo','Bg', 'contested', 'minimum-jerk predicts shape dependence (Huh and Sejnowski, 2015)'
    'Geo','X',  'def',       ''
    'Dev','N',  'obs',       'platform noise differs; estimated alpha rises as fs falls 240 to 60 Hz (Finding #22); Dhieb pointer acceleration'
    'Dev','R',  'def',       'sampling rate, internal smoothing'
    'T',  'X',  'def',       ''
    'Bg', 'X',  'def',       ''
    'T',  'Bg', 'contested', 'genuine behaviour, or measurement only?'
    'T',  'N',  'obs',       'within-dataset f0 to alpha: Cook r +.75 to +.81, Zarandi +.39, Dhieb flat (April analyses); may be estimator, not biology (T2b)'
    'N',  'E',  'def',       ''
    'E',  'R',  'def',       ''
    'X',  'R',  'def',       ''
    'R',  'Bo', 'sim',       'carries the filter pull and the noise pull (#233, #226, #234, #249)'
    'P',  'Bo', 'def',       ''
    'R',  'Th', 'def',       'f0 = FFT peak (0.1-5 Hz) of x and y of the SAME recorded position (runner estimateF0_local); spec T2b'
    'R',  'Nh', 'def',       'template-subtracted IRASA'
    'Th', 'Nh', 'def',       'template is subtracted at the ESTIMATED f0 (runner estimateF0_local), so f0 error enters the noise estimate'
    'Bo', 'Bs', 'def',       ''
    'Th', 'Bs', 'def',       ''
    'Nh', 'Bs', 'def',       ''
    'Bo', 'S',  'def',       'gate: map locally monotonic at (Th, Nh) around Bo'
    'Th', 'S',  'def',       ''
    'Nh', 'S',  'def',       '' };
eFrom = string(edgeDefs(:, 1));  eTo = string(edgeDefs(:, 2));
eCls  = string(edgeDefs(:, 3));  eNote = string(edgeDefs(:, 4));

G = digraph(eFrom, eTo);
if ~isdag(G), error('tempoDAG_v001:cyclic', 'Edge list contains a cycle.'); end
[ok, loc] = ismember(string(G.Nodes.Name), nName);
if ~all(ok), error('tempoDAG_v001:nodes', 'Edge list uses a node missing from nodeDefs.'); end
if any(~ismember(nName, string(G.Nodes.Name))), error('tempoDAG_v001:orphan', 'Node in nodeDefs has no edge.'); end
obsNodes = nName(nObs);

%% 1. Edge table
fprintf('\n== 1. Edges (%d) ==\n', height(G.Edges));
for e = 1:numel(eFrom)
    fprintf('%-4s -> %-4s  %-9s %s\n', eFrom(e), eTo(e), eCls(e), eNote(e));
end

%% 2. Directed routes from true tempo T to the inverted estimate Bs
fprintf('\n== 2. Directed routes T -> Bs ==\n');
P = allpaths(G, 'T', 'Bs');             % named-node graph: paths come back as name cells
for k = 1:numel(P)
    nm = string(P{k});
    nm = nm(:)';
    if any(nm == "Bg"), tag = 'behavioural (through beta_gen)';
    elseif any(nm == "Th"), tag = 'estimator (through f0 estimate)';
    else, tag = 'measurement (through the record to beta_obs)'; end
    fprintf('%-34s %s\n', strjoin(nm, ' -> '), tag);
end
fprintf('Inversion works by cancelling the measurement and estimator routes. That is a\n');
fprintf('modelling claim (spec T2), not a d-separation statement.\n');

%% 3. Estimands and minimal backdoor adjustment sets among observed nodes
est = {
    'G', 'Bo', 'raw group contrast of beta_obs'
    'G', 'Bg', 'group effect on the latent exponent (target)'
    'T', 'Bg', 'behavioural tempo effect (the contested edge)' };
fprintf('\n== 3. Backdoor adjustment sets (observed: %s) ==\n', strjoin(obsNodes, ', '));
for k = 1:size(est, 1)
    S = dagAdjustmentSets_v001(G, est{k, 1}, est{k, 2}, obsNodes);
    fprintf('%s -> %s  (%s)\n', est{k, 1}, est{k, 2}, est{k, 3});
    if isempty(S)
        fprintf('    NOT IDENTIFIED from observed nodes alone\n');
    else
        for q = 1:numel(S), fprintf('    adjust for {%s}\n', strjoin(S{q}, ', ')); end
    end
end

%% 4. Does adjusting for the ESTIMATED tempo block the tempo route?
% Keep only the tempo route: G -> T -> X -> R -> Bo, with R -> Th and R -> Nh.
Gt = rmnode(G, {'U','Pdm','Geo','Dev','Bg','N','E','P','Bs','S'});
fprintf('\n== 4. Tempo route G -> T -> X -> R -> Bo ==\n');
fprintf('G _||_ Bo | {T}      : %d  (true tempo blocks it)\n',  dagDSeparated_v001(Gt, 'G', 'Bo', "T"));
fprintf('G _||_ Bo | {Th}     : %d  (estimated tempo does not: R is a parent of both Th and Bo)\n', ...
    dagDSeparated_v001(Gt, 'G', 'Bo', "Th"));
fprintf('G _||_ Bo | {Th, Nh} : %d\n', dagDSeparated_v001(Gt, 'G', 'Bo', ["Th"; "Nh"]));

%% 5. X is a collider on T and Bg: selection on any descendant opens T - Bg
% Remove the contested T -> Bg and Geo -> Bg edges, then condition within group and subject.
c = eCls == "contested" & ~(eFrom == "G");
G0 = rmedge(G, findedge(G, cellstr(eFrom(c)), cellstr(eTo(c))));
base = ["U"; "G"];
fprintf('\n== 5. With no behavioural T -> Bg edge, conditioning within group and subject ==\n');
fprintf('T _||_ Bg | {U, G}          : %d\n', dagDSeparated_v001(G0, 'T', 'Bg', base));
for z = ["S"; "Bo"; "Th"; "Bs"; "R"]'
    fprintf('T _||_ Bg | {U, G, %-2s}     : %d\n', z, dagDSeparated_v001(G0, 'T', 'Bg', [base; z]));
end
fprintf('0 = dependent: the gate (S) or any other descendant of X manufactures a T-Bg association.\n');

%% 6. Implied marginal independencies (checks that the graph says what we mean)
fprintf('\n== 6. Marginal independencies implied ==\n');
chk = {'Dev','Bg'; 'P','T'; 'P','Bg'; 'P','N'; 'P','Th'};
for k = 1:size(chk, 1)
    fprintf('%-3s _||_ %-3s : %d\n', chk{k, 1}, chk{k, 2}, dagDSeparated_v001(G, chk{k, 1}, chk{k, 2}));
end

%% 7. Figure
if DRAW_FIG
    if ~isfolder(figDir), error('tempoDAG_v001:figDir', 'Missing folder: %s', figDir); end
    xy = nXY(loc, :);
    fig = figure('Color', 'w', 'Position', [60 60 1800 820]);
    p = plot(G, 'XData', xy(:, 1), 'YData', xy(:, 2), 'NodeLabel', {}, 'ArrowSize', 11, ...
        'MarkerSize', 9, 'NodeColor', [.15 .15 .15], 'EdgeColor', [.45 .45 .45]);
    isObs = nObs(loc);
    highlight(p, find(~isObs), 'Marker', 'o', 'NodeColor', [.15 .15 .15]);
    highlight(p, find(isObs),  'Marker', 's', 'NodeColor', [.15 .15 .15]);
    hold on
    for e = 1:numel(eFrom)
        st = edgeStyle(eCls(e));
        highlight(p, find(string(G.Nodes.Name) == eFrom(e)), find(string(G.Nodes.Name) == eTo(e)), ...
            'EdgeColor', st.col, 'LineStyle', st.ls, 'LineWidth', st.lw);
    end
    for k = 1:numnodes(G)
        text(xy(k, 1), xy(k, 2) + 0.42, nLab(loc(k)), 'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'Interpreter', 'tex');
    end
    cls = ["def" "sim" "obs" "lit" "assumed" "contested"];
    lab = ["by construction" "simulation-demonstrated" "observed association" ...
           "literature, unverified" "assumed" "contested or unidentified"];
    hh = gobjects(1, numel(cls));
    for k = 1:numel(cls)
        st = edgeStyle(cls(k));
        hh(k) = plot(nan, nan, 'Color', st.col, 'LineStyle', st.ls, 'LineWidth', st.lw);
    end
    legend(hh, lab, 'Location', 'eastoutside', 'Box', 'off');
    title({'Tempo DAG: behaviour (left) and measurement chain (right)', ...
           'circle: latent    square: observed'}, 'FontWeight', 'normal');
    axis off
    exportgraphics(fig, figPng, 'Resolution', 300);
    fprintf('\nFigure: %s\n', figPng);
end

function st = edgeStyle(cls)
switch string(cls)
    case "def",       st = struct('col', [.45 .45 .45],  'ls', '-',  'lw', 1.0);
    case "sim",       st = struct('col', [0 .447 .698],  'ls', '-',  'lw', 2.4);
    case "obs",       st = struct('col', [0 .62 .451],   'ls', '-',  'lw', 1.6);
    case "lit",       st = struct('col', [.8 .475 .655], 'ls', '--', 'lw', 1.4);
    case "assumed",   st = struct('col', [.45 .45 .45],  'ls', '--', 'lw', 1.0);
    case "contested", st = struct('col', [.835 .369 0],  'ls', '--', 'lw', 2.4);
    otherwise
        error('tempoDAG_v001:class', 'Unknown edge class "%s".', cls);
end
end
