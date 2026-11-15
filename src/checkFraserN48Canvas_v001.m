function T = checkFraserN48Canvas_v001()
%CHECKFRASERN48CANVAS_V001  Raw x,y of every Fraser N=48 participant inside the UIKit point canvas (2026-11-03).
%   PLAN_PencilScaleCorrection_v001 section 1.3: a file in pixels, or a different device, must stop the analysis.
%   checkPencilFraction_v001 (2026-11-02) showed the N=45 corpus fills the 1376 x 1032 pt canvas and never exceeds it;
%   ID401-403 arrived after it. This reads raw x,y from every -enhanced.csv in data/fraser/automatic (no importer
%   scaling, negation or de-duplication) and errors if any value is below 0 or above the canvas. Reports per-participant
%   maxima so a participant recorded on another device would show. Writes results/checkFraserN48Canvas_v001.csv.
%   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  rootDir = fileparts(srcDir);
CANVAS = [1376 1032];   % landscape point canvas, iPad Pro 13" M4 (2752 x 2064 px at UIKit scale 2)
files = dir(fullfile(rootDir, 'data', 'fraser', 'automatic', 'ID*', '*-enhanced.csv'));
if isempty(files), error('checkFraserN48Canvas:noFiles', '%s', 'No -enhanced.csv files under data/fraser/automatic.'); end
ids = strings(numel(files), 1);  mn = NaN(numel(files), 2);  mx = NaN(numel(files), 2);  n = zeros(numel(files), 1);
for k = 1:numel(files)
    o = detectImportOptions(fullfile(files(k).folder, files(k).name));
    o.SelectedVariableNames = {'x', 'y'};
    D = readtable(fullfile(files(k).folder, files(k).name), o);
    [~, ids(k)] = fileparts(files(k).folder);
    v = [D.x D.y];  v = v(all(isfinite(v), 2), :);
    n(k) = size(v, 1);
    if n(k) > 0, mn(k, :) = min(v, [], 1);  mx(k, :) = max(v, [], 1); end
end
[g, pid] = findgroups(ids);
T = table(pid, splitapply(@sum, n, g), splitapply(@(a) min(a), mn(:, 1), g), splitapply(@(a) max(a), mx(:, 1), g), ...
    splitapply(@(a) min(a), mn(:, 2), g), splitapply(@(a) max(a), mx(:, 2), g), ...
    'VariableNames', ["participant" "nSamples" "xMin" "xMax" "yMin" "yMax"]);
writetable(T, fullfile(rootDir, 'results', 'checkFraserN48Canvas_v001.csv'));
disp(T(ismember(T.participant, ["ID401" "ID402" "ID403"]), :))
fprintf('%d files, %d participants. Corpus x %.1f-%.1f, y %.1f-%.1f (canvas %d x %d pt).\n', numel(files), height(T), ...
    min(T.xMin), max(T.xMax), min(T.yMin), max(T.yMax), CANVAS);
bad = T(T.xMin < 0 | T.yMin < 0 | T.xMax > CANVAS(1) | T.yMax > CANVAS(2), :);
if height(bad) > 0, disp(bad); error('checkFraserN48Canvas:outside', '%d participant(s) have x,y outside the point canvas: not points on this device.', height(bad)); end
fprintf('=== All %d participants inside the %d x %d point canvas: units are points ===\n', height(T), CANVAS);
end
