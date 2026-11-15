% checkPencilScale_v001  Per-trial traced extent of the Pilot and Fraser
% ellipse trials in raw CSV units, and the physical size it implies under
% the two unit hypotheses. Companion to checkPencilFraction_v001.
%
% Reads raw x,y from every -enhanced.csv (no importer scaling, negation or
% de-duplication), keeps Ellipse trials (tracingScene_object1_imageNumber
% == 3 in the paired tracingSection.csv) with at least MIN_SAMP samples,
% and records max-min extent in x and y per trial.
%
% Hypotheses (iPad Pro 13" M4, 2752 x 2064 px, 264.615 ppi, UIKit scale 2):
%   points : 264.615/25.4/2 = 5.208965 units/mm
%   pixels : 9.73 (importDB_pilot_v001) or 10.41793 (importDB_fraser_v001)
% Reference: Dagmar's ruler measurement of the on-screen ellipse stimulus's
% major axis, ~120 mm (SESSION_LOG Session 72, 2026-07-24).
%
% OUTPUT
%   results/checkPencilScale_v001.csv          one row per trial
%   results/checkPencilScale_v001_summary.csv  one row per dataset
%
% Scripts run from src/; paths resolve from this file's location.

clearvars
srcDir  = fileparts(mfilename("fullpath"));
rootDir = fileparts(srcDir);
resDir  = fullfile(rootDir, "results");

PT_PER_MM   = 264.615 / 25.4 / 2;
PX_PER_MM   = struct("pilot", 9.73, "fraser", 10.41793);
ELLIPSE     = 3;
MIN_SAMP    = 100;
datasets    = ["pilot", "fraser"];

rows = {};
for ds = datasets
    files = dir(fullfile(rootDir, "data", ds, "automatic", "ID*", "*tracingSection-enhanced.csv"));
    if isempty(files)
        error("checkPencilScale:noFiles", "%s", "No enhanced CSVs found for " + ds);
    end
    for k = 1:numel(files)
        fe = fullfile(files(k).folder, files(k).name);
        fm = strrep(fe, "-enhanced.csv", ".csv");
        if ~isfile(fm)
            error("checkPencilScale:noMeta", "%s", "Missing metadata file " + fm);
        end
        m = readtable(fm);
        t = readtable(fe);
        [~, subj] = fileparts(files(k).folder);
        for r = m.trial(m.tracingScene_object1_imageNumber == ELLIPSE)'
            sel = t.response_index == r;
            if nnz(sel) < MIN_SAMP
                continue
            end
            x = double(t.x(sel)); y = double(t.y(sel));
            rows(end+1, :) = {ds, string(subj), string(files(k).name), r, nnz(sel), ...
                max(x)-min(x), max(y)-min(y)}; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, "VariableNames", ["dataset", "subject", "file", "trial", ...
    "nSamples", "extentX_units", "extentY_units"]);
T.mmX_ifPoints = T.extentX_units / PT_PER_MM;
T.mmY_ifPoints = T.extentY_units / PT_PER_MM;
pxmm = arrayfun(@(d) PX_PER_MM.(d), T.dataset);
T.mmX_ifPixels = T.extentX_units ./ pxmm;
T.mmY_ifPixels = T.extentY_units ./ pxmm;
writetable(T, fullfile(resDir, "checkPencilScale_v001.csv"));

S = groupsummary(T, "dataset", ["median", "min", "max"], ...
    ["extentX_units", "extentY_units", "mmX_ifPoints", "mmX_ifPixels"]);
nSubj = groupsummary(T, ["dataset", "subject"]);
nSubj = groupsummary(nSubj, "dataset");
S.nSubjects = nSubj.GroupCount;
disp(S)
writetable(S, fullfile(resDir, "checkPencilScale_v001_summary.csv"));
