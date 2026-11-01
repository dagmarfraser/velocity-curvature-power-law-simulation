function T = printDatasetNoiseMedians_v001()
% printDatasetNoiseMedians_v001  Per-trial IRASA alpha and template sigma (mm): dataset medians [IQR] on the v015 corpus (N_REPS = 200) and on the first-reported basis (N_REPS = 20: v007, Fraser v008).
% Closes the OPEN flag in the manuscript's empirical-values bullet (v015 medians quoted from PATCH_DesignRecs_TempoMatching_v002.md; Fraser alpha 4.38, sigma 1.71 mm).
% alpha = alphaIRA and sigma = sigmaMM of each trial's result record, finite values only. A missing file is reported (NaN row, printed), never skipped silently; the five core datasets must exist.
% Saves results/printDatasetNoiseMedians_v001.mat and .csv.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');
core = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];  ds = [core "Zarandi" "Dhieb"];
v0 = containers.Map(cellstr(ds), cellstr(["v008" "v007" "v007" "v007" "v007" "v007" "v007"]));
rows = cell(0, 12);
for d = ds
    for b = ["v015" "N20"]
        ver = b;  if b == "N20", ver = v0(char(d)); end
        f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', d, ver));
        if ~isfile(f)
            if ismember(d, core), error('noiseMed:missing', 'Missing core file %s', f); end
            fprintf('NOT FOUND (reported, not skipped): %s\n', f);  rows(end+1, :) = {d, b, ver, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            continue
        end
        R = load(f, 'results').results;
        a = arrayfun(@(r) double(r.alphaIRA), R(:));  s = arrayfun(@(r) double(r.sigmaMM), R(:));
        a = a(isfinite(a));  s = s(isfinite(s));
        rows(end+1, :) = {d, b, ver, numel(R), numel(a), median(a), prctile(a, 25), prctile(a, 75), numel(s), median(s), prctile(s, 25), prctile(s, 75)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'dataset','basis','file','nTrials','nAlpha','alphaMed','alphaQ25','alphaQ75','nSigma','sigmaMed','sigmaQ25','sigmaQ75'});
T.dataset = string(T.dataset);  T.basis = string(T.basis);  T.file = string(T.file);
disp(T(:, {'dataset','basis','file','nTrials','alphaMed','alphaQ25','alphaQ75','sigmaMed'}));
fprintf('sigma IQR [Q25 Q75] by row:\n');  for i = 1:height(T), if ~isnan(T.sigmaMed(i)), fprintf('  %-13s %-5s %.2f [%.2f, %.2f]\n', T.dataset(i), T.basis(i), T.sigmaMed(i), T.sigmaQ25(i), T.sigmaQ75(i)); end, end
save(fullfile(resDir, 'printDatasetNoiseMedians_v001.mat'), 'T');
end
