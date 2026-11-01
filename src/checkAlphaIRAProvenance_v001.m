function T = checkAlphaIRAProvenance_v001()
% checkAlphaIRAProvenance_v001  Why the stored per-trial alphaIRA differs between the N_REPS = 20 result files (v007; Fraser v008) and v015.
% Compares, per dataset and trial (matched by trialID), alphaIRA in each file against the CURRENT noise characterisation (noiseCharacterisation_<name>.mat, ira_alphaMean),
% and alphaMaj / alphaMin (the fields that drive the surrogate noise) between the two files. Finding (2026-10-29): v015 equals the characterisation exactly; the N_REPS = 20 alphaIRA is a stale
% field (the retired v003 characterisation, SESSION_LOG Session 49 'Pain point 2'); alphaMaj/alphaMin are identical, so no simulation result depends on it. A trial missing from the characterisation errors.
% Saves results/checkAlphaIRAProvenance_v001.mat.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');
map = {'Cook_CTRL','noiseCharacterisation_cook.mat','v007'; 'Cook_ASD','noiseCharacterisation_cookASD.mat','v007'; 'Hickman_PLAC','noiseCharacterisation_hickmanPLAC.mat','v007'; ...
       'Hickman_HALO','noiseCharacterisation_hickmanHALO.mat','v007'; 'Zarandi','noiseCharacterisation_zarandi.mat','v007'; 'Fraser','noiseCharacterisation_fraser.mat','v008'};
rows = cell(0, 8);
for i = 1:size(map, 1)
    d = map{i, 1};  B = load(fullfile(srcDir, map{i, 2})).bioResults;  if istable(B), B = table2struct(B); end;  cid = string({B.trialID}');  ira = [B.ira_alphaMean]';
    a15 = load(fullfile(srcDir, "loopClosureResults_" + d + "_all_shaped_xu_v015.mat"), 'results').results;  a07 = load(fullfile(srcDir, "loopClosureResults_" + d + "_all_shaped_xu_" + map{i, 3} + ".mat"), 'results').results;
    [t15, l15] = ismember(string({a15.trialID}'), cid);  [t07, l07] = ismember(string({a07.trialID}'), cid);
    if ~all(t15) || ~all(t07), error('alphaProv:trials', '%s: trials missing from the noise characterisation', d); end
    if ~isequal(string({a15.trialID}), string({a07.trialID})), error('alphaProv:ids', '%s: trial IDs differ between files', d); end
    A15 = [a15.alphaIRA]';  A07 = [a07.alphaIRA]';
    rows(end+1, :) = {string(d), numel(a15), max(abs(A15 - ira(l15))), max(abs(A07 - ira(l07))), median(A07, 'omitnan'), median(A15, 'omitnan'), max(abs([a07.alphaMaj] - [a15.alphaMaj])), max(abs([a07.alphaMin] - [a15.alphaMin]))}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'dataset','nTrials','maxDiff_v015_vs_char','maxDiff_N20_vs_char','median_N20','median_v015','maxDiff_alphaMaj','maxDiff_alphaMin'});
disp(T);  save(fullfile(resDir, 'checkAlphaIRAProvenance_v001.mat'), 'T');
end
