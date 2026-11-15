%% checkMedianShareV015Spec_v001  Deep Review GD13: family shares of the constellation median on the v015 corpus (Table 6b basis). Saves results/checkMedianShareV015Spec_v001.mat.
r = checkConstellationMedianPipelineContribution_v002();   %#ok<NASGU>  (prints; r is a table with olsShare, lmlsShare, irlsShare)
ROOT = fileparts(fileparts(mfilename("fullpath")));
out = struct("dataset", string(r.dataset), "N", r.N, "olsShare", r.olsShare, "lmlsShare", r.lmlsShare, "irlsShare", r.irlsShare);
save(fullfile(ROOT, "results", "checkMedianShareV015Spec_v001.mat"), "-struct", "out");
