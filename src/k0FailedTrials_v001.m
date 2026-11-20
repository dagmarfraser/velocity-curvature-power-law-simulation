% k0FailedTrials_v001  Tabulate emStatus failure reasons in the K=0 (nH0 v033) result files
% and compare with K=4 / K=1 files for the same trialIDs where those are available locally.
% Read-only: loads .mat files, writes nothing. Created 2026-11-20.
% Local K=0 copies: BRM/outputs_K0. K=4/K=1 comparators: old tree PowerLawSimulationPreReg/src (not copied).
b  = '/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee';
k0 = fullfile(b,'PowerLaw_Repro_2026/BRM/outputs_K0');
old = fullfile(b,'PowerLawSimulationPreReg/src');
ds = {'Hickman_PLAC','Hickman_HALO'};
for i = 1:numel(ds)
    d = ds{i};
    r0 = load(fullfile(k0,sprintf('loopClosureResults_%s_all_endMatched_segdetect_nH0_v033.mat',d)),'results').results;
    id0 = string({r0.trialID}'); st0 = string(arrayfun(@(x) x.emStatus(1), r0))';
    fail0 = st0 ~= "ok";
    fprintf('\n##### %s K=0: %d trials, %d status ok, %d not ok\n', d, numel(r0), sum(~fail0), sum(fail0));
    % normalise: replace digit runs by N, group
    key = regexprep(st0, '\d+(\.\d+)?', 'N');
    [u,~,g] = unique(key(fail0)); c = accumarray(g,1);
    [c,o] = sort(c,'descend'); u = u(o);
    for j = 1:numel(u), fprintf('  %3d  %s\n', c(j), u(j)); end
    % numeric lengths for the "end-matched length N below minimum" class
    t = regexp(st0(fail0),'end-matched length (\d+)','tokens','once');
    t = t(~cellfun(@isempty,t)); if ~isempty(t), L = cellfun(@(x) str2double(x{1}),t);
        fprintf('  end-matched lengths: median %g, min %g, max %g, n=%d\n', median(L), min(L), max(L), numel(L)); end
    % comparators
    cmp = struct('K4',sprintf('loopClosureResults_%s_all_endMatched_segdetect_v028.mat',d), ...
                 'K1',sprintf('loopClosureResults_%s_all_endMatched_segdetect_nH1_v028.mat',d));
    for K = ["K4","K1"]
        f = fullfile(old,cmp.(K));
        if ~isfile(f), fprintf('  %s: file not available locally\n',K); continue; end
        rk = load(f,'results').results;
        idk = string({rk.trialID}'); stk = string(arrayfun(@(x) x.emStatus(1), rk))';
        [tf,loc] = ismember(id0, idk);
        failK = false(size(id0)); failK(tf) = stk(loc(tf)) ~= "ok";
        fprintf('  %s: %d trials, %d not ok; same trialIDs matched %d/%d\n', K, numel(rk), sum(stk~="ok"), sum(tf), numel(id0));
        fprintf('    K0 fail & %s fail: %d | K0 fail & %s ok: %d | K0 ok & %s fail: %d | trial only in K0: %d\n', ...
            K, sum(fail0&tf&failK), K, sum(fail0&tf&~failK), K, sum(~fail0&tf&failK), sum(~tf));
        kk = regexprep(stk,'\d+(\.\d+)?','N'); [u2,~,g2]=unique(kk(stk~="ok")); c2=accumarray(g2,1);
        for j=1:numel(u2), fprintf('      %s: %3d  %s\n',K,c2(j),u2(j)); end
        % reasons among K0 fails when also failing / being ok in comparator
        idx = find(fail0&tf&~failK);
        if ~isempty(idx), fprintf('    K0-only failures (ok in %s): reason counts: ', K);
            [u3,~,g3]=unique(key(idx)); c3=accumarray(g3,1); for j=1:numel(u3), fprintf('[%d] %s; ',c3(j),u3(j)); end; fprintf('\n'); end
    end
end
