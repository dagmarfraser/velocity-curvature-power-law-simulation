% acc33Compare_v001  Content comparison of acc33 files: RDS copies (scp'd to scratch) vs old-tree src copies.
% Read-only. Created 2026-11-20. Local RDS copies were placed in the session scratchpad (not in Dropbox).
b   = '/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src';
rd  = '/private/tmp/claude-502/-Users-d-s-fraser-Library-CloudStorage-Dropbox/4f862f87-deac-4410-bda4-88123502e3a7/scratchpad/acc';
nm = {'loopClosureResults_Cook_CTRL_all_endMatched_segdetect_acc33_v028.mat', ...
      'loopClosureResults_Cook_CTRL_all_endMatched_segdetect_acc33_v033.mat', ...
      'loopClosureResults_Cook_CTRL_all_endMatched_segdetect_nH0_acc33_v033.mat'};
for i = 1:numel(nm)
    A = load(fullfile(b,nm{i})); B = load(fullfile(rd,['rds_' nm{i}]));
    fprintf('\n== %s\n', nm{i});
    fprintf('  isequaln(all variables) = %d\n', isequaln(A,B));
    fn = union(fieldnames(A),fieldnames(B));
    for j = 1:numel(fn)
        if ~isfield(A,fn{j}) || ~isfield(B,fn{j}), fprintf('  %s: only in one file\n',fn{j}); continue; end
        if ~isequaln(A.(fn{j}),B.(fn{j})), fprintf('  DIFF in %s\n',fn{j}); end
    end
    if isfield(A,'results') && isfield(B,'results') && ~isequaln(A.results,B.results)
        fa = fieldnames(A.results); n = numel(A.results);
        for j = 1:numel(fa)
            nd = 0; for t = 1:n, if ~isequaln(A.results(t).(fa{j}),B.results(t).(fa{j})), nd = nd+1; end, end
            if nd, fprintf('    results.%s differs in %d/%d trials\n',fa{j},nd,n); end
        end
    end
    if isfield(A,'config') && isfield(B,'config') && ~isequaln(A.config,B.config)
        cf = fieldnames(A.config);
        for j=1:numel(cf), if ~isequaln(A.config.(cf{j}),B.config.(cf{j})), fprintf('    config.%s: %s vs %s\n',cf{j},string(jsonencode(A.config.(cf{j}))),string(jsonencode(B.config.(cf{j})))); end, end
    end
    if isfield(A,'runDate'), fprintf('  runDate old=%s rds=%s\n', string(A.runDate), string(B.runDate)); end
end
