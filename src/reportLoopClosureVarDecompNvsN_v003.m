function report = reportLoopClosureVarDecompNvsN_v003(opts)
% REPORTLOOPCLOSUREVARDECOMPNVSN_V003  N=20 vs N=200 beta_gen* report,
% against the CI-convention-corrected data (Findings #195/#196).
%
% WHY THIS VERSION EXISTS (2026-09-13). v001/v002 read
% loopClosureVarDecomp_v011's own saved N20/N200 outputs. Finding #196
% found v011's own ciWidth = ciLo-ciHi formula is correct only for legacy
% v007/v008 data and silently wrong-signed for the N200 branch (which
% loads v013/v014, a post-v009-format corpus using the opposite ciLo<=ciHi
% convention) -- confirmed directly against the saved
% loopClosureVarDecomp_v011_all6_N200.mat (pctMono for SG-LMLS = 0.0% for
% every dataset, the same degenerate signature Finding #196 traces to
% source). Point estimates (this script's own median/cluster_median/
% cluster_primary fields) are NOT affected -- they never read ciLo/ciHi --
% but any SE/CI width derived from the N200 branch's per-trial inversion
% uncertainty term is understated there, most consequentially for the
% SG-LMLS ("primary") columns v002 added, which have no cross-pipeline
% floor to fall back on (unlike constellation median).
%
% THE FIX, the only change from v002: default N20File/N200File point at
% loopClosureVarDecomp_v012_gated.m's own outputs
% (loopClosureVarDecomp_v012_gated_all6_legacy.mat /
% _confirmatory.mat) instead of v011's. That script's ciWidth already
% branches correctly on corpus convention (Finding #195/#196). Everything
% else -- the self-check, the pooled k=3/4/5/6 table, both per-dataset
% tables (constellation median and SG-LMLS) -- is byte-identical to v002.
%
% v011/v001/v002 and their existing .mat outputs are left completely
% untouched, per this project's version-lineage-immutability convention.
%
% Fraser, D.S. (2026)  v003

    arguments
        opts.N20File  (1,1) string = "loopClosureVarDecomp_v012_gated_all6_legacy.mat"
        opts.N200File (1,1) string = "loopClosureVarDecomp_v012_gated_all6_confirmatory.mat"
    end

    srcDir = fileparts(mfilename("fullpath"));
    cd(srcDir);

    n20Path  = fullfile(srcDir, opts.N20File);
    n200Path = fullfile(srcDir, opts.N200File);
    if ~isfile(n20Path)
        error("reportLoopClosureVarDecompNvsN_v003:noN20File", "%s", ...
            sprintf("FAILED PATH: %s not found. Run loopClosureVarDecomp_v012_gated(''DataSource'',''legacy'') first.", n20Path));
    end
    if ~isfile(n200Path)
        error("reportLoopClosureVarDecompNvsN_v003:noN200File", "%s", ...
            sprintf("FAILED PATH: %s not found. Run loopClosureVarDecomp_v012_gated(''DataSource'',''confirmatory'') first.", n200Path));
    end

    S20  = load(n20Path,  "out", "synthMed", "synthMed_k3", "synthMed_k4", "synthMed_k6");
    S200 = load(n200Path, "out", "synthMed", "synthMed_k3", "synthMed_k4", "synthMed_k6");

    %% --- Self-check against Finding #140's already-published number -------
    expected = struct("mPoolRE", 0.3102, "ciHKSJ", [0.2548, 0.3656], "dfHKSJ", 2);
    tol = 1e-3;
    k3 = S20.synthMed_k3;
    ok = abs(k3.mPoolRE - expected.mPoolRE) < tol && ...
         abs(k3.ciHKSJ(1) - expected.ciHKSJ(1)) < tol && ...
         abs(k3.ciHKSJ(2) - expected.ciHKSJ(2)) < tol && ...
         k3.dfHKSJ == expected.dfHKSJ;
    if ~ok
        error("reportLoopClosureVarDecompNvsN_v003:selfCheckFailed", "%s", sprintf( ...
            "%s does NOT reproduce Finding #140's published k=3 HKSJ number. " + ...
            "Got mPoolRE=%.4f, CI=[%.4f,%.4f], df=%d; expected mPoolRE=%.4f, " + ...
            "CI=[%.4f,%.4f], df=%d. Do not trust the N=200 comparison below until " + ...
            "this is understood -- either the N=20 file is stale/regenerated " + ...
            "differently, or something upstream changed.", ...
            opts.N20File, k3.mPoolRE, k3.ciHKSJ(1), k3.ciHKSJ(2), k3.dfHKSJ, ...
            expected.mPoolRE, expected.ciHKSJ(1), expected.ciHKSJ(2), expected.dfHKSJ));
    end
    fprintf("*** SELF-CHECK PASSED: %s exactly reproduces Finding #140's published ", opts.N20File);
    fprintf("k=3 HKSJ number (mPoolRE=%.4f, CI=[%.4f,%.4f], df=%d). ***\n\n", ...
        k3.mPoolRE, k3.ciHKSJ(1), k3.ciHKSJ(2), k3.dfHKSJ);

    %% --- Pooled comparison: k=3/4/5/6, HKSJ ---------------------------------
    views = struct( ...
        "label",    {"k=3 (headline candidate)", "k=4 (Hickman merged)", ...
                      "k=5 (historical, Finding #137 era)", "k=6 (+Zarandi)"}, ...
        "fieldMed", {"synthMed_k3", "synthMed_k4", "synthMed", "synthMed_k6"});

    fprintf("%s\n", repmat('=', 1, 100));
    fprintf("POOLED beta_gen*, HKSJ method, constellation median (primary) -- N=20 vs N=200 [CI-CORRECTED]\n");
    fprintf("%s\n", repmat('=', 1, 100));
    fprintf("%-38s %8s %20s %4s   %8s %20s %4s  %s\n", ...
        "View", "N20 est", "N20 95%% CI", "df", "N200 est", "N200 95%% CI", "df", "excludes 1/3? (N20->N200)");
    fprintf("%s\n", repmat('-', 1, 100));

    poolRows = cell(numel(views), 1);
    for v = 1:numel(views)
        a = S20.(views(v).fieldMed);
        b = S200.(views(v).fieldMed);
        yn20  = "no "; if a.hksjExcludes,  yn20  = "YES"; end
        yn200 = "no "; if b.hksjExcludes,  yn200 = "YES"; end
        flag = ""; if a.hksjExcludes ~= b.hksjExcludes, flag = "  <-- CHANGED"; end
        fprintf("%-38s %8.4f %20s %4d   %8.4f %20s %4d  %s -> %s%s\n", ...
            views(v).label, a.mPoolRE, ciStr_local(a.ciHKSJ), a.dfHKSJ, ...
            b.mPoolRE, ciStr_local(b.ciHKSJ), b.dfHKSJ, yn20, yn200, flag);
        poolRows{v} = struct("view", string(views(v).label), ...
            "estN20", a.mPoolRE, "seN20", a.sePoolHKSJ, "ciN20", a.ciHKSJ, "dfN20", a.dfHKSJ, ...
            "pN20", a.poolP_HKSJ, "excludesN20", a.hksjExcludes, ...
            "estN200", b.mPoolRE, "seN200", b.sePoolHKSJ, "ciN200", b.ciHKSJ, "dfN200", b.dfHKSJ, ...
            "pN200", b.poolP_HKSJ, "excludesN200", b.hksjExcludes, ...
            "verdictChanged", a.hksjExcludes ~= b.hksjExcludes);
    end
    report.pooled = struct2table([poolRows{:}]);

    %% --- Per-dataset comparison: TWO estimators --------------------------
    nDS20  = numel(S20.out);
    nDS200 = numel(S200.out);
    if nDS20 ~= nDS200
        error("reportLoopClosureVarDecompNvsN_v003:datasetCountMismatch", "%s", sprintf( ...
            "N20 has %d datasets, N200 has %d -- cannot align per-dataset rows.", nDS20, nDS200));
    end

    dsRows = cell(nDS20, 1);
    for k = 1:nDS20
        lab20  = string(S20.out(k).label);
        lab200 = string(S200.out(k).label);
        if lab20 ~= lab200
            error("reportLoopClosureVarDecompNvsN_v003:datasetOrderMismatch", "%s", sprintf( ...
                "Dataset order mismatch at index %d: N20='%s' vs N200='%s'. Do not trust " + ...
                "row-by-row alignment below this point.", k, lab20, lab200));
        end
        dsRows{k} = struct("dataset", lab20, ...
            "medEstN20",  S20.out(k).cluster_median.pointEst,  "medCiN20",  S20.out(k).test_median.ciFM, ...
            "medEstN200", S200.out(k).cluster_median.pointEst, "medCiN200", S200.out(k).test_median.ciFM, ...
            "medDelta",   S200.out(k).cluster_median.pointEst - S20.out(k).cluster_median.pointEst, ...
            "sglmlsEstN20",  S20.out(k).cluster_primary.pointEst,  "sglmlsCiN20",  S20.out(k).test_primary.ciFM, ...
            "sglmlsEstN200", S200.out(k).cluster_primary.pointEst, "sglmlsCiN200", S200.out(k).test_primary.ciFM, ...
            "sglmlsDelta",   S200.out(k).cluster_primary.pointEst - S20.out(k).cluster_primary.pointEst);
    end
    report.perDataset = struct2table([dsRows{:}]);

    fprintf("\n%s\n", repmat('=', 1, 100));
    fprintf("PER-DATASET beta_gen*, CONSTELLATION MEDIAN (across all six pipelines) -- N=20 vs N=200 [CI-CORRECTED]\n");
    fprintf("%s\n", repmat('=', 1, 100));
    fprintf("%-16s %8s %20s   %8s %20s   %8s\n", ...
        "Dataset", "N20 est", "N20 t(M-1) 95%% CI", "N200 est", "N200 t(M-1) 95%% CI", "delta");
    fprintf("%s\n", repmat('-', 1, 90));
    for k = 1:nDS20
        r = dsRows{k};
        fprintf("%-16s %8.4f %20s   %8.4f %20s   %+8.4f\n", ...
            r.dataset, r.medEstN20, ciStr_local(r.medCiN20), r.medEstN200, ciStr_local(r.medCiN200), r.medDelta);
    end

    fprintf("\n%s\n", repmat('=', 1, 100));
    fprintf("PER-DATASET beta_gen*, SG-LMLS (paper's own deployment ''primary'' pipeline, Sec 7) -- N=20 vs N=200 [CI-CORRECTED]\n");
    fprintf("%s\n", repmat('=', 1, 100));
    fprintf("%-16s %8s %20s   %8s %20s   %8s\n", ...
        "Dataset", "N20 est", "N20 t(M-1) 95%% CI", "N200 est", "N200 t(M-1) 95%% CI", "delta");
    fprintf("%s\n", repmat('-', 1, 90));
    for k = 1:nDS20
        r = dsRows{k};
        fprintf("%-16s %8.4f %20s   %8.4f %20s   %+8.4f\n", ...
            r.dataset, r.sglmlsEstN20, ciStr_local(r.sglmlsCiN20), r.sglmlsEstN200, ciStr_local(r.sglmlsCiN200), r.sglmlsDelta);
    end

    fprintf("\n%s\n", repmat('=', 1, 100));

end

function s = ciStr_local(ci)
    s = sprintf("[%.4f, %.4f]", ci(1), ci(2));
end
