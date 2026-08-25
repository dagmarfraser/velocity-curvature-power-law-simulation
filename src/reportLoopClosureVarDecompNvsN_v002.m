function report = reportLoopClosureVarDecompNvsN_v002(opts)
% REPORTLOOPCLOSUREVARDECOMPNVSN_V002  N=20 vs N=200 beta_gen* report.
%
% v002 adds a second per-dataset table using cluster_primary/test_primary
% (SG-LMLS, PP_PRIMARY=4 in loopClosureVarDecomp_v011's own indexing)
% alongside v001's constellation-median table. Prompted directly by v001's
% own output: Hickman HALO's constellation-median point estimate (0.2823)
% sits ~0.017 away from Table 6b's own reported value (0.299), a much
% larger gap than any other dataset showed (all within ~0.003) -- worth
% resolving properly rather than assuming it's noise. The paper's Table 6b
% uses SG-LMLS specifically as its deployment "primary" pipeline (Sec 7,
% "SG-LMLS as primary, three justifications"), a different quantity from
% this script's own internal "constellation median (primary)" label.
% Printing both columns side by side settles which one Table 6b's existing
% numbers actually are, rather than guessing.
%
% Everything else is unchanged from v001.
%
% Fraser, D.S. (2026)

    arguments
        opts.N20File  (1,1) string = "loopClosureVarDecomp_v011_all6_N20.mat"
        opts.N200File (1,1) string = "loopClosureVarDecomp_v011_all6_N200.mat"
    end

    srcDir = fileparts(mfilename("fullpath"));
    cd(srcDir);

    n20Path  = fullfile(srcDir, opts.N20File);
    n200Path = fullfile(srcDir, opts.N200File);
    if ~isfile(n20Path)
        error("reportLoopClosureVarDecompNvsN:noN20File", "%s", ...
            sprintf("FAILED PATH: %s not found. Run loopClosureVarDecomp_v011(''UseConfirmatoryCorpus'',false) first.", n20Path));
    end
    if ~isfile(n200Path)
        error("reportLoopClosureVarDecompNvsN:noN200File", "%s", ...
            sprintf("FAILED PATH: %s not found. Run loopClosureVarDecomp_v011(''UseConfirmatoryCorpus'',true) first.", n200Path));
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
        error("reportLoopClosureVarDecompNvsN:selfCheckFailed", "%s", sprintf( ...
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
    fprintf("POOLED beta_gen*, HKSJ method, constellation median (primary) -- N=20 vs N=200\n");
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
    % "constellation median" = this script's own internal PP_PRIMARY=4
    % naming quirk aside, this is the median across all six pipelines
    % (out(k).cluster_median / test_median).
    % "SG-LMLS" = out(k).cluster_primary / test_primary (PP_PRIMARY index 4
    % = SG-LMLS in PP_NAMES_V007), which is what the PAPER itself calls its
    % primary/deployment pipeline (Sec 7). Printed separately and clearly
    % labeled because these are genuinely different quantities that happen
    % to share the word "primary" in two different naming schemes -- Table
    % 6b's own numbers are checked against both below.
    nDS20  = numel(S20.out);
    nDS200 = numel(S200.out);
    if nDS20 ~= nDS200
        error("reportLoopClosureVarDecompNvsN:datasetCountMismatch", "%s", sprintf( ...
            "N20 has %d datasets, N200 has %d -- cannot align per-dataset rows.", nDS20, nDS200));
    end

    dsRows = cell(nDS20, 1);
    for k = 1:nDS20
        lab20  = string(S20.out(k).label);
        lab200 = string(S200.out(k).label);
        if lab20 ~= lab200
            error("reportLoopClosureVarDecompNvsN:datasetOrderMismatch", "%s", sprintf( ...
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
    fprintf("PER-DATASET beta_gen*, CONSTELLATION MEDIAN (across all six pipelines) -- N=20 vs N=200\n");
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
    fprintf("PER-DATASET beta_gen*, SG-LMLS (paper's own deployment ''primary'' pipeline, Sec 7) -- N=20 vs N=200\n");
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
