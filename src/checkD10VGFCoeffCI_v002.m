function checkD10VGFCoeffCI_v002(options)
% CHECKD10VGFCOEFFCI_V002  Verify D10's claimed VGF LMM coefficient CI half-width.
%
% Fixes two problems found running v001 on BlueBEAR:
%   1. The full-struct load took 842 s (~14 min) and then hit a validation
%      error, meaning that entire load had to be paid for again just to
%      iterate on the validation logic. v002 caches the raw `coeffs`
%      variable to a small local .mat file IMMEDIATELY after loading, before
%      any validation runs, so every subsequent run (this one included, if
%      it needs a v003) uses the cache and finishes in seconds.
%   2. v001 assumed `results.coefficients` is a `table` with columns
%      exactly {Name, Estimate, Lower, Upper} and errored otherwise, without
%      saying WHICH assumption failed. `LinearMixedModel.Coefficients` has
%      been a `dataset` array rather than a `table` in some MATLAB releases,
%      and this project's own D1 printed output (inspectLMMCoefficients_v001.m)
%      used CI_lo/CI_hi rather than Lower/Upper. v002 detects the actual
%      class and column names, tries known equivalents, and prints exactly
%      what it found before doing anything else -- so if even this doesn't
%      match, the diagnostic output tells you what does, from the cache,
%      without re-touching the ~16 GB source.
%
% docs/LMM_VGF_Desiderata_v003.md D6 states "maxCI = 1.85 mm/s" for the VGF
% secondary LMM's largest coefficient confidence half-width, cited without a
% Source/Script line (unlike every other D-item in that file). D10 in the
% manuscript's deviation table repeats this figure to justify treating
% Criterion 1 as PASS for the VGF model. Audit ledger, Session 114: never
% traced to `vgfLMM_v001_latest.mat` itself. This script does that trace.
%
% BlueBEAR usage (RDS is local there, not a network mount -- the Mac RDS
% mount has hung MATLAB on partial matfile reads). Only needs the full
% allocation on a cache-miss (first run, or forceReload=true):
%   sinteractive --ntasks=4 --time=00:30:00 --mem=64G
%   module load MATLAB/2024b && matlab
%   >> addpath(genpath('src')); checkD10VGFCoeffCI_v002()
%
% Output: checkD10VGFCoeffCI_v002_rawCoeffs.mat (cache, written right after
% load, before validation) and checkD10VGFCoeffCI_v002_results.mat (final
% verdict, written only on success) + console report.
%
% Author: Fraser, D.S. / Claude (2026), Round 3 D10 trace item.

    arguments
        options.matFile (1,1) string = fullfile( ...
            '/rds/projects/f/fraserds-mpo-evaluation/2026_prereg', ...
            'velocity-curvature-power-law-simulation-main', ...
            'velocity-curvature-power-law-simulation-main', 'src', ...
            'vgfLMM_v001_latest.mat')
        options.claimedMaxCI (1,1) double = 1.85
        options.vgfGenRangeMin (1,1) double = 90.0
        options.vgfGenRangeMax (1,1) double = 330.3
        options.tolerance (1,1) double = 0.01
        options.forceReload (1,1) logical = false
    end

    scriptDir = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(scriptDir, 'checkD10VGFCoeffCI_v002_rawCoeffs.mat');

    %% 1. Get coeffs -- from cache if available, else pay the expensive load once
    if isfile(cacheFile) && ~options.forceReload
        fprintf('Cache hit: %s\n', cacheFile);
        fprintf('(pass forceReload=true to bypass this and re-load the full source file)\n\n');
        cached = load(cacheFile, 'coeffs');
        coeffs = cached.coeffs;
    else
        if ~isfile(options.matFile)
            error('checkD10VGFCoeffCI:NoFile', ...
                'FAILED PATH: %s not found.', options.matFile);
        end
        fprintf('Cache miss. Loading full results struct from:\n  %s\n', options.matFile);
        fprintf('No per-field partial load is possible for a MAT-file struct (matfile\n');
        fprintf('objects support only one level of indexing) -- this pulls tableTrue\n');
        fprintf('(~17M rows) into memory too, as an unavoidable side effect. Expect\n');
        fprintf('several minutes on a ~16 GB file.\n');
        tLoad = tic;
        S = load(options.matFile, 'results');
        if ~isfield(S, 'results') || ~isfield(S.results, 'coefficients')
            error('checkD10VGFCoeffCI:NoCoefficients', '%s', sprintf( ...
                '%s loaded but has no results.coefficients field.', options.matFile));
        end
        coeffs = S.results.coefficients;
        clear S   % drop tableTrue and everything else immediately
        fprintf('  Loaded in %.1f s\n', toc(tLoad));

        % Checkpoint BEFORE validation: whatever coeffs turns out to be, it is
        % now cheap to re-inspect regardless of what the checks below find.
        save(cacheFile, 'coeffs');
        fprintf('  Cached to %s\n\n', cacheFile);
    end

    %% 2. Diagnose actual type and column names -- print before assuming anything
    fprintf('class(coeffs): %s\n', class(coeffs));
    if istable(coeffs)
        varNames = coeffs.Properties.VariableNames;
    elseif isa(coeffs, 'dataset')
        varNames = coeffs.Properties.VarNames;
        coeffs   = dataset2table(coeffs);   % normalise for uniform handling below
    else
        error('checkD10VGFCoeffCI:UnknownClass', '%s', sprintf( ...
            ['results.coefficients is class %s, neither table nor dataset. ' ...
             'Raw value is cached at %s -- inspect it directly (e.g. load it ' ...
             'and try fieldnames/properties on it) rather than re-running this ' ...
             'script blind.'], class(coeffs), cacheFile));
    end
    fprintf('Columns found: %s\n\n', strjoin(varNames, ', '));

    % Column names for fitlme Coefficients tables are usually
    % {Name, Estimate, SE, tStat, DF, pValue, Lower, Upper}, but this
    % project's own D1 printed output (inspectLMMCoefficients_v001.m) used
    % CI_lo/CI_hi -- try known equivalents rather than hard-failing on the
    % first guess. Disclosed here, not silent: the "Using columns" line
    % below states exactly which names were matched.
    nameCol  = firstMatch_local(varNames, {'Name','Term'});
    estCol   = firstMatch_local(varNames, {'Estimate'});
    lowerCol = firstMatch_local(varNames, {'Lower','CI_lo','CILower'});
    upperCol = firstMatch_local(varNames, {'Upper','CI_hi','CIUpper'});

    if isempty(nameCol) || isempty(estCol) || isempty(lowerCol) || isempty(upperCol)
        error('checkD10VGFCoeffCI:NoRecognisedColumns', '%s', sprintf( ...
            ['Could not find recognised Name/Estimate/Lower/Upper-equivalent ' ...
             'columns among: %s. Raw coefficients are cached at %s -- inspect ' ...
             'them directly and extend the candidate lists in this script ' ...
             'rather than re-running blind.'], strjoin(varNames, ', '), cacheFile));
    end
    fprintf('Using columns: name=%s, estimate=%s, lower=%s, upper=%s\n\n', ...
        nameCol, estCol, lowerCol, upperCol);

    %% 3. Compute and report
    halfWidth   = (coeffs.(upperCol) - coeffs.(lowerCol)) / 2;
    [maxHW, ix] = max(halfWidth);
    matches     = abs(maxHW - options.claimedMaxCI) < options.tolerance;
    vgfRange    = options.vgfGenRangeMax - options.vgfGenRangeMin;
    nameVals    = coeffs.(nameCol);

    fprintf('N coefficients:        %d\n', height(coeffs));
    fprintf('Max CI half-width:     %.4f mm/s\n', maxHW);
    fprintf('  Term:     %s\n', string(nameVals(ix)));
    fprintf('  Estimate: %.4f  Lower: %.4f  Upper: %.4f\n', ...
        coeffs.(estCol)(ix), coeffs.(lowerCol)(ix), coeffs.(upperCol)(ix));
    fprintf('\nD6/D10 claimed maxCI:  %.2f mm/s\n', options.claimedMaxCI);
    fprintf('Match (tol %.2f):      %d\n', options.tolerance, matches);
    fprintf('\nRelative uncertainty (max CI / VGF_gen range %.1f-%.1f mm/s): %.2f%%\n', ...
        options.vgfGenRangeMin, options.vgfGenRangeMax, 100*maxHW/vgfRange);

    [sortedHW, sortIx] = sort(halfWidth, 'descend');
    nTop = min(5, numel(sortedHW));
    fprintf('\nTop %d coefficients by CI half-width:\n', nTop);
    for k = 1:nTop
        i = sortIx(k);
        fprintf('  %-45s  halfWidth=%.4f  Estimate=%.4f\n', ...
            string(nameVals(i)), halfWidth(i), coeffs.(estCol)(i));
    end

    if ~matches
        warning('checkD10VGFCoeffCI:Mismatch', '%s', sprintf( ...
            ['Recomputed max CI half-width (%.4f mm/s) does not match D6/D10''s ' ...
             'claimed %.2f mm/s. Do not propagate the claimed figure further ' ...
             'until this is reconciled -- check whether D6''s number came from ' ...
             'a different coefficient set or model run.'], maxHW, options.claimedMaxCI));
    end

    outFile = fullfile(scriptDir, 'checkD10VGFCoeffCI_v002_results.mat');
    save(outFile, 'coeffs', 'maxHW', 'ix', 'matches', 'nameCol', 'estCol', ...
        'lowerCol', 'upperCol', 'options');
    fprintf('\nSaved: %s\n', outFile);

end

function m = firstMatch_local(varNames, candidates)
% Returns the first name in `candidates` present in `varNames`, or "" if none.
    m = "";
    for k = 1:numel(candidates)
        if any(strcmp(varNames, candidates{k}))
            m = candidates{k};
            return
        end
    end
end
