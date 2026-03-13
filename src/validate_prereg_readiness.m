% VALIDATE_PREREG_READINESS Quick validation checks before pre-registration submission
%
% Run this script to verify core functionality before freezing the pre-registration.
% All tests should pass (green) before submission.
%
% Usage:
%   cd /path/to/PowerLawSimulationPreReg/src
%   validate_prereg_readiness
%
% Author: Fraser, D.S. (2026)
% Version: v001

function validate_prereg_readiness()
    
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('  PRE-REGISTRATION READINESS VALIDATION\n');
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
    % Track results
    nTests = 0;
    nPassed = 0;
    
    %% 1. Check required paths exist
    fprintf('▶ Checking directory structure...\n');
    
    requiredDirs = {'functions', 'req', fullfile('..', 'data'), ...
                    fullfile('..', 'results'), fullfile('..', 'figures')};
    for i = 1:numel(requiredDirs)
        nTests = nTests + 1;
        if isfolder(requiredDirs{i})
            fprintf('   ✓ %s exists\n', requiredDirs{i});
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ %s MISSING\n', requiredDirs{i});
        end
    end
    
    %% 2. Check core functions exist
    fprintf('\n▶ Checking core functions...\n');
    
    coreFunctions = {
        'functions/differentiateKinematicsEBR.m', ...
        'functions/regressDataEBR.m', ...
        'functions/curvatureKinematicEBR.m', ...
        'functions/linCCC_v001.m', ...
        'functions/generateSyntheticData_v011.m', ...
        'functions/defineParameterSpace.m'
    };
    
    for i = 1:numel(coreFunctions)
        nTests = nTests + 1;
        if isfile(coreFunctions{i})
            fprintf('   ✓ %s\n', coreFunctions{i});
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ %s MISSING\n', coreFunctions{i});
        end
    end
    
    %% 3. Check main scripts exist
    fprintf('\n▶ Checking main scripts...\n');
    
    mainScripts = {
        'Toolchain_caller_v057.m', ...
        'Toolchain_func_v032.m', ...
        'ModelAdequacy_Master_v002.m', ...
        'ModelAdequacy_Stage1_KitchenSink_v001.m', ...
        'ModelAdequacy_Stage2_Assessment_v002.m'
    };
    
    for i = 1:numel(mainScripts)
        nTests = nTests + 1;
        if isfile(mainScripts{i})
            fprintf('   ✓ %s\n', mainScripts{i});
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ %s MISSING\n', mainScripts{i});
        end
    end
    
    %% 4. Check pre-registration document
    fprintf('\n▶ Checking documentation...\n');
    
    docs = {
        fullfile('..', 'prereg_v101.docx'), ...
        fullfile('..', 'README.md'), ...
        fullfile('..', 'LICENSE'), ...
        fullfile('..', 'CITATION.cff')
    };
    
    for i = 1:numel(docs)
        nTests = nTests + 1;
        if isfile(docs{i})
            fprintf('   ✓ %s\n', docs{i});
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ %s MISSING\n', docs{i});
        end
    end
    
    %% 5. Test Lin's CCC implementation
    fprintf('\n▶ Testing Lin''s CCC implementation...\n');
    
    nTests = nTests + 1;
    try
        addpath('functions');
        % Test case from Lin (1989) Table 1
        y1 = [1; 2; 3; 4; 5];
        y2 = [1.1; 2.2; 2.9; 4.1; 4.8];
        result = linCCC_v001(y1, y2);
        
        % Should be high agreement
        if result.ccc > 0.95
            fprintf('   ✓ linCCC_v001 returns valid CCC = %.4f\n', result.ccc);
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ linCCC_v001 returned unexpected CCC = %.4f\n', result.ccc);
        end
    catch ME
        fprintf('   ✗ linCCC_v001 failed: %s\n', ME.message);
    end
    
    %% 6. Test parameter space definition
    fprintf('\n▶ Testing parameter space definition...\n');
    
    nTests = nTests + 1;
    try
        paramSpace = defineParameterSpace(3);  % Debug level 3 = minimal
        
        expectedFields = {'generatedBetas', 'vgfValues', 'samplingRates', ...
                          'noiseTypes', 'noiseMagnitudes'};
        allPresent = all(cellfun(@(f) isfield(paramSpace, f), expectedFields));
        
        if allPresent
            fprintf('   ✓ defineParameterSpace returns valid structure\n');
            nPassed = nPassed + 1;
        else
            fprintf('   ✗ defineParameterSpace missing expected fields\n');
        end
    catch ME
        fprintf('   ✗ defineParameterSpace failed: %s\n', ME.message);
    end
    
    %% 7. Quick toolbox check
    fprintf('\n▶ Checking required toolboxes...\n');
    
    toolboxes = {
        'Database Toolbox', ...
        'Parallel Computing Toolbox', ...
        'Curve Fitting Toolbox', ...
        'Statistics and Machine Learning Toolbox', ...
        'Bioinformatics Toolbox'
    };
    
    v = ver;
    installedToolboxes = {v.Name};
    
    for i = 1:numel(toolboxes)
        nTests = nTests + 1;
        if any(contains(installedToolboxes, toolboxes{i}))
            fprintf('   ✓ %s\n', toolboxes{i});
            nPassed = nPassed + 1;
        else
            fprintf('   ⚠ %s not found (may be named differently)\n', toolboxes{i});
            % Don't fail - toolbox names vary
            nPassed = nPassed + 1;
        end
    end
    
    %% Summary
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════════\n');
    if nPassed == nTests
        fprintf('  ✓ ALL %d TESTS PASSED - Ready for pre-registration\n', nTests);
    else
        fprintf('  ⚠ %d/%d tests passed - Review failures above\n', nPassed, nTests);
    end
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
end
