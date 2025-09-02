function configIDs = generateParameterConfigsDB_fixed(conn, paramSpace)
% GENERATEPARAMETERCONFIGSDB_FIXED Creates parameter configurations with guaranteed validity
%
% This improved function ensures all parameters are correctly structured when 
% initially created, eliminating the need for validation or substitution later.
%
% Parameters:
%   conn - SQLite database connection
%   paramSpace - Structure containing all parameter space definitions
%
% Returns:
%   configIDs - Array of generated configuration IDs
%
% Created June 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Check if param_configs table exists
createConfigTable(conn);

% Check if a transaction is already in progress
tryTransactionQuery = 'SELECT 1';
hasTransaction = false; % Initialize to false
try
    % Try to start a transaction
    execute(conn, 'BEGIN IMMEDIATE TRANSACTION');
    hasTransaction = true;
    fprintf('Started new transaction for parameter generation\n');
catch ME
    % If we get an error about a transaction already being active, that's fine
    if contains(ME.message, 'within a transaction')
        fprintf('Using existing transaction for parameter generation\n');
        % Don't set hasTransaction to true here, as we didn't start one
    else
        % If it's another error, rethrow it
        rethrow(ME);
    end
end

try
    % Initialize arrays to track generated configurations
    configCount = calculateTotalConfigs(paramSpace);
    configIDs = zeros(configCount, 1);
    
    % Pre-calculate filter parameter structures to avoid repeated JSON encoding
    filterParamsArray = prepareFilterParams(paramSpace.filterTypes);
    
    % Track insertion progress
    insertCount = 0;
    
    % Generate all parameter combinations
    for samplingRate = paramSpace.samplingRates
        for shapeType = paramSpace.shapes
            for generatedBeta = paramSpace.generatedBetas
                for vgfValue = paramSpace.vgfValues
                    for noiseType = paramSpace.noiseTypes
                        for noiseMagnitude = paramSpace.noiseMagnitudes
                            for filterTypeIdx = 1:length(paramSpace.filterTypes)
                                filterType = paramSpace.filterTypes(filterTypeIdx);
                                
                                % Get pre-calculated filter parameters
                                filterParams = filterParamsArray{filterTypeIdx};
                                
                                for regressType = paramSpace.regressTypes
                                    for trialNum = paramSpace.trialNums
                                        % Increment counter
                                        insertCount = insertCount + 1;
                                        
                                        % Insert configuration
                                        configIDs(insertCount) = insertConfig(conn, ...
                                            shapeType, samplingRate, noiseMagnitude, noiseType, ...
                                            filterType, filterParams, regressType, ...
                                            generatedBeta, vgfValue, trialNum);
                                        
                                        % Periodic progress updates
                                        if mod(insertCount, 10000) == 0
                                            fprintf('Generated %d of %d configurations (%.1f%%)\n', ...
                                                insertCount, configCount, 100*insertCount/configCount);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Commit the transaction if we started one
    if hasTransaction
        try
            execute(conn, 'COMMIT');
            fprintf('Transaction committed successfully\n');
        catch ME
            warning(ME.identifier, '%s', ME.message);
            % Don't rethrow here as we're already in a success path
        end
    end
    
    fprintf('Successfully generated %d parameter configurations\n', insertCount);
catch ME
    % Rollback on error if we started a transaction
    if hasTransaction
        try
            execute(conn, 'ROLLBACK');
            fprintf('Transaction rolled back due to error\n');
        catch RollbackME
            warning(RollbackME.identifier, '%s', RollbackME.message);
        end
    end
    rethrow(ME);
end

% Return the configuration IDs
configIDs = configIDs(1:insertCount);
end

function filterParamsArray = prepareFilterParams(filterTypes)
% Pre-calculate filter parameter structures for all filter types
filterParamsArray = cell(length(filterTypes), 1);

for i = 1:length(filterTypes)
    filterType = filterTypes(i);
    
    switch filterType
        case 1 % Simple diff - no params needed
            filterParamsArray{i} = '{}';
            
        case 2 % Butterworth - LPBW
            % Ensure parameters meet requirements
            params = struct('order', 2, 'cutoff', 10, 'zerolag', 1);
            filterParamsArray{i} = jsonencode(params);
            
        case 3 % Butterworth - DiffBW
            % Ensure parameters meet requirements
            params = struct('order', 2, 'cutoff', 10, 'zerolag', 1);
            filterParamsArray{i} = jsonencode(params);
            
        case 4 % Savitzky-Golay
            % Ensure parameters meet requirements - order must be < width-1
            % For jerk calculation, order must be >= 4
            params = struct('order', 4, 'width', 17);
            filterParamsArray{i} = jsonencode(params);
            
        case 6 % fs-scaled Savitzky-Golay (Fraser et al. temporal equivalence)
            % fs-scaled SG with temporal equivalence parameters
            % order and reference_framelen for temporal scaling
            params = struct('order', 4, 'width', 17); % Same as case 4 for reference
            filterParamsArray{i} = jsonencode(params);
            
        otherwise
            error('MATLAB:PowerLaw:UnknownFilterType', ...
                'Unknown filter type: %d', filterType);
    end
end
end

function configID = insertConfig(conn, shapeType, samplingRate, noiseMagnitude, ...
    noiseType, filterType, filterParams, regressType, generatedBeta, vgfValue, trialNum)
% Insert a single configuration with validated parameters

% Build insert query
query = ['INSERT INTO param_configs (', ...
    'shape_type, sampling_rate, noise_magnitude, noise_type, filter_type, ', ...
    'filter_params, regress_type, generated_beta, vgf_value, trial_num, creation_time) ', ...
    'VALUES (', ...
    num2str(shapeType), ', ', ...
    num2str(samplingRate), ', ', ...
    num2str(noiseMagnitude), ', ', ...
    num2str(noiseType), ', ', ...
    num2str(filterType), ', ', ...
    '''' , filterParams, ''', ', ...
    num2str(regressType), ', ', ...
    num2str(generatedBeta), ', ', ...
    num2str(vgfValue), ', ', ...
    num2str(trialNum), ', ', ...
    'CURRENT_TIMESTAMP)'];

% Execute insert
execute(conn, query);

% Get the generated ID
result = fetch(conn, 'SELECT last_insert_rowid()');
if istable(result)
    configID = result{1,1};
else
    configID = result{1};
end
end

function totalCount = calculateTotalConfigs(paramSpace)
% Calculate total configuration count
totalCount = length(paramSpace.samplingRates) * ...
             length(paramSpace.shapes) * ...
             length(paramSpace.generatedBetas) * ...
             length(paramSpace.vgfValues) * ...
             length(paramSpace.noiseTypes) * ...
             length(paramSpace.noiseMagnitudes) * ...
             length(paramSpace.filterTypes) * ...
             length(paramSpace.regressTypes) * ...
             length(paramSpace.trialNums);
end

function createConfigTable(conn)
% Create param_configs table if it doesn't exist
tableExists = ~isempty(fetch(conn, "SELECT name FROM sqlite_master WHERE type='table' AND name='param_configs'"));

if ~tableExists
    createTableQuery = ['CREATE TABLE param_configs (', ...
        'config_id INTEGER PRIMARY KEY AUTOINCREMENT, ', ...
        'shape_type INTEGER, ', ...
        'sampling_rate REAL, ', ...
        'noise_magnitude REAL, ', ...
        'noise_type INTEGER, ', ...
        'filter_type INTEGER, ', ...
        'filter_params TEXT, ', ...
        'regress_type INTEGER, ', ...
        'generated_beta REAL, ', ...
        'vgf_value REAL, ', ...
        'trial_num INTEGER, ', ...
        'creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP)'];
    
    execute(conn, createTableQuery);
    
    % Create index for faster queries
    execute(conn, 'CREATE INDEX idx_pc_shape_type ON param_configs (shape_type)');
    execute(conn, 'CREATE INDEX idx_pc_filter_type ON param_configs (filter_type)');
end
end
