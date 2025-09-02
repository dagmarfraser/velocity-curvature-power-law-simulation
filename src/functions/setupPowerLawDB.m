function setupPowerLawDB(dbFile)
% SETUPPOWERLAWDB Creates a SQLite database for power law analysis results
%
% This function initializes a SQLite database with the necessary tables
% for storing power law analysis configurations and results. If the database
% already exists, it opens it without recreating to preserve existing data.
%
% Parameters:
%   dbFile - Path to the SQLite database file to create or use
%
% Created May 2025
% Updated May 2025 - Removed noise_color column, added trial_num
% Updated May 2025 - Modified to preserve existing database
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Check if Database Toolbox is available
if ~license('test', 'Database_Toolbox')
    error('Database Toolbox is required for database implementation');
end

% Create directory if it doesn't exist
[dbDir, ~, ~] = fileparts(dbFile);
if ~isempty(dbDir) && ~exist(dbDir, 'dir')
    mkdir(dbDir);
end

try
    % Create database or connect to existing
    if exist(dbFile, 'file')
        conn = sqlite(dbFile);
        fprintf('Connected to existing SQLite database at: %s\n', dbFile);
    else
        conn = sqlite(dbFile, 'create');
        fprintf('Created new SQLite database at: %s\n', dbFile);
    end
    
    % Create tables for storing parameter configurations (if they don't exist)
    sqlquery = ['CREATE TABLE IF NOT EXISTS param_configs (' ...
                'config_id INTEGER PRIMARY KEY AUTOINCREMENT, ' ...
                'shape_type INTEGER, ' ...           % Angular frequency (phi)
                'sampling_rate REAL, ' ...           % Hz
                'noise_magnitude REAL, ' ...         % mm standard deviation
                'noise_type INTEGER, ' ...           % Noise spectral color, as in v023
                'filter_type INTEGER, ' ...          % Differentiation and filtering method
                'filter_params TEXT, ' ...           % JSON string of filter parameters
                'regress_type INTEGER, ' ...         % Regression method
                'generated_beta REAL, ' ...          % Beta value used for trajectory generation
                'vgf_value REAL, ' ...               % VGF value used for trajectory generation
                'trial_num INTEGER, ' ...            % Trial number for repeatability
                'creation_time TEXT)'];              % When the configuration was created
    execute(conn, sqlquery);
    
    % Create table for storing results
    sqlquery = ['CREATE TABLE IF NOT EXISTS results (' ...
                'result_id INTEGER PRIMARY KEY AUTOINCREMENT, ' ...
                'config_id INTEGER, ' ...            % Foreign key to param_configs
                'worker_id INTEGER, ' ...            % ID of the worker that processed this config
                'beta REAL, ' ...                    % Calculated beta value
                'vgf REAL, ' ...                     % Calculated VGF value
                'duration REAL, ' ...                % Duration of the trajectory
                'err_madirolas REAL, ' ...           % Error relative to Madirolas method
                'err_curvature REAL, ' ...           % Curvature calculation error
                'r_squared REAL, ' ...               % R-squared value for the fit
                'success INTEGER, ' ...              % Boolean indicating success (1) or failure (0)
                'processing_time REAL, ' ...         % Time taken to process this configuration
                'error_message TEXT, ' ...           % Error message if processing failed
                'processed_time TEXT, ' ...          % When the result was computed
                'FOREIGN KEY (config_id) REFERENCES param_configs(config_id))'];
    execute(conn, sqlquery);
    
    % Create checkpoint table for tracking job progress
    sqlquery = ['CREATE TABLE IF NOT EXISTS job_checkpoints (' ...
                'job_id TEXT PRIMARY KEY, ' ...            % Job identifier
                'config_ids TEXT, ' ...                    % Serialized config IDs
                'last_completed_idx INTEGER, ' ...         % Last completed configuration index
                'total_configs INTEGER, ' ...              % Total number of configurations
                'start_time TEXT, ' ...                    % When the job started
                'last_update_time TEXT, ' ...              % When the job was last updated
                'completion_time TEXT, ' ...               % When the job was completed
                'status TEXT)'];                           % Job status (running, completed, failed)
    execute(conn, sqlquery);
    
    % Create indices for improved query performance
    execute(conn, 'CREATE INDEX IF NOT EXISTS idx_config_id ON results(config_id)');
    execute(conn, 'CREATE INDEX IF NOT EXISTS idx_shape_type ON param_configs(shape_type)');
    execute(conn, 'CREATE INDEX IF NOT EXISTS idx_filter_type ON param_configs(filter_type)');
    execute(conn, 'CREATE INDEX IF NOT EXISTS idx_regress_type ON param_configs(regress_type)');
    execute(conn, 'CREATE INDEX IF NOT EXISTS idx_trial_num ON param_configs(trial_num)');
    
    % Ensure any transactions are committed - only if needed
    try
        % Check if a transaction is active before committing
        status = fetch(conn, 'PRAGMA transaction_active');
        if istable(status) && status{1,1} == 1 || iscell(status) && status{1} == 1 || status(1) == 1
            execute(conn, 'COMMIT');
        end
    catch
        % Ignore errors checking transaction status
    end
    
    % Close connection
    close(conn);
catch ME
    % Ensure connection is closed even if an error occurs
    if exist('conn', 'var') && ~isempty(conn)
        try
            % Attempt to rollback any transactions
            execute(conn, 'ROLLBACK');
        catch
            % Ignore errors during rollback
        end
        
        try
            close(conn);
        catch
            % Ignore errors during close
        end
    end
    
    % Re-throw the original error
    rethrow(ME);
end
end