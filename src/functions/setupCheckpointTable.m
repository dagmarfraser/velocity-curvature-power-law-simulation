function setupCheckpointTable(conn)
% SETUPCHECKPOINTTABLE Create checkpoint table if it doesn't exist
%
% This function creates the job_checkpoints table in the SQLite database
% if it doesn't already exist. This table is used to track job progress
% for checkpoint and resume functionality.
%
% Parameters:
%   conn - SQLite database connection
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Create checkpoint table if it doesn't exist
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
end