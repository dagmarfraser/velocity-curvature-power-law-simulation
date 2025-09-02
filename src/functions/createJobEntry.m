function createJobEntry(conn, jobIdentifier, configIDs)
% CREATEJOBENTRY Create a new job entry in the checkpoint table
%
% This function adds a new job entry to the job_checkpoints table
% with the provided job identifier and configuration IDs.
%
% Parameters:
%   conn - SQLite database connection
%   jobIdentifier - Unique identifier for the job
%   configIDs - Array of configuration IDs to process in this job
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Serialize configIDs as comma-separated string
configIDsStr = strjoin(cellfun(@num2str, num2cell(configIDs), 'UniformOutput', false), ',');

% Insert job entry
sqlquery = ['INSERT INTO job_checkpoints ' ...
            '(job_id, config_ids, last_completed_idx, total_configs, start_time, status) ' ...
            'VALUES (''' jobIdentifier ''', ''' configIDsStr ''', 0, ' ...
            num2str(length(configIDs)) ', datetime(''now''), ''running'')'];
execute(conn, sqlquery);
end