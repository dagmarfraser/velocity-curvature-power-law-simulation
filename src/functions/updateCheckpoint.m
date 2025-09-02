function updateCheckpoint(conn, jobIdentifier, lastCompletedIdx)
% UPDATECHECKPOINT Update the checkpoint with the last completed index
%
% This function updates the job_checkpoints table with the last completed
% configuration index and the current timestamp.
%
% Parameters:
%   conn - SQLite database connection
%   jobIdentifier - Unique identifier for the job
%   lastCompletedIdx - Index of the last successfully completed configuration
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Update the checkpoint with the last completed configuration index
sqlquery = ['UPDATE job_checkpoints ' ...
            'SET last_completed_idx = ' num2str(lastCompletedIdx) ', ' ...
            'last_update_time = datetime(''now'') ' ...
            'WHERE job_id = ''' jobIdentifier ''''];
execute(conn, sqlquery);
end