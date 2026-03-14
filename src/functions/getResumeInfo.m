function [configIDs, lastCompletedIdx] = getResumeInfo(conn, jobIdentifier)
% GETRESUMEINFO Get information needed to resume a job
%
% This function retrieves the necessary information from the job_checkpoints
% table to resume a previously interrupted job.
%
% Parameters:
%   conn - SQLite database connection
%   jobIdentifier - Unique identifier for the job
%
% Returns:
%   configIDs - Array of configuration IDs for this job
%   lastCompletedIdx - Index of the last successfully completed configuration
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Get information needed to resume a job
sqlquery = ['SELECT config_ids, last_completed_idx FROM job_checkpoints ' ...
            'WHERE job_id = ''' jobIdentifier ''''];
result = fetch(conn, sqlquery);

% Extract results based on the result format
if istable(result)
    configIDsStr = result{1, 'config_ids'};
    lastCompletedIdx = double(result{1, 'last_completed_idx'});
elseif iscell(result) && size(result, 2) >= 2
    configIDsStr = result{1, 1};
    lastCompletedIdx = double(result{1, 2});
else
    configIDsStr = result(1, 1);
    lastCompletedIdx = double(result(1, 2));
end

% Parse the comma-separated config IDs string into a numeric array
configIDs = cellfun(@str2double, strsplit(configIDsStr, ','));
end
