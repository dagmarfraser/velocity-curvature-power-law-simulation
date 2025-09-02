function jobExists = checkJobExists(conn, jobIdentifier)
% CHECKJOBEXISTS Check if a job exists in the checkpoint table
%
% This function checks if a job with the provided identifier
% exists in the job_checkpoints table.
%
% Parameters:
%   conn - SQLite database connection
%   jobIdentifier - Unique identifier for the job
%
% Returns:
%   jobExists - Boolean indicating if the job exists (true) or not (false)
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Check if a job exists in the checkpoint table
sqlquery = ['SELECT COUNT(*) FROM job_checkpoints ' ...
            'WHERE job_id = ''' jobIdentifier ''''];
result = fetch(conn, sqlquery);

% Handle different result formats
if istable(result)
    count_val = result{1,1}; % For table format
elseif iscell(result)
    count_val = result{1}; % For cell array format
else
    count_val = result(1); % For numeric array format
end

jobExists = count_val > 0;
end