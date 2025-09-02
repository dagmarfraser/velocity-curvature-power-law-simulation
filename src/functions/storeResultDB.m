function storeResultDB(conn, configID, workerID, results)
% STORERESULTDB Stores power law analysis results in the database
%
% This function inserts the results of power law analysis for a specific
% parameter configuration into the database.
%
% Parameters:
%   conn - SQLite database connection
%   configID - Configuration ID for these results
%   workerID - Worker ID that processed this configuration
%   results - Structure containing:
%     .beta - Calculated beta value
%     .vgf - Calculated VGF value
%     .duration - Duration of the trajectory
%     .err_madirolas - Error relative to Madirolas method
%     .err_curvature - Curvature calculation error
%     .success - Boolean indicating success or failure
%     .processing_time - Time taken to process this configuration (optional)
%     .error_message - Error message if processing failed (optional)
%
% Created May 2025
% Updated May 2025 - Updated for v027 to match v023 parameter outputs
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Prepare values, handling NaN and missing fields
if isfield(results, 'beta') && ~isempty(results.beta) && ~isnan(results.beta)
    beta = num2str(results.beta);
else
    beta = 'NULL';
end

if isfield(results, 'vgf') && ~isempty(results.vgf) && ~isnan(results.vgf)
    vgf = num2str(results.vgf);
else
    vgf = 'NULL';
end

if isfield(results, 'duration') && ~isempty(results.duration) && ~isnan(results.duration)
    duration = num2str(results.duration);
else
    duration = 'NULL';
end

if isfield(results, 'err_madirolas') && ~isempty(results.err_madirolas) && ~isnan(results.err_madirolas)
    err_madirolas = num2str(results.err_madirolas);
else
    err_madirolas = 'NULL';
end

if isfield(results, 'err_curvature') && ~isempty(results.err_curvature) && ~isnan(results.err_curvature)
    err_curvature = num2str(results.err_curvature);
else
    err_curvature = 'NULL';
end

if isfield(results, 'r_squared') && ~isempty(results.r_squared) && ~isnan(results.r_squared)
    r_squared = num2str(results.r_squared);
else
    r_squared = 'NULL';
end

if isfield(results, 'success')
    if islogical(results.success)
        success = num2str(double(results.success));
    else
        success = num2str(results.success);
    end
else
    success = '0';
end

if isfield(results, 'processing_time') && ~isempty(results.processing_time) && ~isnan(results.processing_time)
    processing_time = num2str(results.processing_time);
else
    processing_time = 'NULL';
end

if isfield(results, 'error_message') && ~isempty(results.error_message)
    % Escape any single quotes to prevent SQL injection
    error_message = strrep(results.error_message, '''', '''''');
    error_message = ['''' error_message ''''];
else
    error_message = 'NULL';
end

% Construct the SQL query
sqlquery = ['INSERT INTO results ' ...
          '(config_id, worker_id, beta, vgf, duration, ' ...
          'err_madirolas, err_curvature, r_squared, success, ' ...
          'processing_time, error_message, processed_time) ' ...
          'VALUES (' num2str(configID) ', ' num2str(workerID) ', ' ...
          beta ', ' vgf ', ' duration ', ' ...
          err_madirolas ', ' err_curvature ', ' r_squared ', ' ...
          success ', ' processing_time ', ' ...
          error_message ', datetime(''now''))'];

% Execute the query
execute(conn, sqlquery);
end
