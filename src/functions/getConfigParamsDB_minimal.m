function params = getConfigParamsDB_minimal(conn, configID)
% GETCONFIGPARAMSDB_MINIMAL Retrieves parameters without validation or substitution
%
% This lightweight function simply retrieves parameters as they are stored
% in the database without any modification or validation, assuming they
% were created correctly.
%
% Parameters:
%   conn - SQLite database connection
%   configID - Configuration ID to retrieve
%
% Returns:
%   params - Structure containing parameter values exactly as stored
%
% Created June 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Fetch the configuration row
query = ['SELECT * FROM param_configs WHERE config_id = ' num2str(configID)];
result = fetch(conn, query);

if isempty(result)
    error('MATLAB:PowerLaw:MissingConfig', 'Configuration ID %d not found in database.', configID);
end

% Convert to struct for easier access - no validation
params = struct();

if istable(result)
    % Convert table to struct
    params = table2struct(result);
    params = params(1); % Get first row
else
    % Convert cell array to struct
    params.config_id = result{1, 1};
    params.shape_type = result{1, 2};
    params.sampling_rate = result{1, 3};
    params.noise_magnitude = result{1, 4};
    params.noise_type = result{1, 5};
    params.filter_type = result{1, 6};
    params.filter_params = result{1, 7};
    params.regress_type = result{1, 8};
    params.generated_beta = result{1, 9};
    params.vgf_value = result{1, 10};
    params.trial_num = result{1, 11};
    params.creation_time = result{1, 12};
end

% Parse filter parameters JSON only if necessary
if isfield(params, 'filter_params') && ~isempty(params.filter_params) && ...
        (ischar(params.filter_params) || isstring(params.filter_params))
    % Minimal JSON parsing - no validation
    try
        params.filter_params = jsondecode(params.filter_params);
    catch
        % Don't crash - just return the string and let the calling function handle it
    end
end

% Ensure numeric fields are double precision
params.filter_type = double(params.filter_type);
params.sampling_rate = double(params.sampling_rate);
params.noise_magnitude = double(params.noise_magnitude);
params.noise_type = double(params.noise_type);
params.regress_type = double(params.regress_type);
params.generated_beta = double(params.generated_beta);
params.vgf_value = double(params.vgf_value);
params.trial_num = double(params.trial_num);
end
