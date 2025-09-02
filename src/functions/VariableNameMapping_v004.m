function varargout = VariableNameMapping_v004(varargin)
% VARIABLENAMEMAPPING_V004 - Enhanced Variable Name Mapping with Safe Comparison Functionality
%
% **NEW IN v004**: 
% Added safeComparison functionality to handle nominal vs numeric comparison issues
% throughout the ToolchainCaller framework. Builds on proven type detection patterns
% from v002-v003 that fixed the 47.1% validation failure rate.
%
% **CRITICAL FIX v003**: 
% FIXED variable 'i' collision issue that was overwriting MATLAB's complex number constant.
% This resolves the "All input arguments must be tables" error in pipeline context.
%
% **v003 CHANGES**:
%   - Replaced all 'for i =' loops with 'for idx =' to avoid complex number collision
%   - Maintains all v002 functionality while fixing workspace contamination
%   - Resolves Fraser's "Colon Theory" pipeline failure issue
%
% **CORE UTILITY FOR HIERARCHICAL DRILL-DOWN PIPELINE**
% Provides robust mapping between LME model coefficient names and actual data 
% variable names, with special handling for categorical variable encoding issues.
%
% **v002 ENHANCEMENT**: 
% FIXED validation logic to properly handle NOMINAL variables in addition to categorical.
% Resolves 47.1% validation failure rate identified in large variable space testing.
%
% **SUCCESSOR RESEARCH ENHANCEMENT**: 
% Addresses categorical variable naming challenges identified in large variable 
% space synthetic data analysis for Fraser et al. (2025) successor research.
%
% **KEY FUNCTIONALITY**:
%   - Maps LME coefficient names to data table variable names
%   - Handles categorical AND NOMINAL variable encoding complexities
%   - Provides validation and error handling for unmappable variables
%   - Supports interaction term decomposition
%   - Enables consistent variable naming across pipeline stages
%   - **NEW: Safe variable comparisons handling nominal/categorical vs numeric**
%
% **USAGE PATTERNS**:
%   % Basic coefficient to variable mapping
%   data_var = VariableNameMapping_v004('mapCoefficientToVariable', coeff_name, data_table);
%   
%   % Interaction term decomposition
%   components = VariableNameMapping_v004('decomposeInteraction', interaction_term, data_table);
%   
%   % Variable validation
%   isValid = VariableNameMapping_v004('validateMapping', coeff_name, data_table);
%   
%   % **NEW: Safe variable comparisons**
%   indices = VariableNameMapping_v004('safeComparison', data_table.filterType, 1, 'equal');
%   indices = VariableNameMapping_v004('safeComparison', data_table.regressionType, 2, 'equal');
%   
%   % Get mapping configuration
%   config = VariableNameMapping_v004('getMappingConfig');
%
% **INTEGRATION**: Designed for use in Stage 2+ of hierarchical drill-down pipeline
% to resolve coefficient naming issues in segmentation logic and nominal vs numeric
% comparison problems throughout ModelAdequacy framework.
%
% Author: Fraser, D.S. (2025)
% Version: v004 - Added safe comparison functionality for nominal vs numeric fixes

%% INPUT VALIDATION AND DISPATCH
if nargin < 1
    error('VariableNameMapping_v004:InvalidInput', 'Function requires at least one input argument');
end

action = varargin{1};

% Dispatch to appropriate function based on action
switch lower(action)
    case 'mapcoefficienttovariable'
        [varargout{1:nargout}] = mapCoefficientToVariable(varargin{2:end});
        
    case 'decomposeinteraction'
        [varargout{1:nargout}] = decomposeInteraction(varargin{2:end});
        
    case 'validatemapping'
        [varargout{1:nargout}] = validateVariableMapping(varargin{2:end});
        
    case 'getmappingconfig'
        [varargout{1:nargout}] = getMappingConfig(varargin{2:end});
        
    case 'extractvariablefromcoeff'
        [varargout{1:nargout}] = extractVariableFromCoeff(varargin{2:end});
        
    case 'handlecategoricalencoding'
        [varargout{1:nargout}] = handleCategoricalEncoding(varargin{2:end});
        
    case 'createfullmappingtable'
        [varargout{1:nargout}] = createFullMappingTable(varargin{2:end});
        
    case 'debugmapping'
        [varargout{1:nargout}] = debugMapping(varargin{2:end});
        
    case 'safecomparison'
        [varargout{1:nargout}] = safeVariableComparison(varargin{2:end});
        
    otherwise
        error('VariableNameMapping_v004:UnknownAction', ...
            'Unknown action: %s. Valid actions: mapCoefficientToVariable, decomposeInteraction, validateMapping, getMappingConfig, safeComparison', action);
end

end

%% NEW: SAFE COMPARISON FUNCTIONALITY

function indices = safeVariableComparison(variable, target_value, comparison_type, varargin)
% SAFEVARIABLECOMPARISON - Handle nominal/categorical vs numeric comparisons
%
% **ROBUST TYPE DETECTION**: Uses the same proven type detection patterns 
% from validateVariableMapping that fixed the 47.1% validation failure rate.
% Critical precedence: checks nominal FIRST since nominal vars return true for both checks.
%
% **PURPOSE**: Eliminates nominal vs numeric comparison errors throughout
% the ToolchainCaller framework by providing type-safe comparison operations.
%
% INPUTS:
%   variable - Variable to compare (nominal, categorical, or numeric)
%   target_value - Value to compare against  
%   comparison_type - 'equal', 'not_equal', 'greater', 'less', etc.
%   varargin - Optional parameters
%
% OUTPUT:
%   indices - Numeric indices where comparison is true
%
% USAGE:
%   indices = VariableNameMapping_v004('safeComparison', data_table.filterType, 1, 'equal');
%   indices = VariableNameMapping_v004('safeComparison', data_table.regressionType, 2, 'not_equal');

    % Input validation
    if nargin < 3
        error('VariableNameMapping_v004:InsufficientInputs', ...
            'safeVariableComparison requires variable, target_value, and comparison_type');
    end
    
    if isempty(variable)
        error('VariableNameMapping_v004:EmptyVariable', ...
            'Variable cannot be empty');
    end
    
    if ~ischar(comparison_type) && ~isstring(comparison_type)
        error('VariableNameMapping_v004:InvalidComparisonType', ...
            'comparison_type must be a string');
    end
    
    % Convert comparison_type to lowercase for consistent handling
    comparison_type = lower(char(comparison_type));
    
    % Use the same robust type detection as validateVariableMapping
    % CRITICAL: Check nominal FIRST since nominal vars return true for both checks
    if isa(variable, 'nominal')
        % Handle nominal variables (check nominal FIRST - critical precedence)
        try
            switch comparison_type
                case 'equal'
                    indices = find(variable == nominal(target_value));
                case 'not_equal'
                    indices = find(variable ~= nominal(target_value));
                otherwise
                    error('VariableNameMapping_v004:UnsupportedNominalComparison', ...
                        'Unsupported comparison "%s" for nominal variables. Supported: equal, not_equal', comparison_type);
            end
        catch ME
            if contains(ME.message, 'Undefined function')
                % Handle case where nominal conversion fails
                error('VariableNameMapping_v004:NominalConversionFailed', ...
                    'Cannot convert target_value %s to nominal for comparison', num2str(target_value));
            else
                rethrow(ME);
            end
        end
        
    elseif iscategorical(variable)
        % Handle categorical variables (non-nominal)
        try
            switch comparison_type
                case 'equal'
                    indices = find(variable == categorical(target_value));
                case 'not_equal'
                    indices = find(variable ~= categorical(target_value));
                otherwise
                    error('VariableNameMapping_v004:UnsupportedCategoricalComparison', ...
                        'Unsupported comparison "%s" for categorical variables. Supported: equal, not_equal', comparison_type);
            end
        catch ME
            if contains(ME.message, 'Undefined function')
                % Handle case where categorical conversion fails
                error('VariableNameMapping_v004:CategoricalConversionFailed', ...
                    'Cannot convert target_value %s to categorical for comparison', num2str(target_value));
            else
                rethrow(ME);
            end
        end
        
    elseif isnumeric(variable)
        % Handle numeric variables
        switch comparison_type
            case 'equal'
                indices = find(variable == target_value);
            case 'not_equal'
                indices = find(variable ~= target_value);
            case 'greater'
                indices = find(variable > target_value);
            case 'less'
                indices = find(variable < target_value);
            case 'greater_equal'
                indices = find(variable >= target_value);
            case 'less_equal'
                indices = find(variable <= target_value);
            otherwise
                error('VariableNameMapping_v004:UnsupportedNumericComparison', ...
                    'Unsupported comparison type: %s. Supported: equal, not_equal, greater, less, greater_equal, less_equal', comparison_type);
        end
    else
        error('VariableNameMapping_v004:UnsupportedDataType', ...
            'Unsupported data type for comparison: %s. Supported: nominal, categorical, numeric', class(variable));
    end
    
    % Ensure indices are returned as column vector for consistency
    if isrow(indices)
        indices = indices';
    end
    
    % Optional: Add debugging information if requested
    if nargin > 3 && isstruct(varargin{1}) && isfield(varargin{1}, 'debug') && varargin{1}.debug
        fprintf('SafeComparison Debug:\n');
        fprintf('  Variable type: %s\n', class(variable));
        fprintf('  Target value: %s\n', num2str(target_value));
        fprintf('  Comparison: %s\n', comparison_type);
        fprintf('  Found indices: %d matches\n', length(indices));
    end
end

%% CORE MAPPING FUNCTIONS (UNCHANGED FROM v003)

function mapped_variable = mapCoefficientToVariable(coefficient_name, data_table, varargin)
% Map LME model coefficient name to actual data table variable name
%
% **ENHANCED LOGIC**: Handles categorical variable encoding patterns and 
% interaction terms with robust error handling and validation.
%
% INPUTS:
%   coefficient_name - String, coefficient name from LME model
%   data_table - Table, the data table with actual variable names
%   varargin - Optional parameters for mapping configuration
%
% OUTPUT:
%   mapped_variable - String, corresponding data table variable name

    % Input validation
    if nargin < 2
        error('VariableNameMapping_v004:InsufficientInputs', ...
            'mapCoefficientToVariable requires coefficient_name and data_table');
    end
    
    if isempty(coefficient_name) || ~ischar(coefficient_name) && ~isstring(coefficient_name)
        error('VariableNameMapping_v004:InvalidCoefficientName', ...
            'coefficient_name must be a non-empty string or character array');
    end
    
    if ~istable(data_table)
        error('VariableNameMapping_v004:InvalidDataTable', ...
            'data_table must be a MATLAB table');
    end
    
    % Convert to string for consistent handling
    coefficient_name = char(coefficient_name);
    
    % Get mapping configuration
    config = getMappingConfig();
    
    % Handle intercept term
    if strcmp(coefficient_name, '(Intercept)')
        mapped_variable = 'intercept';
        return;
    end
    
    % Check for direct mapping first (exact variable name match)
    data_variables = data_table.Properties.VariableNames;
    if ismember(coefficient_name, data_variables)
        mapped_variable = coefficient_name;
        return;
    end
    
    % Handle interaction terms
    if contains(coefficient_name, ':')
        % For interaction terms, return the primary variable (first component)
        interaction_components = decomposeInteraction(coefficient_name, data_table);
        if ~isempty(interaction_components)
            mapped_variable = interaction_components{1};
            return;
        end
    end
    
    % Handle categorical variable encoding patterns
    mapped_variable = handleCategoricalEncoding(coefficient_name, data_table, config);
    
    % If still no mapping found, try pattern-based extraction
    if isempty(mapped_variable)
        mapped_variable = extractVariableFromCoeff(coefficient_name, data_table, config);
    end
    
    % Final validation
    if isempty(mapped_variable)
        % Create detailed error message with suggestions
        available_vars = strjoin(data_variables, ', ');
        error('VariableNameMapping_v004:UnmappableCoefficient', ...
            'Cannot map coefficient "%s" to data variable.\nAvailable variables: %s\nConsider checking categorical variable encoding.', ...
            coefficient_name, available_vars);
    end
    
    % Validate the mapped variable exists in data table
    if ~ismember(mapped_variable, data_variables)
        error('VariableNameMapping_v004:MappedVariableNotFound', ...
            'Mapped variable "%s" not found in data table. Available variables: %s', ...
            mapped_variable, strjoin(data_variables, ', '));
    end
end

function interaction_components = decomposeInteraction(interaction_term, data_table, varargin)
% Decompose interaction coefficient name into component variable names
%
% **ENHANCED DECOMPOSITION**: Handles complex interaction patterns with 
% categorical variable encoding and provides robust component mapping.
% **v003 FIX**: Replaced 'for i =' with 'for idx =' to avoid complex number collision.
%
% INPUTS:
%   interaction_term - String, interaction coefficient name (contains ':')
%   data_table - Table, the data table with actual variable names
%   varargin - Optional parameters
%
% OUTPUT:
%   interaction_components - Cell array of mapped component variable names

    % Input validation
    if nargin < 2
        error('VariableNameMapping_v004:InsufficientInputs', ...
            'decomposeInteraction requires interaction_term and data_table');
    end
    
    if ~contains(interaction_term, ':')
        % Not an interaction term - return as single component
        mapped_var = mapCoefficientToVariable(interaction_term, data_table);
        interaction_components = {mapped_var};
        return;
    end
    
    % Split interaction term
    raw_components = strsplit(interaction_term, ':');
    interaction_components = cell(size(raw_components));
    
    % Map each component - FIXED v003: Use 'idx' instead of 'i'
    for idx = 1:length(raw_components)
        component = strtrim(raw_components{idx});
        
        try
            % Try to map component to variable
            mapped_component = mapCoefficientToVariable(component, data_table);
            interaction_components{idx} = mapped_component;
            
        catch ME
            % If mapping fails, try to extract base variable name
            try
                % **FIX v002**: Provide required config parameter
                config = getMappingConfig();
                extracted_var = extractVariableFromCoeff(component, data_table, config);
                if ~isempty(extracted_var)
                    interaction_components{idx} = extracted_var;
                else
                    % As last resort, use the component as-is with warning
                    warning('VariableNameMapping_v004:UnmappableComponent', ...
                        'Could not map interaction component "%s". Using as-is.', component);
                    interaction_components{idx} = component;
                end
            catch
                % Final fallback
                interaction_components{idx} = component;
            end
        end
    end
    
    % Remove empty components
    interaction_components = interaction_components(~cellfun(@isempty, interaction_components));
    
    % Remove duplicates while preserving order
    [~, unique_idx] = unique(interaction_components, 'stable');
    interaction_components = interaction_components(unique_idx);
end

function isValid = validateVariableMapping(coefficient_name, data_table, varargin)
% Validate that a coefficient can be successfully mapped to a data variable
%
% **v002 ENHANCED VALIDATION**: FIXED to handle both categorical AND nominal variables
% This resolves the 47.1% validation failure rate from v001.
%
% INPUTS:
%   coefficient_name - String, coefficient name to validate
%   data_table - Table, the data table
%   varargin - Optional parameters
%
% OUTPUT:
%   isValid - Logical, true if mapping is possible

    isValid = false;
    
    try
        % Handle intercept term specially
        if strcmp(coefficient_name, '(Intercept)')
            % Intercept is always valid if the table is non-empty
            isValid = height(data_table) > 0;
            return;
        end
        
        % Attempt mapping
        mapped_var = mapCoefficientToVariable(coefficient_name, data_table);
        
        % Additional validation: check if variable contains sufficient data for segmentation
        if ismember(mapped_var, data_table.Properties.VariableNames)
            variable_data = data_table.(mapped_var);
            
            % v002 FIX: Check for reasonable variation in the variable
            % FIXED: Now handles both categorical AND nominal variables
            % CRITICAL FIX: Check nominal FIRST since nominal vars return true for both checks
            if isa(variable_data, 'nominal')
                % Handle nominal variables - FIXED v002
                try
                    % For nominal variables, check the number of unique levels
                    unique_codes = unique(double(variable_data));
                    % Remove NaN values (undefined nominal values)
                    unique_codes = unique_codes(~isnan(unique_codes));
                    isValid = length(unique_codes) >= 2;
                catch
                    % Fallback: if nominal processing fails, assume valid
                    isValid = true;
                end
                
            elseif iscategorical(variable_data)
                % Handle categorical variables (non-nominal)
                unique_categories = categories(variable_data);
                valid_categories = sum(~isundefined(unique_categories));
                isValid = valid_categories >= 2;  % Need at least 2 categories for segmentation
                
            elseif isnumeric(variable_data)
                unique_values = unique(variable_data(~isnan(variable_data)));
                isValid = length(unique_values) >= 2;  % Need at least 2 values for segmentation
                
            else
                isValid = true;  % Other data types assumed valid
            end
        end
        
    catch
        % If any error occurs during mapping, it's not valid
        isValid = false;
    end
end

%% SPECIALIZED MAPPING FUNCTIONS (UNCHANGED FROM v003)

function mapped_variable = handleCategoricalEncoding(coefficient_name, data_table, config)
% Handle categorical variable encoding patterns from LME models
%
% **CATEGORICAL ENCODING LOGIC**: MATLAB's fitlme function applies specific
% encoding patterns to categorical variables, which this function decodes.
% v002: Enhanced to handle both categorical and nominal variables.
% **v003 FIX**: Replaced 'for i =' with 'for idx =' to avoid complex number collision.
%
% INPUTS:
%   coefficient_name - String, encoded coefficient name
%   data_table - Table, original data table
%   config - Struct, mapping configuration
%
% OUTPUT:
%   mapped_variable - String, decoded variable name (empty if not found)

    mapped_variable = '';
    data_variables = data_table.Properties.VariableNames;
    
    % Get categorical variables from data table (including nominal)
    categorical_vars = {};
    for idx = 1:length(data_variables)  % FIXED v003: Use 'idx' instead of 'i'
        var_name = data_variables{idx};
        if iscategorical(data_table.(var_name)) || isa(data_table.(var_name), 'nominal')
            categorical_vars{end+1} = var_name;
        end
    end
    
    % Try matching against categorical variable patterns
    for idx = 1:length(categorical_vars)  % FIXED v003: Use 'idx' instead of 'i'
        cat_var = categorical_vars{idx};
        
        % Pattern 1: Variable name appears at start of coefficient
        if startsWith(coefficient_name, cat_var)
            mapped_variable = cat_var;
            return;
        end
        
        % Pattern 2: Variable name with suffix/encoding (e.g., "filterType_2")
        pattern = sprintf('^%s[_\\.]', cat_var);
        if ~isempty(regexp(coefficient_name, pattern, 'once'))
            mapped_variable = cat_var;
            return;
        end
        
        % Pattern 3: Check if coefficient contains the variable name
        if contains(coefficient_name, cat_var)
            % Additional validation: ensure it's not a substring of another variable
            is_unique_match = true;
            for j_idx = 1:length(data_variables)  % Use j_idx to avoid any confusion
                other_var = data_variables{j_idx};
                if ~strcmp(other_var, cat_var) && contains(other_var, cat_var) && contains(coefficient_name, other_var)
                    is_unique_match = false;
                    break;
                end
            end
            
            if is_unique_match
                mapped_variable = cat_var;
                return;
            end
        end
    end
    
    % Pattern 4: Handle common MATLAB categorical encodings
    % Remove common suffixes and prefixes that MATLAB adds
    cleaned_coeff = coefficient_name;
    
    % Remove numeric suffixes (e.g., "_1", "_2", ".1", ".2")
    cleaned_coeff = regexprep(cleaned_coeff, '[_\.]?\d+$', '');
    
    % Remove parenthetical suffixes (e.g., "(1)", "(2)")
    cleaned_coeff = regexprep(cleaned_coeff, '\(\d+\)$', '');
    
    % Try exact match with cleaned name
    if ismember(cleaned_coeff, data_variables)
        mapped_variable = cleaned_coeff;
        return;
    end
    
    % Try partial matches with cleaned name
    for idx = 1:length(data_variables)  % FIXED v003: Use 'idx' instead of 'i'
        var_name = data_variables{idx};
        if startsWith(var_name, cleaned_coeff) || startsWith(cleaned_coeff, var_name)
            mapped_variable = var_name;
            return;
        end
    end
end

function extracted_variable = extractVariableFromCoeff(coefficient_name, data_table, config)
% Extract base variable name from complex coefficient names
%
% **PATTERN-BASED EXTRACTION**: Uses configurable patterns to extract
% variable names from encoded coefficient names.
% **v003 FIX**: Replaced 'for i =' with 'for idx =' to avoid complex number collision.
%
% INPUTS:
%   coefficient_name - String, coefficient name to analyze
%   data_table - Table, data table for validation
%   config - Struct, extraction configuration
%
% OUTPUT:
%   extracted_variable - String, extracted variable name (empty if not found)

    extracted_variable = '';
    data_variables = data_table.Properties.VariableNames;
    
    % Apply extraction patterns from configuration
    extraction_patterns = config.extraction_patterns;
    
    for idx = 1:length(extraction_patterns)  % FIXED v003: Use 'idx' instead of 'i'
        pattern = extraction_patterns{idx};
        
        % Apply pattern to coefficient name
        tokens = regexp(coefficient_name, pattern, 'tokens');
        
        if ~isempty(tokens) && ~isempty(tokens{1})
            candidate_var = tokens{1}{1};
            
            % Validate extracted variable exists in data table
            if ismember(candidate_var, data_variables)
                extracted_variable = candidate_var;
                return;
            end
        end
    end
    
    % Fallback: try each data variable as a potential base
    for idx = 1:length(data_variables)  % FIXED v003: Use 'idx' instead of 'i'
        var_name = data_variables{idx};
        
        % Check if variable name is contained in coefficient name
        if contains(coefficient_name, var_name)
            % Additional validation to avoid false positives
            % Ensure it's not just a substring of another variable
            is_unique = true;
            for j_idx = 1:length(data_variables)  % Use j_idx to avoid any confusion
                other_var = data_variables{j_idx};
                if ~strcmp(other_var, var_name) && contains(other_var, var_name) && contains(coefficient_name, other_var)
                    is_unique = false;
                    break;
                end
            end
            
            if is_unique
                extracted_variable = var_name;
                return;
            end
        end
    end
end

%% CONFIGURATION AND UTILITY FUNCTIONS (UNCHANGED FROM v003)

function config = getMappingConfig(varargin)
% Get configuration for variable name mapping
%
% **CONFIGURABLE MAPPING**: Provides centralized configuration for all
% mapping operations, allowing customization for different data structures.
%
% OUTPUT:
%   config - Struct with mapping configuration

    config = struct();
    
    % Extraction patterns for different coefficient naming conventions
    config.extraction_patterns = {
        '^([a-zA-Z][a-zA-Z0-9_]*)[_\.].*',  % Variable name followed by underscore/dot
        '^([a-zA-Z][a-zA-Z0-9_]*)\(',       % Variable name followed by parenthesis
        '^([a-zA-Z][a-zA-Z0-9_]*)\d+$',     % Variable name followed by digits
        '([a-zA-Z][a-zA-Z0-9_]*)',          % Any valid variable name
    };
    
    % Common categorical variable encodings
    config.categorical_patterns = {
        '%s_\d+',      % Variable_Number
        '%s\.\d+',     % Variable.Number  
        '%s\(\d+\)',   % Variable(Number)
    };
    
    % Special handling for interaction terms
    config.interaction_separator = ':';
    config.handle_nested_interactions = true;
    
    % Error handling configuration
    config.strict_validation = true;
    config.allow_fallback_mapping = true;
    config.warn_on_ambiguous_mapping = true;
    
    % Apply any custom configuration from input arguments
    if nargin > 0 && isstruct(varargin{1})
        custom_config = varargin{1};
        field_names = fieldnames(custom_config);
        for idx = 1:length(field_names)  % FIXED v003: Use 'idx' instead of 'i'
            config.(field_names{idx}) = custom_config.(field_names{idx});
        end
    end
end

function mapping_table = createFullMappingTable(lme_model, data_table, varargin)
% Create comprehensive mapping table for all coefficients in an LME model
%
% **COMPREHENSIVE MAPPING**: Maps all coefficients to data variables with
% validation and diagnostics for troubleshooting.
% **v003 FIX**: Replaced 'for i =' with 'for idx =' to avoid complex number collision.
%
% INPUTS:
%   lme_model - LinearMixedModel, fitted LME model
%   data_table - Table, original data table
%   varargin - Optional parameters
%
% OUTPUT:
%   mapping_table - Table with columns: CoefficientName, MappedVariable, IsValid, Notes

    if nargin < 2
        error('VariableNameMapping_v004:InsufficientInputs', ...
            'createFullMappingTable requires lme_model and data_table');
    end
    
    % Extract coefficient names from model
    if isstruct(lme_model) && isfield(lme_model, 'Coefficients')
        coefficients = lme_model.Coefficients;
    elseif isa(lme_model, 'LinearMixedModel')
        coefficients = lme_model.Coefficients;
    else
        error('VariableNameMapping_v004:InvalidModel', ...
            'lme_model must be a LinearMixedModel or struct with Coefficients field');
    end
    
    coeff_names = coefficients.Name;
    n_coeffs = length(coeff_names);
    
    % Initialize mapping table
    mapped_variables = cell(n_coeffs, 1);
    is_valid = false(n_coeffs, 1);
    notes = cell(n_coeffs, 1);
    
    % Map each coefficient - FIXED v003: Use 'idx' instead of 'i'
    for idx = 1:n_coeffs
        coeff_name = coeff_names{idx};
        
        try
            mapped_var = mapCoefficientToVariable(coeff_name, data_table);
            mapped_variables{idx} = mapped_var;
            is_valid(idx) = validateVariableMapping(coeff_name, data_table);
            
            if is_valid(idx)
                notes{idx} = 'Successfully mapped and validated';
            else
                notes{idx} = 'Mapped but validation failed';
            end
            
        catch ME
            mapped_variables{idx} = '';
            is_valid(idx) = false;
            notes{idx} = sprintf('Mapping failed: %s', ME.message);
        end
    end
    
    % Create mapping table
    mapping_table = table(coeff_names, mapped_variables, is_valid, notes, ...
        'VariableNames', {'CoefficientName', 'MappedVariable', 'IsValid', 'Notes'});
    
    % Add summary information
    n_successful = sum(is_valid);
    n_failed = sum(~is_valid);
    
    fprintf('\n=== VARIABLE MAPPING SUMMARY (v004) ===\n');
    fprintf('Total coefficients: %d\n', n_coeffs);
    fprintf('Successfully mapped: %d (%.1f%%)\n', n_successful, 100*n_successful/n_coeffs);
    fprintf('Failed mappings: %d (%.1f%%)\n', n_failed, 100*n_failed/n_coeffs);
    
    if n_failed > 0
        fprintf('\nFailed mappings:\n');
        failed_idx = find(~is_valid);
        for idx = 1:length(failed_idx)  % FIXED v003: Use 'idx' instead of 'i'
            coeff_idx = failed_idx(idx);
            fprintf('  %s -> %s\n', coeff_names{coeff_idx}, notes{coeff_idx});
        end
    end
    
    fprintf('=======================================\n\n');
end

function debugMapping(coefficient_name, data_table, varargin)
% Debug variable mapping for troubleshooting
%
% **DEBUGGING UTILITY**: Provides detailed analysis of mapping process
% for troubleshooting coefficient naming issues.
% **v003 FIX**: Replaced 'for i =' with 'for idx =' to avoid complex number collision.
%
% INPUTS:
%   coefficient_name - String, coefficient to debug
%   data_table - Table, data table
%   varargin - Optional parameters

    fprintf('\n=== VARIABLE MAPPING DEBUG (v004) ===\n');
    fprintf('Coefficient: "%s"\n', coefficient_name);
    fprintf('Data table variables: %s\n', strjoin(data_table.Properties.VariableNames, ', '));
    
    % Check variable types in data table - FIXED v003: Use 'idx' instead of 'i'
    fprintf('Variable types in data table:\n');
    for idx = 1:length(data_table.Properties.VariableNames)
        var_name = data_table.Properties.VariableNames{idx};
        var_data = data_table.(var_name);
        if iscategorical(var_data)
            fprintf('  %s: categorical\n', var_name);
        elseif isa(var_data, 'nominal')
            fprintf('  %s: nominal\n', var_name);
        elseif isnumeric(var_data)
            fprintf('  %s: numeric\n', var_name);
        else
            fprintf('  %s: %s\n', var_name, class(var_data));
        end
    end
    
    % Test direct mapping
    data_variables = data_table.Properties.VariableNames;
    if ismember(coefficient_name, data_variables)
        fprintf('✓ Direct mapping found: %s\n', coefficient_name);
    else
        fprintf('✗ No direct mapping\n');
    end
    
    % Test interaction decomposition
    if contains(coefficient_name, ':')
        fprintf('→ Interaction term detected\n');
        try
            components = decomposeInteraction(coefficient_name, data_table);
            fprintf('  Components: %s\n', strjoin(components, ', '));
        catch ME
            fprintf('  Decomposition failed: %s\n', ME.message);
        end
    end
    
    % Test categorical encoding
    fprintf('→ Testing categorical encoding patterns...\n');
    config = getMappingConfig();
    try
        cat_mapped = handleCategoricalEncoding(coefficient_name, data_table, config);
        if ~isempty(cat_mapped)
            fprintf('  ✓ Categorical mapping: %s\n', cat_mapped);
        else
            fprintf('  ✗ No categorical mapping\n');
        end
    catch ME
        fprintf('  ✗ Categorical mapping error: %s\n', ME.message);
    end
    
    % Test pattern extraction
    fprintf('→ Testing pattern extraction...\n');
    try
        extracted = extractVariableFromCoeff(coefficient_name, data_table, config);
        if ~isempty(extracted)
            fprintf('  ✓ Pattern extraction: %s\n', extracted);
        else
            fprintf('  ✗ No pattern extraction\n');
        end
    catch ME
        fprintf('  ✗ Pattern extraction error: %s\n', ME.message);
    end
    
    % Final mapping attempt
    fprintf('→ Final mapping attempt...\n');
    try
        final_mapped = mapCoefficientToVariable(coefficient_name, data_table);
        fprintf('  ✓ Final mapping: %s\n', final_mapped);
        
        % Test validation
        is_valid = validateVariableMapping(coefficient_name, data_table);
        fprintf('  Validation: %s\n', ternary(is_valid, 'PASSED', 'FAILED'));
        
        % Additional validation details for v002
        if ismember(final_mapped, data_table.Properties.VariableNames)
            variable_data = data_table.(final_mapped);
            fprintf('  Variable type: %s\n', class(variable_data));
            
            if iscategorical(variable_data)
                unique_cats = categories(variable_data);
                fprintf('  Categories: %d (%s)\n', length(unique_cats), strjoin(unique_cats, ', '));
            elseif isa(variable_data, 'nominal')
                unique_vals = unique(variable_data);
                labels = getlabels(unique_vals);
                fprintf('  Nominal levels: %d (%s)\n', length(labels), strjoin(labels, ', '));
            elseif isnumeric(variable_data)
                unique_vals = unique(variable_data(~isnan(variable_data)));
                fprintf('  Unique numeric values: %d\n', length(unique_vals));
            end
        end
        
    catch ME
        fprintf('  ✗ Final mapping failed: %s\n', ME.message);
    end
    
    fprintf('=====================================\n\n');
end

%% UTILITY FUNCTIONS

function result = ternary(condition, true_val, false_val)
% Simple ternary operator for cleaner conditional assignment
    if condition
        result = true_val;
    else
        result = false_val;
    end
end