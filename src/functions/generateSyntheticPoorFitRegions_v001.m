function synthetic_poor_regions = generateSyntheticPoorFitRegions_v001(data_table, config)
% GENERATESYNTHETICPOORFITREGIONS_V001 Create synthetic poor fit regions for Stage 3 testing
%
% **ENHANCED**: Now uses VariableNameMapping_v004 safe comparison framework
% to eliminate all nominal vs numeric comparison errors.
%
% **PURPOSE**: Generate artificial poor fit regions to test Stage 3 conditional analysis
% functionality without requiring genuine inadequacies from Stage 2 assessment.
%
% **SYNTHETIC REGION TYPES**:
%   - VGF boundary regions (high vs low VGF values)
%   - Noise magnitude extremes (high noise vs low noise)
%   - Filter type interactions (differential filter performance)
%   - Beta generation edge cases (extreme beta values)
%   - Cross-parameter interaction zones
%
% **USAGE**: For testing Stage 3 pipeline functionality only
%
% INPUT:
%   data_table - Stage 1 data table with camelCase variables
%   config - Configuration structure with minRegionSize
%
% OUTPUT:
%   synthetic_poor_regions - Cell array of synthetic poor fit region structures

synthetic_poor_regions = {};

% Validate input data structure
required_vars = {'deltaBeta', 'betaGenerated', 'VGF', 'noiseMagnitude', 'noiseColor', 'filterType'};
for i = 1:length(required_vars)
    if ~ismember(required_vars{i}, data_table.Properties.VariableNames)
        error('Missing required variable: %s', required_vars{i});
    end
end

n_total = height(data_table);
min_region_size = config.minRegionSize;

% BLUEBEAR MEMORY SAFETY: Set reasonable region size limits for 288GB system
max_region_size = min(200000, n_total * 0.3); % Max 30% of data or 200K obs
fprintf('=== SYNTHETIC POOR FIT REGIONS FOR TESTING ===\n');
fprintf('Synthetic data mode: Generating artificial poor fit regions for Stage 3 testing...\n');
fprintf('Generating synthetic poor fit regions for testing...\n');
fprintf('  Total observations: %d\n', n_total);
fprintf('  Minimum region size: %d\n', min_region_size);
fprintf('  Maximum region size (BlueBear): %d\n', max_region_size);
fprintf('  💾 BlueBear: 288GB system - using moderate region limits\n');

%% Region 1: High VGF Values (Boundary Effects)
vgf_values = data_table.VGF;
vgf_threshold = quantile(vgf_values, 0.75); % Top 25% of VGF values
high_vgf_indices = find(vgf_values >= vgf_threshold);

% BLUEBEAR SAFETY: Cap region size if too large
if length(high_vgf_indices) > max_region_size
    fprintf('  📊 Large VGF region (%d obs) - sampling %d for BlueBear\n', length(high_vgf_indices), max_region_size);
    high_vgf_indices = randsample(high_vgf_indices, max_region_size);
end

if length(high_vgf_indices) >= min_region_size
    region1 = struct();
    region1.description = 'VGF_HIGH_boundary_effects';
    region1.indices = high_vgf_indices;
    region1.n_obs = length(high_vgf_indices);
    region1.inadequacy_type = 'boundary_nonlinearity';
    region1.synthetic_rationale = 'High VGF values where linear assumptions break down';
    synthetic_poor_regions{end+1} = region1;
    
    fprintf('  ✓ Region 1: High VGF boundary (n=%d)\n', region1.n_obs);
end

%% Region 2: High Noise Magnitude (Measurement Breakdown)
noise_values = data_table.noiseMagnitude;
noise_threshold = quantile(noise_values, 0.8); % Top 20% of noise values
high_noise_indices = find(noise_values >= noise_threshold);

% BLUEBEAR SAFETY: Cap region size if too large
if length(high_noise_indices) > max_region_size
    fprintf('  📊 Large noise region (%d obs) - sampling %d for BlueBear\n', length(high_noise_indices), max_region_size);
    high_noise_indices = randsample(high_noise_indices, max_region_size);
end

if length(high_noise_indices) >= min_region_size
    region2 = struct();
    region2.description = 'noiseMagnitude_HIGH_breakdown';
    region2.indices = high_noise_indices;
    region2.n_obs = length(high_noise_indices);
    region2.inadequacy_type = 'measurement_saturation';
    region2.synthetic_rationale = 'High noise conditions causing measurement system breakdown';
    synthetic_poor_regions{end+1} = region2;
    
    fprintf('  ✓ Region 2: High noise breakdown (n=%d)\n', region2.n_obs);
end

%% Region 3: Filter Type 1 with High Noise (Interaction Effect)
% ENHANCED SAFE COMPARISON: Use VariableNameMapping_v004 for maximum robustness
fprintf('  🔧 Filtering filterType==1 using safe comparison framework...\n');

% Try VariableNameMapping_v004 safe comparison first (most robust)
if exist('VariableNameMapping_v004.m', 'file')
    try
        filter1_indices = VariableNameMapping_v004('safeComparison', data_table.filterType, 1, 'equal');
        fprintf('  ✓ Safe comparison successful: Found %d filterType==1 observations\n', length(filter1_indices));
    catch ME
        fprintf('  ⚠ Safe comparison failed: %s\n', ME.message);
        fprintf('  🔧 Falling back to direct nominal handling...\n');
        filter1_indices = [];
    end
else
    fprintf('  ⚠ VariableNameMapping_v004 not found, using fallback approach\n');
    filter1_indices = [];
end

% Fallback to direct nominal/categorical handling if safe comparison fails
if isempty(filter1_indices)
    try
        if isa(data_table.filterType, 'nominal') || isa(data_table.filterType, 'categorical')
            filter1_indices = find(data_table.filterType == nominal(1));
            fprintf('  ✓ Fallback nominal comparison successful: Found %d observations\n', length(filter1_indices));
        else
            filter1_indices = find(data_table.filterType == 1);
            fprintf('  ✓ Fallback numeric comparison successful: Found %d observations\n', length(filter1_indices));
        end
    catch ME
        fprintf('  ❌ All comparison methods failed: %s\n', ME.message);
        fprintf('  🔧 Skipping Region 3 due to filterType comparison issues\n');
        filter1_indices = [];
    end
end

% Create intersection with high noise indices if we successfully got filter1_indices
if ~isempty(filter1_indices)
    filter1_high_noise = intersect(filter1_indices, high_noise_indices);
    
    if length(filter1_high_noise) >= min_region_size
        region3 = struct();
        region3.description = 'filterType_1_high_noise_interaction';
        region3.indices = filter1_high_noise;
        region3.n_obs = length(filter1_high_noise);
        region3.inadequacy_type = 'filter_noise_interaction';
        region3.synthetic_rationale = 'Filter Type 1 performance degradation under high noise';
        synthetic_poor_regions{end+1} = region3;
        
        fprintf('  ✓ Region 3: Filter 1 + high noise (n=%d)\n', region3.n_obs);
    else
        fprintf('  ⚠ Region 3: Insufficient observations (%d < %d required)\n', length(filter1_high_noise), min_region_size);
    end
else
    fprintf('  ⚠ Region 3: Skipped due to filterType filtering failure\n');
end

%% Region 4: Extreme Beta Values (Edge Case Performance)
beta_values = data_table.betaGenerated;
extreme_beta_indices = find(beta_values >= quantile(beta_values, 0.9) | beta_values <= quantile(beta_values, 0.1));

if length(extreme_beta_indices) >= min_region_size
    region4 = struct();
    region4.description = 'betaGenerated_extreme_values';
    region4.indices = extreme_beta_indices;
    region4.n_obs = length(extreme_beta_indices);
    region4.inadequacy_type = 'parameter_edge_effects';
    region4.synthetic_rationale = 'Extreme beta values where model assumptions become questionable';
    synthetic_poor_regions{end+1} = region4;
    
    fprintf('  ✓ Region 4: Extreme beta values (n=%d)\n', region4.n_obs);
end

%% Region 5: Low VGF with High Noise Color (Cross-Parameter Stress)
low_vgf_indices = find(vgf_values <= quantile(vgf_values, 0.25));
high_noise_color_indices = find(data_table.noiseColor >= quantile(data_table.noiseColor, 0.8));
stress_condition_indices = intersect(low_vgf_indices, high_noise_color_indices);

if length(stress_condition_indices) >= min_region_size
    region5 = struct();
    region5.description = 'low_VGF_high_noiseColor_stress';
    region5.indices = stress_condition_indices;
    region5.n_obs = length(stress_condition_indices);
    region5.inadequacy_type = 'cross_parameter_stress';
    region5.synthetic_rationale = 'Low VGF with high noise color creating measurement stress';
    synthetic_poor_regions{end+1} = region5;
    
    fprintf('  ✓ Region 5: VGF-noise color stress (n=%d)\n', region5.n_obs);
end

%% OPTIONAL: Test additional filterType values for robustness
% Only if we have VariableNameMapping_v004 available and filterType filtering worked
if exist('VariableNameMapping_v004.m', 'file') && ~isempty(filter1_indices)
    try
        % Test filtering for filterType != 1 to create additional diverse regions
        filter_not_1_indices = VariableNameMapping_v004('safeComparison', data_table.filterType, 1, 'not_equal');
        
        if length(filter_not_1_indices) >= min_region_size
            region6 = struct();
            region6.description = 'filterType_NOT_1_alternative_methods';
            region6.indices = filter_not_1_indices;
            region6.n_obs = length(filter_not_1_indices);
            region6.inadequacy_type = 'alternative_filter_performance';
            region6.synthetic_rationale = 'Non-primary filter types may show different adequacy patterns';
            synthetic_poor_regions{end+1} = region6;
            
            fprintf('  ✓ Region 6: Alternative filter types (n=%d)\n', region6.n_obs);
        end
    catch ME
        % Silently continue if additional testing fails
        fprintf('  ℹ Additional filterType testing skipped: %s\n', ME.message);
    end
end

%% Add synthetic inadequacy metrics to each region
for i = 1:length(synthetic_poor_regions)
    region = synthetic_poor_regions{i};
    
    % Simulate Stage 2 adequacy metrics that would trigger conditional analysis
    region.residual_pattern_cohens_d = 0.6 + 0.4 * rand(); % Above 0.5 threshold
    region.coefficient_instability = 0.04 + 0.02 * rand(); % Above 0.03 threshold
    region.cv_performance_degradation = 0.18 + 0.05 * rand(); % Above 15% threshold
    
    % Add fake adequacy assessment details
    region.adequacy_failed_criteria = {'residual_patterns', 'coefficient_stability'};
    region.global_model_r_squared_in_region = 0.7 + 0.2 * rand(); % Moderate fit
    region.expected_improvement_potential = 0.05 + 0.1 * rand(); % 5-15% improvement potential
    
    synthetic_poor_regions{i} = region;
end

% BLUEBEAR SAFETY: Limit total number of regions for efficient processing
max_regions = 4; % Reasonable for BlueBear's 288GB - allows thorough testing without excess
if length(synthetic_poor_regions) > max_regions
    fprintf('  🔧 Limiting synthetic regions from %d to %d for BlueBear efficiency\n', length(synthetic_poor_regions), max_regions);
    % Keep the most diverse regions (first few types)
    synthetic_poor_regions = synthetic_poor_regions(1:max_regions);
end

fprintf('Generated %d synthetic poor fit regions for Stage 3 testing\n', length(synthetic_poor_regions));
fprintf('  💾 BlueBear optimized: Moderate limits for 288GB system\n');
fprintf('  🔧 Safe comparison framework: Enhanced robustness for variable filtering\n');

% Final validation of all generated regions
fprintf('  🔍 Validating generated regions...\n');
for i = 1:length(synthetic_poor_regions)
    region = synthetic_poor_regions{i};
    
    % Basic validation
    if any(region.indices < 1) || any(region.indices > n_total)
        warning('Region %d (%s) has indices outside valid range', i, region.description);
    end
    
    % Check for duplicate indices
    if length(unique(region.indices)) ~= length(region.indices)
        warning('Region %d (%s) has duplicate indices - cleaning up', i, region.description);
        region.indices = unique(region.indices);
        region.n_obs = length(region.indices);
        synthetic_poor_regions{i} = region;
    end
end

fprintf('  ✅ All regions validated and ready for Stage 3 testing\n');

end