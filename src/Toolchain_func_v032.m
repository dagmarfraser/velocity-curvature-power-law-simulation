function [DATA, beta, VGF, genDuration, errorMadirolas, errorCurvature] = Toolchain_func_v032(parLoop, loopCount, cfg)
% Toolchain main file - v032 with MATLAB 2022a HPC compatibility 
% Dagmar Scott Fraser d.s.fraser@bham.ac.uk
% v000  November 2022
% v025 Updated to use differentiateKinematicsEBR 
% v026 Removed excessive debug output but kept essential progress information
% v027 Updated to use generateCustomNoise_v001 
% v028 Updated to use generateCustomNoise_v002 for zero noise handling
% v029 Fixed array dimension issue: pass column vectors to differentiateKinematicsEBR
% v030 Added noise validation and sanitization to prevent filter failures with colored noise
% v031 CRITICAL FIX: Updated to use generateCustomNoise_v003 (v002 was corrupted)
%      Removed excessive noise sanitization as v003 generates clean noise
%      FIXED: Multiple undefined variables and improved variable naming
%      - Line 231 (hashing): noiseMagnitude → noiseStdDev
%      - Lines 258-259 (noise generation): undefined noiseType → defined noiseType
%      - Line 485 (display output): undefined noiseType → defined noiseType
%      - Renamed confusing 'MaType' to clear 'noiseType' throughout
% v032 MATLAB 2022a HPC compatibility - function structure verified for HPC environments
%      No nested functions - already compatible with MATLAB 2022a HPC requirements
%      Updated version tracking for consistency with caller v056
%
% Note on beta values: Throughout this toolchain, beta values are consistently treated 
% as positive (~1/3). The power law relation is implemented as v = VGF × κ^(-β) where
% the negative sign is incorporated in the equation. This provides consistency with 
% the literature referring to the "one-third power law" while avoiding sign ambiguity.

% Convert inputs to the expected types
parLoop = double(parLoop);
loopCount = double(loopCount);

% Initialize critical variables
localPath = pwd;
versionTc = cfg.versionTc;
cfg.Serialise = parLoop;

% Ensure displayGraphs settings are available and properly typed
if isfield(cfg, 'display') && length(cfg.display) >= 3
    displayGraphs = double(cfg.display(1));
    displayNoisyData = double(cfg.display(2));
    displayRegress = double(cfg.display(3));
else
    displayGraphs = 0;
    displayNoisyData = 0;
    displayRegress = 0;
end

% Ensure critical cfg fields are the correct type
if ~isfield(cfg, 'TrialNum') || isempty(cfg.TrialNum)
    error('MATLAB:PowerLaw:MissingParameter', 'Trial number (cfg.TrialNum) is missing. This is a required parameter.');
else
    cfg.TrialNum = double(cfg.TrialNum);
end

%% Shape and Trajectory Generation
shapeStr = {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'};
canvas = cfg.canvas;
fs = cfg.fs;

f0=@(nu) 2/3*( 1+nu.^2/2 )./(1+nu.^2+nu.^4/15); % beta0 estimation given the shape.

LMSeeds = [];

for shapesNum = 1:length(shapeStr)
    % calculate true Beta from H&S
    beta0calc(shapesNum)=f0(str2num(shapeStr{shapesNum}));
end

% Handle potential type mismatches
shapeChoice = double(cfg.ShapeChoice);
fs = double(cfg.fs);
variantPowerLaw = double(cfg.variantPowerLaw); 
powerLaw = double(cfg.powerLaw);
yGain = double(cfg.yGain);
resample = double(cfg.resample);
noiseType = double(cfg.noiseType);  % Noise spectral type: 0=white, 1=pink, 2=brown, etc.
noiseStdDev = double(cfg.noiseStdDev);
filterType = double(cfg.filterType);
filterParams = cfg.filterParams;
% For struct fields, ensure they're also double
if isstruct(filterParams)
    fieldNames = fieldnames(filterParams);
    for i = 1:length(fieldNames)
        if isnumeric(filterParams.(fieldNames{i}))
            filterParams.(fieldNames{i}) = double(filterParams.(fieldNames{i}));
        end
    end
    % Convert struct to array format for differentiateKinematicsEBR compatibility
    if isfield(filterParams, 'order') && isfield(filterParams, 'width')
        % SG filters (cases 4 & 6)
        filterParams = [double(filterParams.order), double(filterParams.width)];
    elseif isfield(filterParams, 'order') && isfield(filterParams, 'cutoff') && isfield(filterParams, 'zerolag')
        % Butterworth filters (cases 2 & 3)
        filterParams = [double(filterParams.order), double(filterParams.cutoff), double(filterParams.zerolag)];
    else
        % Empty struct for case 1 (simple diff) - convert to empty array
        filterParams = [];
    end
end
regressType = double(cfg.regressType);
edgeClip = double(cfg.edgeClip);
curvatureChoice = double(cfg.curvatureChoice);
limitBreak = double(cfg.limitBreak);
displayGraphs = double(cfg.displayGraphs);
LMSeeds = double([yGain powerLaw]); % LM seeds
orbitCount = double(cfg.orbitCount);

try
    % Generate Synthetic Data!  given a shape, a power law, and a yGain! (and
    % fs and canvas size.
    forceRegen = 0;
    genDebug = 0;
    
    % Adjust shape index for the Frequency array in generateSyntheticData_v011
    % The array is {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'}
    shapeIndex = shapeChoice;
    
    % Ensure the shape index is within bounds of the array
    % The ellipse is at index 6 (value '2')
    if shapeIndex > 9 || shapeIndex < 1
        error('Shape index %d is out of bounds for the shape array (valid range: 1-9). Aborting.', shapeIndex);
    end
    
    % Check if we can safely access the shape string
    if shapeIndex > length(shapeStr)
        error('Shape index %d exceeds the length of shape string array (%d). Aborting.', shapeIndex, length(shapeStr));
    end
    
    % Calculate angular frequency from shape string for reference
    shapeAngFreq = str2num(shapeStr{shapeIndex});
    
    % Generate synthetic data using the improved v011 function
    [xt, yt, k_local, v_local] = generateSyntheticData_v011(shapeIndex, canvas, fs, powerLaw, yGain, orbitCount, resample, forceRegen, genDebug);

    genDuration = length(xt)/fs;

    %% reverse yt to get counter clockwise motion which replicates the experiment and makes the Polar Form work out trivially
    % this requires the data to be symmetrical about 0.
    yt = -yt;

    if displayGraphs
        figure(201)
        for idx = 1:fs/10:length(xt)
            scatter(xt(idx), yt(idx))
            hold on
            drawnow
        end
    end

    % deterministic up to this point
    % reproducible noise and Least Squares
    % ENHANCED: Hash-based unique seeding for maximum diversity
    % Create deterministic hash from parameter combination
    paramStr = sprintf('%.0f_%.3f_%.6f_%.3f_%.1f_%.3f_%d_%d', ...
        double(shapeChoice), double(fs), double(powerLaw), double(yGain), ...
        double(noiseType), double(noiseStdDev), double(filterType), double(regressType));
    
    % Simple deterministic hash function
    hashValue = 0;
    for i = 1:length(paramStr)
        hashValue = mod(hashValue * 31 + double(paramStr(i)), 2^31-1);
    end
    
    % Create unique seed for this parameter combo + trial
    uniqueSeed = mod(hashValue + cfg.TrialNum * 100000, 2^31-1);
    rng(uniqueSeed, "twister");
    
    % CRITICAL FIX: Use generateCustomNoise_v003 (the working version)
    % v002 was corrupted and causing "Arrays have incompatible sizes" errors
    xNoise = xt + generateCustomNoise_v003(length(xt), noiseType, noiseStdDev, fs)';
    yNoise = yt + generateCustomNoise_v003(length(yt), noiseType, noiseStdDev, fs)';

    % Ensure arrays are proper column vectors to prevent size incompatibility
    xNoise = xNoise(:);
    yNoise = yNoise(:);
    
    % Verify noise arrays have the expected length
    if length(xNoise) ~= length(xt) || length(yNoise) ~= length(yt)
        error('Noise arrays have incorrect length. Expected %d and %d, got %d and %d', ...
            length(xt), length(yt), length(xNoise), length(yNoise));
    end

    % SIMPLIFIED: Basic validation only (v003 generates clean noise)
    % Only check for obvious problems since v003 is validated
    if any(~isfinite(xNoise)) || any(~isfinite(yNoise))
        warning('Non-finite values detected in noise arrays - replacing with trajectory values');
        xNoise(~isfinite(xNoise)) = xt(~isfinite(xNoise));
        yNoise(~isfinite(yNoise)) = yt(~isfinite(yNoise));
    end

    if displayNoisyData
        if 0
            skipper = fs / 6;
            figure(202)
            for idx = 1:skipper:length(xNoise)
                plot(xt, yt, '.'); hold on
                scatter(xNoise(idx), yNoise(idx), 500, 'o')
                title([num2str(length(xNoise)),' datalength at ',num2str(cfg.fs),'Hz - frame num ', num2str(idx)]);
                hold off
                drawnow
            end
        else % maker a movie of the data
            % Calculate frame skip to achieve 30 fps
            frameSkip = round(fs / 30);
            numFrames = floor(length(xNoise) / frameSkip);

            % Create invisible figure with smaller size
            fig = figure('Visible', 'off', 'Position', [100, 100, 400, 300]);

            % Preallocate movie frames
            movieFrames = struct('cdata', cell(1, numFrames), 'colormap', cell(1, numFrames));

            % Set up the plot once
            plot(xt, yt, '-', 'Color', [0.7, 0.7, 0.7]);
            hold on;

            % Create trail objects to store all points
            h_trail = plot(NaN(1, length(xNoise)), NaN(1, length(yNoise)), 'o', ...
                'Color', [0, 0.7, 0.9], 'MarkerSize', 4, 'MarkerFaceColor', [0, 0.7, 0.9], 'MarkerEdgeColor', 'none');
            h_current = plot(NaN, NaN, 'o', 'Color', [0, 0.5, 1], 'MarkerSize', 8, 'MarkerFaceColor', [0, 0.5, 1]);
            title_obj = title('');
            axisFiddle = 50;
            % Set axis limits (adjust as needed)
            xlim([min(xt) - axisFiddle, max(xt) + axisFiddle]);
            ylim([min(yt) - axisFiddle, max(yt) + axisFiddle]);

            % Main loop for creating the animation
            tic
            for frameIdx = 1:numFrames
                disp(['frame ', num2str(frameIdx), ' of ', num2str(numFrames), ' after time ', num2str(toc)]);
                idx = frameIdx * frameSkip;

                % Update trail with all previous points
                set(h_trail, 'XData', xNoise(1:idx), 'YData', yNoise(1:idx));

                % Update current position
                set(h_current, 'XData', xNoise(idx), 'YData', yNoise(idx));

                % Update title
                set(title_obj, 'String', sprintf('%d pts at %d Hz - Time: %.2f s', length(xNoise), fs, idx / fs));

                % Capture the frame
                drawnow;
                movieFrames(frameIdx) = getframe(fig);
            end

            % Create and write the video
            movieName = ['trajectory', num2str(parLoop)];
            v = VideoWriter([movieName,'.mp4'], 'MPEG-4');
            v.FrameRate = 30;  % Set to 30 fps
            open(v);
            writeVideo(v, movieFrames);
            close(v);

            % Close the figure
            close(fig);

            disp('Movie creation complete. Output saved as trajectoryXXXXX.mp4');
        end
    end

    if cfg.MaticSpline
        t = (1/cfg.fs:1/cfg.fs:(length(xNoise)/cfg.fs));
        xNoise = spline(t, xNoise, t);
        yNoise = spline(t, yNoise, t);
    end

    %% Transform trajectory into velocity in x and y
    % filterType - chooses filter and differentiation
    %   1 MATLAB diff, scaled by the sample rate to give an approximation of the
    %   2 Nth Order Fp Hz Low pass filter, filtfilt for zero lag followed by Finite Differences
    %   3 Finite Differences followed by Nth Order Fp Hz Low pass filter filtfilt for zero lag
    %   4 Savitzky-Golay smoothing differential filter.

    % CRITICAL FIX: Pass column vectors (not transposed) to differentiateKinematicsEBR
    % The function expects column vectors for proper matrix operations in Butterworth cases
    [dx, dy] = differentiateKinematicsEBR(xNoise, yNoise, filterType, filterParams, cfg.fs);

    % Verify the arrays are properly shaped
    if isempty(dx) || isempty(dy)
        error('Differentiation returned empty arrays');
    end

    % Add detailed diagnostic output for troubleshooting
    if size(dx, 1) <= 2*edgeClip || size(dy, 1) <= 2*edgeClip
        fprintf('ARRAY SIZE DIAGNOSTIC:\\n');
        fprintf('  xNoise: %dx%d (length=%d)\\n', size(xNoise, 1), size(xNoise, 2), length(xNoise));
        fprintf('  yNoise: %dx%d (length=%d)\\n', size(yNoise, 1), size(yNoise, 2), length(yNoise));
        fprintf('  dx: %dx%d\\n', size(dx, 1), size(dx, 2));
        fprintf('  dy: %dx%d\\n', size(dy, 1), size(dy, 2));
        fprintf('  edgeClip: %d\\n', edgeClip);
        fprintf('  filterType: %d\\n', filterType);
        error('Not enough data points after edge clipping. dx size: %dx%d, dy size: %dx%d, edgeClip: %d', ...
            size(dx, 1), size(dx, 2), size(dy, 1), size(dy, 2), edgeClip);
    end

    velocityX = dx(edgeClip:end-edgeClip,2);
    velocityY = dy(edgeClip:end-edgeClip,2);

    velocity = ( ( velocityX.^2 + velocityY.^2 ) .^0.5 ); %

    accelerationX = dx(edgeClip:end-edgeClip,3);
    accelerationY = dy(edgeClip:end-edgeClip,3);

    if curvatureChoice == 1
        curvature = curvatureKinematicEBR(velocityX, velocityY, accelerationX, accelerationY);
    else
        curvature = zeros(length(accelerationY), 1);
        for ks = 1:length(accelerationY)
            tripletXY = [ xNoise(edgeClip+ks-1:edgeClip+ks+1) yNoise(edgeClip+ks-1:edgeClip+ks+1) ]';
            % we must take the absolute Menger Curvature to
            % match the Kinematic Curvature
            curvature(ks) = abs(curvatureMengerEBR(tripletXY));
        end
    end

    % Ensure responses and predictors are column vectors
    responses = velocity(:);
    predictors = curvature(:);

    % Check if there are valid data points
    validIdx = ~isnan(responses) & ~isnan(predictors) & ~isinf(responses) & ~isinf(predictors) & predictors > 0;
    
    if sum(validIdx) < 10
        error('Not enough valid data points for regression. Valid points: %d, total points: %d', ...
            sum(validIdx), length(responses));
    end
    
    responses = responses(validIdx);
    predictors = predictors(validIdx);

    % we couch this in a try catch as sometimes iterative
    % regressions do not reach a result within the limited
    % iterations permitted.
    [beta, VGF, stats] = regressDataEBR(responses, predictors, regressType, LMSeeds, displayRegress, limitBreak);
    DATA = true;

    %% we only have curvature for the trimmed trajectory for calculating Madirolas Error
    trajX = dx(edgeClip:end-edgeClip,1);
    trajY = dy(edgeClip:end-edgeClip,1);

    % Make sure shapeChoice is valid for baselineShapesFunc_v003
    [xyDelta, curvatureDelta, returnk, returnShape, divergedTrajCount] = baselineShapesFunc_v003([trajX trajY], curvature, shapeIndex, 0);

    errorMadirolas = nanmean(abs(xyDelta));
    errorCurvature = nanmean(curvatureDelta);
    
catch ME
    % Handle errors more gracefully with detailed information
    fprintf('ERROR in Toolchain_func_v032: %s\\n', ME.message);
    
    % Add diagnostic information for common error patterns
    if contains(ME.message, 'Unable to perform assignment') || contains(ME.message, 'Arrays have incompatible sizes')
        fprintf('ASSIGNMENT/ARRAY SIZE DIAGNOSTIC:\\n');
        if exist('dx', 'var')
            fprintf('  dx size: %dx%d\\n', size(dx));
        else
            fprintf('  dx: not created yet\\n');
        end
        if exist('dy', 'var')
            fprintf('  dy size: %dx%d\\n', size(dy));
        else
            fprintf('  dy: not created yet\\n');
        end
        if exist('xNoise', 'var') && exist('yNoise', 'var')
            fprintf('  xNoise: %dx%d (length=%d)\\n', size(xNoise, 1), size(xNoise, 2), length(xNoise));
            fprintf('  yNoise: %dx%d (length=%d)\\n', size(yNoise, 1), size(yNoise, 2), length(yNoise));
        end
        if exist('filterType', 'var')
            fprintf('  filterType: %d\\n', filterType);
        end
        if exist('filterParams', 'var')
            fprintf('  filterParams: %s\\n', mat2str(filterParams));
        end
        if exist('edgeClip', 'var')
            fprintf('  edgeClip: %d\\n', edgeClip);
        end
    end
    
    % Set default error return values
    DATA = false;
    beta = NaN;
    VGF = NaN;
    genDuration = NaN;
    errorMadirolas = NaN;
    errorCurvature = NaN;
    
    % Rethrow the error for upstream handling if configured to do so
    if isfield(cfg, 'rethrowErrors') && cfg.rethrowErrors
        rethrow(ME);
    end
end

padNum = num2str(cfg.Serialise,'%09.f');
% this assumes we are running from the source directory - and we have our
% standard file structure...
if cfg.saveAll
    save(fullfile([localPath(1:end-3) 'data' filesep 'processed' filesep 'synthetic' filesep versionTc '_', padNum, '_serial']));
end

% Get filter and regress type strings for output
filterStr = { 'NON' '_BW' 'BW_' 'S-G' 'NIL' 'SGS'};
regressStr = { 'REGR' 'FIT0' 'FITY' 'LMLS' 'IRLS' };
noiseStr = { 'WHITE' 'PINK'};

% Display detailed per-item progress (keep this useful output)
disp(['betaCalc ',num2str(beta),...
    ' betaSet ',num2str(powerLaw),...
    ' dur ', num2str(genDuration),...
    ' ', filterStr{filterType},'-Filt ',...
    regressStr{regressType}, ' regr ',...
    'trial ',num2str(cfg.TrialNum), ...
    ', 1/f^',num2str(noiseType),' noise, ',...
    num2str(noiseStdDev/cfg.pixelScale),'mm StdDev, ', ...
    num2str(errorMadirolas/cfg.pixelScale),'mm errMad ',...
    num2str(parLoop/loopCount * 100),'% done'...
    ]);

end
