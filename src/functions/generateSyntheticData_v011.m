function [xt, yt, k_local, v_local] = generateSyntheticData_v011(shapeNum, canvas, fs, powerLaw, yGain, orbitCount, resample, regen, debug)
% generateSyntheticData_v011 - Generate synthetic trajectory data for power law analysis
%
% Extended from v010 to include better shape index validation and error handling
% Fixed the issue with handling shape indices
%
% PARAMETERS:
%   shapeNum - Index (1-9) that selects a shape from {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'}
%   canvas - Canvas dimensions [width, height]
%   fs - Sampling frequency in Hz
%   powerLaw - Beta value to enforce (typically ~1/3)
%   yGain - Velocity gain factor
%   orbitCount - Number of orbits around the shape perimeter
%   resample - Resample factor (0 for no resampling)
%   regen - Force regeneration flag (1 to regenerate cached data)
%   debug - Enable debug output (0 for no debug)
%
% RETURNS:
%   xt, yt - x and y coordinates of generated trajectory
%   k_local - Curvature values
%   v_local - Velocity values
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% sampling frequency for data recovery
sampleFrequency = fs;
perimeterOrbits = orbitCount;

% Define shapes array with angular frequencies
Frequency = {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'};

% Validate shape index and map to angular frequency
if shapeNum > 0 && shapeNum <= length(Frequency)
    % Keep shapeNum as index for JCL_generateCurvesFunc which expects an index
    shapeAngFreq = shapeNum;
    
    % For debugging, print which shape we're using
    if debug
        fprintf('Using shape index %d which corresponds to angular frequency %s\n', ...
            shapeNum, Frequency{shapeNum});
    end
else
    error('Shape index %d is out of bounds for the shape array (valid range: 1-%d).', ...
        shapeNum, length(Frequency));
end

% beta forced, we might want to generate outwith the normal Beta for non Power Law compliant motion
beta0 = powerLaw; % [0-1] what beta are we going to enforce on the shape.

% RESAMPLE THE SHAPE
if resample
    resampleShape = 1;
    resampleFactor = resample; % (we will end up with ((original number of points - 1) * factor) +1 for overlap
else
    if debug
        fprintf('No resampling requested\n');
    end
    resampleShape = 0;
end

% vary accuracy - 0 use the values generated, 1 use rounded values
pixelAccurate = 0; % leave off? Tablet returns the centroid of the stylus tip, and not an integer pixel
debugShp = 0; % display the spline

%% generate path XY
rectXY = canvas; %wacomXY
offCentre = 0;

f0=@(nu) 2/3*( 1+nu.^2/2 )./(1+nu.^2+nu.^4/15);

for shapesNum = 1:length(Frequency)
    % calculate true Beta from H&S
    beta0calc(shapesNum)=f0(str2num(Frequency{shapesNum}));
end

%% check if this already exists, otherwise regenerate
saveName = ['baselineShp', num2str(shapeAngFreq),'_', num2str(sampleFrequency),'Hz.mat'];
if isfile(saveName) && ~regen
    if debug
        fprintf('Loading cached shape data from %s\n', saveName);
    end
    load(saveName)
else
    if debug
        fprintf('Generating new shape data (will be cached as %s)\n', saveName);
    end
    %%generate the pure shape as per H&S
    try
        [pathXYtrue, thetaCurve] = JCL_generateCurvesFunc(shapeAngFreq, rectXY, offCentre);
    catch ME
        error('Error in JCL_generateCurvesFunc: %s. Check shapeAngFreq value (%d).', ...
            ME.message, shapeAngFreq);
    end
    
    % functional implementation of JCL_generateCurves.... but the points on
    % this are limited to about 140 for the ellipse... 
    % consider MATLAB ellispoid for arbitrary points...

    pathXYtrue = unique(pathXYtrue','rows','stable');
    % note this needs to be ' to get the unique rows correct shape for
    % ensuing
    if 1
        x = pathXYtrue(:,1); % as spun up
        y = pathXYtrue(:,2); %
    else
        x = flip(pathXYtrue(:,1)); % the points are clockwise
        y = flip(pathXYtrue(:,2)); % we want them in anti clockwise order to match movement
    end

    %% all the shapes are mirror images across the x axis apart from the 3 ang freq triangle

    if resampleShape % https://uk.mathworks.com/help/curvefit/splines-in-the-plane.html
        %do we want to cut off the overlap points
        try
            splX = spline(1:length(x),x,1:1/resampleFactor:length(x));
            splY = spline(1:length(y),y,1:1/resampleFactor:length(y));
        catch ME
            error('Error in spline resampling: %s. Check x, y data and resampleFactor (%d).', ...
                ME.message, resampleFactor);
        end

        % this spline requires overlap to be effective
        if debugShp
            figure(1)
            clf
            plot(splX,splY,'k-')
            hold on
            plot(x,y,'r.')
            hold off
        end

        pathXYresample = [splX ; splY]';
    else
        pathXYresample = [x y];
    end

    pathXY = pathXYresample(); % this must be a full overlapping loop of points.

    %% explore two methods for recovering curvature from points
    for idx = 1:length(pathXY)
        if idx == 1 % then take the point from end
            pathTriplet = [pathXY(end-1,:); pathXY(idx,:); pathXY(idx+1,:)]; % account for END == 1
        elseif idx == length(pathXY)
            pathTriplet = [pathXY(end-1,:); pathXY(end,:); pathXY(1+1,:)]; % account for 1 == END
        else
            pathTriplet = [pathXY(idx-1,:); pathXY(idx,:); pathXY(idx+1,:)];
        end

        % method one
        [R,xcyc] = fit_circle_through_3_points(pathTriplet);

        centreLocal(idx,:) = xcyc';
        K(idx) = 1/R;  % method one
    end

    %method 2 MENGER CURVATURE!
    k = LineCurvature2D(pathXY);  %This method is simplest, and is equivalent to drawing a circle through 3 points.
    % is uses the next 2 values if it is the first value... and the previous 2
    % if it is the last...

    % methods are broadly equivalent
    try
        save(saveName, 'x', 'y', 'pathXY','pathXYresample' , "k", 'K', "splX", "splY");
    catch
        save(saveName, 'x', 'y', 'pathXY','pathXYresample' , "k", 'K');
    end
end

% returnShape = pathXY;
%
% figure(100)
% plot(K);
% hold on;
% plot(k)
% plot(abs(k)-K')
% title('2 METHODS OF CALCULATING CARTESIAN CURVATURE COMPARED')

k = K'; %we use the first method

%% k has a discontinuity error at the ends

k(end-(4*resample):end) = k(end-((4*resample)+1));
k(1:resample) = k(resample+1);

if debug
    figure(101)
    plot(k)
    title('check for discontinuities in k')
    drawnow
end

% verify we have data with the right Beta from the k via this formula
if debug
    log(yGain)
    v = ((yGain)* (k.^-beta0)); % we receive yGain as exp(yGain) - so flip it back
    figure(102)
    scatter(log(k),log(v));
    ylabel('LOG(v)')
    xlabel('LOG(k)')
    title('RAW RESAMPLED DATA')
    beta0
    %[betaR(1),yGainR(1),~] = regressDataEBR(log(v), log(k), 1, [], 1);
    [betaR(2),yGainR(2), ~] = regressDataEBR(log(v), log(k), 2, [], 1)
    [betaR(3),yGainR(3), ~] = regressDataEBR(log(v), log(k), 3, [log(yGain) beta0], 1)

    drawnow
    disp('WAIT FOR KEY PRESS')
    pause
end

%% we now inflate the resampled data to have several perimeter orbits
splX = repmat(splX,[1, perimeterOrbits]);
splY = repmat(splY,[1, perimeterOrbits]);
k = repmat(k,[perimeterOrbits,1]);

for idx = 2:length(splX)
    EucDistResample(idx-1) = ( (splX(idx)-splX(idx-1))^2 + (splY(idx)-splY(idx-1))^2 ) ^0.5;
end

perimeterEucResample = sum(EucDistResample);
% coherence check - this divided by 5 should be less than canvas x

%% so now we have a perimiter we can race around from the first point and find where we are

%% so we are at point 1 in resampled XY..
% we have curvature, Beta and samplingFrequency.. and thereforce
% displacement given the calc of speed.  This will tell us how heaviuly we
% need to resample for a speed at a freq for it not to be comedic.

% what speed should we be at beta0, at the first point on the list
% then where about on the list are we given the displacement of that speed
% given this frequency...
%
% we can then sub in the curvature differences and speed differences as
% white noise?  as similar to the envelope seen organically in the Shapes
% data....

xt(1) = splX(1);
yt(1) = splY(1);

cumulDistResample = cumsum(EucDistResample);
% since we start on the first point we ditch the first point of this as we
% are now targeting the second... need to generalise this for repeated
% loops... or just recall this whole function?

cumulDistResample = cumulDistResample - cumulDistResample(1);

% this will break down if the movement speed is greater than the length of
% the perimeter

% after we have generated we can test it still conforms to the power law.
% and then have observers check them...

%figure(200)
%clf
travelledArcLength = [];
thetaLatch = 0;
%% NEEDS LATCHING FOR TAN as we cross the x axis... different behaviour for overlapping shapes and simple hulls
idx = 0;

while 1
    idx = idx + 1;

    %% UPDATE K!  This uses segmented k from the target points... not analystical k.
    if idx == 1
        k_local(idx) = k(1);
    else
        k_local(idx) = k(whereOnPerim(1));
    end

    % determine current speed to satisfy Power Law
    currentSpeed = (yGain) * k_local(idx)^-beta0; % 7.5 pixels per second, for ellipse...
    v_local(idx) = currentSpeed;
    % which is ~ the EucDist for the initial points, and so an order of magnitude above when we resample at 10x

    %how far does that speed take us along the hull in 1 sample?
    travelledArcLength(idx) = currentSpeed / sampleFrequency;

    % the sum of the travelledArcLength is compared to the cumsum of the
    % EucDistResample.. this tell us where we are on the perimeter

    % cumsum is a list of monotonically increasing distances around the
    % perimeter
    whereOnPerim = find(cumulDistResample > sum(travelledArcLength));

    if isempty(whereOnPerim)
        break % we are out of perimeter!
    end

    try
        excessTravel = sum(travelledArcLength) - cumulDistResample(whereOnPerim(1)-1);
    catch
        keyboard
    end

    whereOnPerim(1);
    targetXY = [splX(whereOnPerim(1)) splY(whereOnPerim(1))];

    theta = atan( (splY(whereOnPerim(1)) - splY(whereOnPerim(1)-1)) ...
        / (splX(whereOnPerim(1)) - splX(whereOnPerim(1)-1))); %inverse sine of opposite / hypo

    thetaKeep(idx) = theta;

    thetaStr = 'DUMMY';

    % we now account for the travel that exceeds a resample point
    deltaX = cos(theta) * excessTravel; %travelledArcLength(idx);
    deltaY = sin(theta) * excessTravel; %travelledArcLength(idx);

    %% tricky trigonometry
    if theta >0 &&idx ~=1
        if (thetaKeep(idx-1)*thetaKeep(idx) < 0)
            thetaLatch = thetaLatch +1;
        end
    end
    if thetaLatch > 1 && shapeNum > 5
        thetaLatch = 0;
    end

    if mod(thetaLatch,2) ==1 %modulo 2 remainder ==1
        xt(idx+1) = splX(whereOnPerim(1)-1) - deltaX;
        yt(idx+1) = splY(whereOnPerim(1)-1) - deltaY;
    else %modulo 2 remainder ==0
        xt(idx+1) = splX(whereOnPerim(1)-1) + deltaX;
        yt(idx+1) = splY(whereOnPerim(1)-1) + deltaY;
    end

    if debug
        if 0
            %% investigate if the points we are trying to reach is somewhere adrift?
            figure(909)
            scatter(splX, splY, '.')
            hold on
            scatter(splX(whereOnPerim(1)-1), splY(whereOnPerim(1)-1), 'o')
            plot([ xt(idx) xt(idx+1)], [ yt(idx) yt(idx+1)] );
            scatter(splX(whereOnPerim(1)), splY(whereOnPerim(1)), 'x')
            title([' .. spline, o prior point, x target point, line the delta'])
            %             disp('PRESS ANY KEY')
            %             pause
        else
            % display only every 13 samples.. otherwise it grinds to a halt
            %    if idx==1 || (mod(idx, sampleFrequency)==0)
            if idx==1 || (mod(idx, 13)==0)
                figure(909)
                scatter(xt(idx+1), yt(idx+1),'o')
                % axis([-200 550 -300 300]);
                hold on
                plot(splX,splY,'.');
                title(['ThetaLatch is ', num2str(thetaLatch)])
                drawnow
            end
        end
    end
end

if debug
    log(yGain)
    %v = ((yGain) * (k.^-beta0)); % 
    figure(102)
    scatter(log(k_local),log(v_local));
    ylabel('LOG(v)')
    xlabel('LOG(k)')
    title('AFTER WE MAKE DATA WITH EXTRA POINTS')
    beta0
    %[betaR(1),yGainR(1),~] = regressData(log(v_local), log(k_local), 1, [], 1);
    [betaR(2),yGainR(2), ~] = regressDataEBR(log(v_local), log(k_local), 2, [], 1)
    [betaR(3),yGainR(3), ~] = regressDataEBR(log(v_local), log(k_local), 3, [log(yGain) beta0], 1)
    drawnow
end

% disp(['Mean diff between points on resampled shape hull - ', num2str(mean(diff(cumulDistResample)))]);
% %mean(diff(cumulDistResample)) % should be an order or two magnitude lower than
% disp(['mean travel between samples... above should be an order lower than -',num2str(mean(travelledArcLength) )]);
% %mean(travelledArcLength) %... to allow for accuracy in looping around
% %disp('the resample value of 1000 is needed for yGain = 1, but as gain and velocities go up, resample can drop')
