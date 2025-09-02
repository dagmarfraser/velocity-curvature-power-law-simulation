function [xyDelta, curvatureDelta, returnk, returnShape, divergedTrajCount] = baselineShapesFunc_v003(rawXY, rawCurvature, shapeNum, debugInput)

% inputs
% rawXY - WACOM tablet points.. generalise to Lydia's later
% rawCurvature = presplined curvature from the Frenet Serret calc
% shapeNum enumerates this list in the original {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'};
% in the shapes experiment    NaN NaN 2/33 2/5 4/5 4/3 2 3 4
% debug - show your working
% outputs
% curvature and xy delta.. xydelta can be absolute... curvature delta will
% be minus for smaller (within hull) and + for larger (outside hull), but
% does this need to be scaled by the actual local curvature?
% returnk - the reported k at the local part of the original curve locally, in v3 spanned
% by the points rawXY n-1, n, n+1

xyDelta = NaN; %euclidean distance from nearest point
%xyHull = NaN; % within hull or without, xyDelta -outside, +inside for
%simple shapes...
curvatureDelta = NaN;
returnk = NaN;
divergedTrajCount= NaN;

% sampling frequency for data recovery
sampleFrequency = 133 ; %fixed in early version, will later generalise to Lydia's ~60
shapeAngFreq = shapeNum;% Ellipse is 6th in the OG Frequency List, 7th in ours DSF

debug = debugInput;

% v1 limitation - just one loop?  nah the points are just points...
% v2 hot pursuit mode around the shape for complex overlapping shapes
pursuitMode = 0; % enable limited look ahead for overlapping shapes
looky = ceil(133/4); % used in pursuit method
% v3 non naive curvature delta spanning adjacent points on the hull
v3Curvature = 1; %0 = nearest, 1 = 3 points
% in v2  curvature was compared just at the nearest point
% but as we resample higher that will approach an infinitely large
% circle..or a straight line
% so we actually want the local true curvature expressed between the n-1, n
% and n+1 rawXY points in a resampledXY triplet

%%  determine tangential speed...
%
% v = γκ^−1/3
%
% find curvature of arbitrary point our pathXY
% find the circle described by the previous and next point
% K = 1/r
% then for arbitray y we have the power law

resampleShape = 1; % 0 use limited and non splined shape points, 1 resampled and splined points.
resampleFactor = 10; % (we will end up with ((original number of points - 1) * factor) +1 for overlap

%% generate path XY
Frequency =   {'NaN' 'NaN' '2/33' '2/5' '4/5' '4/3' '2' '3' '4'};
rectXY = [1920;1080]; %wacomXY
offCentre = 0; %

f0=@(nu) -2/3*( 1+nu.^2/2 )./(1+nu.^2+nu.^4/15);

for shapesNum = 1:length(Frequency)
    % calculate true Beta from H&S
    beta0calc(shapesNum)=f0(str2num(Frequency{shapesNum}));

end

%% check if this already exists, otherwise regenerate
saveName = ['baselineShp', num2str(shapeAngFreq),'_', num2str(sampleFrequency),'Hz.mat'];
if isfile(saveName)
    load(saveName)
else
    %%generate the pure shape as per H&S
    % note the mistmatch is shapeAngFreq...
    % ellipse is 6 in OG list, 7 in our paper implementation... DSF
    [pathXYtrue, thetaCurve] = JCL_generateCurvesFunc(shapeAngFreq-1, rectXY, offCentre);
    % this is the original display function from Huh and Sejnowski

    % functional implementation of JCL_generateCurves.... but the points on
    % this are limited to about 140 for the ellipse... we want more..?
    % consider MATLAB ellispoid for arbitrary points...

    pathXYtrue = unique(pathXYtrue','rows','stable');
    % note this needs to be ' to get the uniqe rows correct shape for
    % ensuing
    if (shapeAngFreq > 5) || ~pursuitMode
        x = (pathXYtrue(:,1)); % 
        y = -(pathXYtrue(:,2)); % 
    else
        x = flip(pathXYtrue(:,1)); % the points are clockwise
        y = flip(pathXYtrue(:,2)); % we want them in anti clockwise order to match movement
    end
    %% all the shapes are mnirror images across the x axis apart from the 3


    if resampleShape % https://uk.mathworks.com/help/curvefit/splines-in-the-plane.html
        %do we want to cut off the overlap points

        splX = spline(1:length(x),x,1:1/resampleFactor:length(x));
        splY = spline(1:length(y),y,1:1/resampleFactor:length(y));

        % this spline requires overlap to be effective
        if debug
            figure(1)
            clf
            plot(splX,splY,'k-')
            hold on
            plot(x,y,'r.')
            hold off
        end

        pathXYresample = [splX ; splY]';

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

    %method 2
    k = LineCurvature2D(pathXY);  %This method is simplest, and is equivalent to drawing a circle through 3 points.
    % is uses the next 2 values if it is the first value... and the previous 2
    % if it is the last...

    % methods are broadly equivalent

    save(saveName, 'pathXY','pathXYresample' , "k", 'K', "splX", "splY");

end

returnShape = pathXY;

% figure(100)
% plot(K);
% hold on;
% plot(k)
% plot(abs(k)-K')
% title('2 METHODS OF CALCULATING CARTESIAN CURVATURE')

k = K'; %we use the first method

% we assume we start on a non crossover point
if (shapeAngFreq > 5) || ~pursuitMode
    % then there are no crossovers.. we can use simple method.
    T = delaunayn(pathXYresample);
    [curvatureIndex, xyDelta] = dsearchn(pathXYresample,T,rawXY); %returns
    % the indices of the closest points in pathXYResample to the query
    % points in rawXY measured in Euclidean distance.
    % this is balanced within and without the hull.
    outind = -1;
    [xyDeltaHull,~] = dsearchn(pathXYresample,T,rawXY, outind);
    outsideHull = find(xyDeltaHull == outind);    
    xyDelta(outsideHull) = -xyDelta(outsideHull); %set outsideHull to -value

    if ~v3Curvature
        % we then use curvatureIndex to get the indice of k to compare with
        % rawCurvature - however this is an approximation that is affected
        % by the resample rate somewhat arbitrarily
        curvatureDelta = k(curvatureIndex) - rawCurvature; % unscaled...
        % a deviation of 0.01 from 0.1 vs 0.01 from 0.2 is not addressed in this measure

        % further is this balanced within and without the hull?
        % if we add abs before we take the average?

        % this is an approximation... rather a good one ultimately... however
        % the below is more exacting.

        returnk = k(curvatureIndex);

    else
        %V3 Curvature using the dsearchn nearest triplet

        idT = 1;%special case - use the current and NEXT TWO
        % find the suggested curvatureIndex of pathXYresample for poin n-1, n
        % and n+1 - then find the local curvature od the true shape to compare
        % to the expressed one.

        prePoint    = pathXYresample(curvatureIndex(idT  ),:);
        onPoint     = pathXYresample(curvatureIndex(idT+1),:);
        postPoint   = pathXYresample(curvatureIndex(idT+2),:);

        kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
        k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
        curvatureDelta(idT) = k(idT) - rawCurvature(idT);

        returnk(idT) = k(idT);

        for idT = 2:(length(rawXY)-1) % noraml case, use either side
            % find the suggested curvatureIndex of pathXYresample for poin n-1, n
            % and n+1 - then find the local curvature od the true shape to compare
            % to the expressed one.

            prePoint    = pathXYresample(curvatureIndex(idT-1),:);
            onPoint     = pathXYresample(curvatureIndex(idT  ),:);
            postPoint   = pathXYresample(curvatureIndex(idT+1),:);

            kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
            k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
            curvatureDelta(idT) = k(idT) - rawCurvature(idT);

            returnk(idT) = k(idT);

        end

        idT = length(rawXY);% special case use CURRENT and previous 2.
        % find the suggested curvatureIndex of pathXYresample for poin n-1, n
        % and n+1 - then find the local curvature od the true shape to compare
        % to the expressed one.

        prePoint    = pathXYresample(curvatureIndex(idT-2),:);
        onPoint     = pathXYresample(curvatureIndex(idT-1),:);
        postPoint   = pathXYresample(curvatureIndex(idT  ),:);

        kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
        k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
        curvatureDelta(idT) = k(idT) - rawCurvature(idT);

        returnk(idT) = k(idT);

    end

    divergedTrajCount = 0;

    if debug
        figure(101)
        clf
        for idx = 1:10:length(rawXY)

            scatter(pathXYresample(:,1), pathXYresample(:,2),'.')
            hold on
            scatter(pathXYresample(curvatureIndex(idx),1), pathXYresample(curvatureIndex(idx),2),'*')
            scatter(rawXY(idx,1), rawXY(idx,2),'o')
            title(['xyDelta ', num2str(xyDelta(idx)), ' for point ', num2str(idx), ' of ', num2str(length(rawXY))]);
            hold off
            drawnow

        end
    end

else
    %disp("COMPLEX SHAPE ROUTE")

    %% we need to repeat this section until Divergence does not occur in the first few samples
    % of course wacky divergence can occur later for terrible trials... but
    % this should zero in on the shape at least initially

    [curvatureIndexNaive, xyDeltaNaive] = dsearchn(pathXYresample,rawXY); %returns
    % the indices of the closest points in pathXYResample to the query
    % points in rawXY measured in Euclidean distance.
    % this is balanced within and without the hull.

    % we then use curvatureIndex to get the indice of k to compare with
    % rawCurvature
    curvatureDeltaNaive = k(curvatureIndexNaive) - rawCurvature; % unscaled...

    % is this balanced within and without the hull?  if we add abs
    % we can just cut out elements where there is crossover?
    % so rawXY values within circles at the crossover points.
    % something something speed from the current nearest point
    % limits the range of pathXYresample.
    % likely requires finding the first points nearest point... and then
    % and using a section of the resampleXY... then when xy is closer to
    % the first point of the next segment we flip over...

    %% first find out where the first point is closest too then look 'looky' points ahead

    lookAheadStart = curvatureIndexNaive(1); % is our start point

    curvatureIndex(1) = curvatureIndexNaive(1); %say point 1692 on the pathXYresample is closest

    xyDelta = xyDeltaNaive(1);

    if debug
        figure(100)
        clf
    end

    escapeXY = 100; % pixels if we are 100 pixels away we are diverging from the shape
    recoveryXY = 10;
    divergedTraj = 0; % boolean, have we diverged within the first 100
    divergedTrajCount = 0; % how much extra lookahead do we utilise
    idx = 1;

    %% INSERT WHILE HERE
    % while idx <= length(rawXY)
    % if divergedTraj, resample for a better start point i.e. look ahead
    % start
    % then set divergedTraj to 0 and continue.  We probably want to just look ahead a few x look ahead
    % as that will probably have escape this local minima...

    while idx<length(rawXY)

        idx = idx + 1; % 2...n

        if divergedTraj % this can only happen in the first 100 values.. arbitrary

            divergedTrajCount = divergedTrajCount + 1;

            lookAheadStart = curvatureIndexNaive(1)+(looky*divergedTrajCount); % is our start point

            curvatureIndex(1) = lookAheadStart;

            idx = 2;

            divergedTraj =0; % we restart

        end

        if divergedTrajCount > 20
            disp(['DIVERGED TRAJECTORY COUNT > 20, = ', num2str(divergedTrajCount)])

            debug = 1;
        end

        if divergedTrajCount > (length(pathXYresample)/looky)
            disp('DIVERGED TRAJECTORY COUNT EXCEEDS LENGTH OF SHAPE')
            disp('likely erratic initial points!')

            figure(102)
            plot(rawXY(:,1),rawXY(:,2));hold on
            scatter(rawXY(1,1),rawXY(1,2), 'x')

            keyboard
        end
        if (lookAheadStart+looky) > length(pathXYresample)

            lookAheadStart = 1 + (looky*divergedTrajCount);
            % we are overlapping with the next loop.. so go past a little
            % however if we are going around again with looky*divTraj - we
            % bail
            if (lookAheadStart+looky) > length(pathXYresample);
                idx = Inf;
                break
            end

        end

        %% original loop start
        %for idx = 2:length(rawXY)

        [curvatureIndexLocal, xyDelta(idx)] = dsearchn(pathXYresample(lookAheadStart:lookAheadStart+looky,:),rawXY(idx,:));

        if debug

            scatter(pathXYresample(:,1), pathXYresample(:,2),'.')
            hold on
            scatter(pathXYresample(lookAheadStart:lookAheadStart+looky,1),pathXYresample(lookAheadStart:lookAheadStart+looky,2),'*')
            scatter(rawXY(idx,1), rawXY(idx,2),'o')

            %         end
            %
            %
            %         if debug

            title(['xyDelta ', num2str(xyDelta(idx)), ' for point ', num2str(idx), ' of ', num2str(length(rawXY))]);
            hold off
            drawnow

        end

        % the local curvatureIndexLocal is only between lookAheadStart and
        % LookAheadStart+looky so we need to add lookAheadStart to get the
        % true index in the master list

        curvatureIndex(idx) = curvatureIndexLocal + lookAheadStart - 1;
        %% update the values for next loop

        lookAheadStart = curvatureIndex(idx);

        %% this would be much easier with a circular buffer... consider https://uk.mathworks.com/matlabcentral/fileexchange/52411-circularbuffer

        if xyDelta(idx) > escapeXY
            if idx < 50
                divergedTraj = 1;
                %           disp('DIVERGENCE of 100 PIXELS at 50 SAMPLES?')
                % this likely means we started on an ambigious point near a
                % cross over so lets start again...
            end

        end

    end

    % if idx == Inf
    % outside the update loop - we do them en masse
    if idx == Inf
        debug = debugInput;
        %then we have a deranged xyDelta so we do not set curvatureDelta
        %below
    else
        if ~v3Curvature
            % we then use curvatureIndex to get the indice of k to compare with
            % rawCurvature - however this is naive
            curvatureDelta = k(curvatureIndex) - rawCurvature; % unscaled...
            % is this balanced within and without the hull?  if we add abs before
            % we take the average?

            % this is an approximation... rather a good one ultimately... however
            % the below is more exacting.

            returnk = k(curvatureIndex);

        else
            %V3 Curvature using the dsearchn nearest triplet

            idT = 1;%special case - use the current and NEXT TWO
            % find the suggested curvatureIndex of pathXYresample for poin n-1, n
            % and n+1 - then find the local curvature od the true shape to compare
            % to the expressed one.

            prePoint    = pathXYresample(curvatureIndex(idT  ),:);
            onPoint     = pathXYresample(curvatureIndex(idT+1),:);
            postPoint   = pathXYresample(curvatureIndex(idT+2),:);

            kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
            k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
            curvatureDelta(idT) = k(idT) - rawCurvature(idT);

            returnk(idT) = k(idT);

            for idT = 2:(length(rawXY)-1) % noraml case, use either side
                % find the suggested curvatureIndex of pathXYresample for poin n-1, n
                % and n+1 - then find the local curvature od the true shape to compare
                % to the expressed one.

                prePoint    = pathXYresample(curvatureIndex(idT-1),:);
                onPoint     = pathXYresample(curvatureIndex(idT  ),:);
                postPoint   = pathXYresample(curvatureIndex(idT+1),:);

                kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
                k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
                curvatureDelta(idT) = k(idT) - rawCurvature(idT);

                returnk(idT) = k(idT);

            end

            idT = length(rawXY);% special case use CURRENT and previous 2.
            % find the suggested curvatureIndex of pathXYresample for poin n-1, n
            % and n+1 - then find the local curvature od the true shape to compare
            % to the expressed one.

            prePoint    = pathXYresample(curvatureIndex(idT-2),:);
            onPoint     = pathXYresample(curvatureIndex(idT-1),:);
            postPoint   = pathXYresample(curvatureIndex(idT  ),:);

            kT      = LineCurvature2D([prePoint; onPoint; postPoint]);
            k(idT) = kT(2);% Line Curvature returns k at all three points, we need the middle one
            curvatureDelta(idT) = k(idT) - rawCurvature(idT);

            returnk(idT) = k(idT);

        end
    end

end






