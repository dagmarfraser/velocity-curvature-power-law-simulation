function [xyDelta, curvatureDelta] = baselineShapesFunc_v001(rawXY, rawCurvature, shapeNum, debug)

% inputs
% rawXY - WACOM tablet points.. generalise to Lydia's later
% rawCurvature = presplined curvature from the Frenet Serret calc
% shapeNum enumerates this list in the original {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'};
% in the shapes experiment    NaN NaN 2/33 2/5 4/5 4/3 2 3 4
% debug - show your working
% outputs
% curvature and xy delta.. xydelta can be absolute... curvature delta will
% be minus for smaller and + for larger, but does this need to be scaled by
% the actual local curvature?

% sampling frequency for data recovery
sampleFrequency = 133 ; %fixed in early version, will later generalise to Lydia's ~60
shapeAngFreq = shapeNum;% Ellipse is 6th in the OG Frequency List, 7th in ours DSF

looky = ceil(133/4);

%v1 limitation - just one loop?  nah the points are just points...

%%  determine tangential speed...
%
% v = γκ^−1/3
%
% find curvature of arbitrary point our pathXY
% find the circle described by the previous and next point
% K = 1/r
% then for arbitray y we have the power law

yGain = 1; % Gain 'constant'
% we imagined we might need this for nearvy curvature calculations.. but we
% will just search at the beginning and end of the movement for nearest
% 3 linked points and derive curvature from them

resampleShape = 1; %
resampleFactor = 10; % (we will end up with ((original number of points - 1) * factor) +1 for overlap
% vary accuracy - 0 use the values generated, 1 use rounded values


%% generate path XY
Frequency =   {'NaN' 'NaN' '2/33' '2/5' '4/5' '4/3' '2' '3' '4'};
rectXY = [1920;1080]; %wacomXY
offCentre = 0; %

f0=@(nu) -2/3*( 1+nu.^2/2 )./(1+nu.^2+nu.^4/15);

for shapesNum = 1:length(Frequency)
    % calculate true Beta from H&S
    beta0calc(shapesNum)=f0(str2num(Frequency{shapesNum}));

end

%% check if this already exists, otherwise regenerate, maybe at differing resample?
saveName = ['baselineShp', num2str(shapeAngFreq),'_', num2str(sampleFrequency),'Hz.mat'];
if isfile(saveName)
    load(saveName)
else
    %%generate the pure shape as per H&S
    % note the mistmatch is shapeAngFreq...
    % ellipse is 6 in OG list, 7 in our paper implementation... DSF
    [pathXYtrue, thetaCurve] = JCL_generateCurvesFunc(shapeAngFreq-1, rectXY, offCentre);
    % functional implementation of JCL_generateCurves.... but the points on
    % this are limited to about 140 for the ellipse... we want more..?
    % consider MATLAB ellispoid for arbitrary points...
    length(pathXYtrue)
    pathXYtrue = unique(pathXYtrue','rows','stable');
    % note this needs to be ' to get the uniqe rows correct

    x = flip(pathXYtrue(:,1)); % ADD 150 for NO 0 CROSSINGS
    y = flip(pathXYtrue(:,2)); % ADD 150 for NO 0 CROSSINGS


    if resampleShape % https://uk.mathworks.com/help/curvefit/splines-in-the-plane.html
        %do we want to cut off the overlap points? not for the spline, before, but
        %certainly after

        splX = spline(1:length(x),x,1:1/resampleFactor:length(x));
        splY = spline(1:length(y),y,1:1/resampleFactor:length(y));

        % this spline requires overlap to be effective

        figure(1)
        clf
        plot(splX,splY,'k-')
        hold on
        plot(x,y,'r.')
        hold off

        pathXYresample = [splX ; splY]';
        %how robust is this to other shapes? works with the 4 leaf clover
        % however we then have far too small changes for the 3 point fit through
        % for the circle... but we can use the old points for that by seeing where
        % we are as each of these points in on that circle...

        % else

    end

    pathXY = pathXYresample(); % this must be a full overlapping loop of points.
    for idx = 1:length(pathXY)
        if idx == 1 % then take the point from end
            pathTriplet = [pathXY(end-1,:); pathXY(idx,:); pathXY(idx+1,:)]; % account for END == 1
        elseif idx == length(pathXY)
            pathTriplet = [pathXY(end-1,:); pathXY(end,:); pathXY(1+1,:)]; % account for 1 == END
        else
            pathTriplet = [pathXY(idx-1,:); pathXY(idx,:); pathXY(idx+1,:)];
        end
        idx;
        [R,xcyc] = fit_circle_through_3_points(pathTriplet);

        centreLocal(idx,:) = xcyc';
        K(idx) = 1/R;  % method one

    end

    k = LineCurvature2D(pathXY);  %This method is simplest, and is equivalent to drawing a circle through 3 points.
    % is uses the next 2 values if it is the first value... and the previous 2
    % if it is the last... good enough.. but not perfect... i.e. it is not an
    % analytical solution.

    % we want to save
    save(saveName, 'pathXY','pathXYresample' , "k", 'K', "splX", "splY");

end
% figure(100)
% plot(K);
% hold on;
% plot(k)
% plot(abs(k)-K')
% title('2 METHODS OF CALCULATING CARTESIAN CURVATURE')

k = K'; %we use the first method

% we assume we start on a non crossover point
if shapeAngFreq > 5 % then there are no crossovers.. we can use simple method.
    [curvatureIndex, xyDelta] = dsearchn(pathXYresample,rawXY); %returns
    % the indices of the closest points in pathXYResample to the query
    % points in rawXY measured in Euclidean distance.
    % this is balanced within and without the hull.

    % we then use curvatureIndex to get the indice of k to compare with
    % rawCurvature
    curvatureDelta = k(curvatureIndex) - rawCurvature; % unscaled...
    % is this balanced within and without the hull?  if we add abs

    if debug
        figure(101)
        clf
        for idx = 1:length(rawXY)

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
    disp("COMPLEX SHAPE ROUTE")
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

    %% first find out where the first point is closest too then look 100 points ahead

    lookAheadStart = curvatureIndexNaive(1); % is our start point

    curvatureIndex = curvatureIndexNaive(1); %say point 1692 on the pathXYresample is closest
    xyDelta = xyDeltaNaive(1);

    if debug
        figure(100)
        clf
    end

    escapeXY = 100; % pixels if we are 100 pixels away we are diverging from the shape
    recoveryXY = 10;
    divergedTraj = 0;

    for idx = 2:length(rawXY)

        if (lookAheadStart+looky) > length(pathXYresample)
            lookAheadStart = 1;
        end

        if debug
            scatter(pathXYresample(:,1), pathXYresample(:,2),'.')
            hold on
            scatter(pathXYresample(lookAheadStart:lookAheadStart+looky,1),pathXYresample(lookAheadStart:lookAheadStart+looky,2),'*')
            scatter(rawXY(idx,1), rawXY(idx,2),'o')

        end

        [curvatureIndexLocal, xyDelta(idx)] = dsearchn(pathXYresample(lookAheadStart:lookAheadStart+looky,:),rawXY(idx,:));
        if debug
            title(['xyDelta ', num2str(xyDelta(idx)), ' for point ', num2str(idx), ' of ', num2str(length(rawXY))]);
            hold off
            drawnow
        end
        curvatureIndex(idx) = curvatureIndexLocal + lookAheadStart - 1;

        lookAheadStart = curvatureIndex(idx);

        %% this would be much easier with a circular buffer... consider https://uk.mathworks.com/matlabcentral/fileexchange/52411-circularbuffer

        if xyDelta(idx) > escapeXY
            divergedTraj = 1;
            % this likely means we started on an ambigious point near a
            % cross over so lets start again...

        end

    end


    curvatureDelta = k(curvatureIndex) - rawCurvature; % unscaled...

    if divergedTraj
        %the points during the first orbit, before it catches up to the first look ahead must be dismissed.
        truePoints = find(xyDelta <recoveryXY)
        xyDelta(1:truePoints(1)) = NaN;
        curvatureDelta(1:truePoints(1)) = NaN;

        if debug

            figure(101)
            plot(xyDelta)
            title(['First ', num2str(truePoints(1)), ' of ', num2str(length(xyDelta)), ' diverged.'])
        end
    end


end


% CAVEAT - this can have points that have nearest, points that precede
% points later... so doing it in a loop seems indicated..





