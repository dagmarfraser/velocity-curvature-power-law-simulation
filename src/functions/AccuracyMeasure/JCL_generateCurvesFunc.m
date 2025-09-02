function [pathXY, thetaCurve]=JCL_generateCurvesFunc(shapeAngFreq, rectXY, offCentre)
% functional implementation of JCL_generateCurves
%% quick program to draw shapes as per the winTAB task file...
% d.s.fraser@bham.ac.uk April 2020
% for the https://jencooklab.com/
% assumes PsychoPhysics Toolbox is installed.
% v000 - add function input and outputs
%
% shapeAngFreq 1-9 for Frequency =   {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'}
% rectXY gives rect [0;0;X,Y]
% offCentre = 1 gives behaviours similar to WACOM in H&S trial
%
commandwindow

Curve_Type =  {'Spiral' 'Spiral' 'Sinusoidal' 'Sinusoidal' 'Sinusoidal' 'Sinusoidal' 'Sinusoidal' 'Sinusoidal' 'Sinusoidal' } 
Frequency =   {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'}
% 2/33 seems to be bugged - will find out - Dagmar

PsychDebugWindowConfiguration

white=[255;255;255];     black=[0;0;0];     red=[255;0;0];     blue=[0;0;255];     green=[0;255;0];  grey = [128; 128; 128];

background_color=grey; % choose background colour here
path_color = black;
widthPath = 1; % choose line width here


screens=Screen('Screens');
if max(screens)==1 || max(screens)==2 || max(screens)==3  % if multiple monitors
    screenNumber=1; % choose the main monitor
elseif max(screens)==0
    screenNumber=0;
else error('Too many monitors?')
end

Resolution=Screen('Resolution', screenNumber);

%% choose one of these lines to frame the shape - 
 rect = [0;0;rectXY(1);rectXY(2)];
%rect=[ 0; 0; Resolution.width; Resolution.height];
%rect=[ 0; 0; 1920; 1080]; %  1.7778 / 1 forcing the screen to be the same size as the WACOM no matter what
% 
% rect=[ 0; 0; 2000; 1200]; % 1.6667/1 forcing the screen to be the same size as the *new* SAMSUNG TABLET
% rect=[ 0; 0; 1000; 600]; % 1.6667/1 forcing the screen to be the same size as the *new* SAMSUNG TABLET
% rect=[ 0; 0; 1200; 780]; % 1.6667/1 forcing the screen to be the same size as the *new* SAMSUNG TABLET
%rect=[ 0; 0; 1024; 768]; % forcing the screen to be 4:3 RECT
if ~offCentre
    chosenCentre  = [ rect(3)/2 rect(4)/2 ];  
else
    chosenCentre = (rect(3:4)/2); % this the centre of the screen about which the Path coords are to be taken.. not OpenGL
end
%[windowPtr,rect]=Screen(?OpenWindow?,windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,numberOfBuffers] [,stereomode] [,multisample][,imagingmode][,specialFlags][,clientRect][,fbOverrideRect][,vrrParams=[]]);
[window,rect]=Screen(screenNumber,'OpenWindow', background_color, rect,[],2);

for iterations = shapeAngFreq:shapeAngFreq %1:length(Frequency)
    
    [pathXY,thetaCurve]=GENERATE_ALL_Curves_SIMPLE(Curve_Type{iterations},Frequency{iterations});
    length(pathXY)
    pathSave{iterations} = pathXY; % store path for each Frequency
    
    %plot(path(1,:),path(2,:)); drawnow

    Screen('DrawLines',window, pathXY, widthPath, path_color, chosenCentre);
    Screen(window,'Flip');
    
    imageArray = Screen('GetImage', window);
    imwrite(imageArray, ['Freq_', num2str(iterations),'.png'])
    pause(1)
    
    %T = array2table(pathXY');
    %Write the table to a CSV file
    %writetable(T,['FreqListRect_', num2str(iterations),'.csv']); %- no worky?

end

Screen('Close')
sca
%save('pathXY','pathXY') 

