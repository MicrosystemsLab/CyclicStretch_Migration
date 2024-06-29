% ========================================================================%
% PIVdisplacement
% Calculates cell migration from PIVlab data. Returns average displacement
% of cells over time.
%
% Stacey Surace || Pruitt Lab
% June 25, 2024
%=========================================================================%
clear

load data % Load PIV data

global numImage timePerframe

% Constants
px2muRatio = 0.59;        % micron to pixel ratio
timePerframe = 5/60;      % time per image [hours]

% Conversion
um_hr = (px2muRatio/timePerframe);       % [microns/hour]

% PIV data grid
rows = size(uData,1);
columns = size(uData,2);
numImage = size(uData,3);

%% Migration Data

% Speed data
v_speedData = abs(vData);  % magnitude of y-direction velocity
u_speedData = abs(uData);  % magnitude of x-direction velocity

% Mean velocity/speed per frame (microns/hr)
% x-direction:
XvelocityMean = um_hr * frameMean(uData);
XspeedMean = um_hr * frameMean(u_speedData);
% y-direction:
YvelocityMean = um_hr * frameMean(vData);
YspeedMean = um_hr * frameMean(v_speedData);

% Displacement from mean velocity per frame (x and y directions)
XdispMean = dispMean(XvelocityMean);
YdispMean = dispMean(YvelocityMean);

% Total velocity magnitude (microns/hr) and total displacment (microns)
MagSpeedMean = um_hr * frameMean(velocitymagnitudeData);
TotDistMean = dispMean(MagSpeedMean);

% Mean direction of velocity per frame
DirectionMean = frameMean(Direction);                    % Degrees
% VortMean = (1/timePerframe) * frameMean(Vorticity);      % 1/hour

% All results
% (x Velocity, x speed, x displacement, y velocity, y speed,y displacement,
% total velocity, total displacement, direction)
AAVelDisp = [XvelocityMean' XspeedMean' XdispMean' YvelocityMean' ...
             YspeedMean' YdispMean' MagSpeedMean' TotDistMean' ...
             DirectionMean'];

%% Functions

function MeanFrame = frameMean(Data)
%-------------------------------------------------------------------------%
% Finds mean of data over each frame (no unit conversion)
% Input: field data
% Output: mean data matrix
%-------------------------------------------------------------------------%

global numImage

MeanFrame(1) = 0; % initially at rest at time=0

% through images
for image = 1:numImage

    MeanFrame(image+1) = mean(Data(:,:,image),'all');

end

end


function MeanDisp = dispMean(velMeanData)
%-------------------------------------------------------------------------%
% Finds displacement (microns) from reference point using mean velocity 
% (microns per hour) of each frame
% Input: mean velocity per frame
% Output: mean displacement matrix
%-------------------------------------------------------------------------%

global numImage timePerframe

MeanDisp(1) = 0;  % initial reference point at time=0

% Cycle through images
for image = 1:numImage

    % Change in position from mean velocity over time of frame
    MeanDisp(image+1) = MeanDisp(image) + timePerframe * velMeanData(image+1);

end

end
