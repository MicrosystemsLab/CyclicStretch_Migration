% CorrelationData:
% Finds spatial correlation length of speeds from velocity magnitude field.
%
% Stacey Surace
%=========================================================================%
clear; clc; clf; close all;

% Constants/Settings
px2muRatio = 0.59;     % microscopy: 0.59 um to px
pxGridPoint = 25;      % Pixels between grid points of speed field
threshold = 0.1;       % threshold for desired correlation

% Reports all files
matFiles = dir('*.mat');       % struct array of matlab files
fileNum = size(matFiles,1);

% Loads files in directory
for i=1:fileNum

    matFileName = matFiles(i).name    % display file name

    % Load file
    str(i,1) = string(matFileName);
    load(matFileName(1:end-4))

    % calculate spatial correlation length
    CorrLength(i,1)=correlationLength(velocitymagnitudeData,px2muRatio,pxGridPoint,threshold);
end

function corrLength = correlationLength(speedField,px2muRatio,pxGridPoint,threshold)
%-------------------------------------------------------------------------%
% Solves for spatial correlation length
% Input: velocity field from PIV data
%-------------------------------------------------------------------------%

% Position grid
xLim = size(speedField,1);
yLim = size(speedField,2);
[x,y] = ndgrid(1:xLim,1:yLim);

maxR = 80;

% Solves corr coefficient over radial distance r for each image, t
for t = 1:size(speedField,3) 

    aveSpeed = mean(speedField(:,:,t),"all"); % average speed of image

    for r = 1:maxR % cycle through radial distances

        % Reset values
        sumNumer = 0;  % sum of distances r
        sumDenom = 0;
        
        % cycle through grid points (i,j)
        for i = 1:xLim
            for j = 1:yLim  
                
                sumspeed = 0;    % reset average speed

                % distance of each point from reference (i,j)
                rMag = sqrt(abs((x - i)).^2 + abs((y - j)).^2);

                % Radial distances are binned in 14.75 microns sections
                % 1 r = 14.75 microns (for PIVlab data)
                % Find indices for bin at radius r
                [row, col] = find(rMag >= r & rMag < (r + 1));

                if isempty(row) == 0  % must have a speed element in bin

                    % add speeds within bin to solve for average
                    for cc = 1: size(row,1)
                        
                        sumspeed = sumspeed + speedField(row(cc),col(cc),t);

                    end

                    aveSpeedR = sumspeed/size(row,1);  % calculate average
                    
                    % Correlation coefficient equation
                    sumNumer = sumNumer + ((speedField(i,j,t) - aveSpeed) ...
                        * (aveSpeedR - aveSpeed));
                    sumDenom = sumDenom + sqrt((speedField(i,j,t) ...
                        - aveSpeed)^2 * (aveSpeedR - aveSpeed)^2);

                end

            end

        end
        
        Cr(r) = sumNumer/(sumDenom); % corr coefficient for each value r

    end

    Crt(t,:) = Cr(1,:); % collect coefficients over all time periods
end

MeanC = mean(Crt); % time-average corr coefficient

% remove any NAN values
cc = ~isnan(MeanC);
MeanC=MeanC(cc);

% Convert length r to microns
microns = px2muRatio .* pxGridPoint .* (1:size(MeanC,2));

%-------------------------------------------------------------------------%
% Plotting
%-------------------------------------------------------------------------%

% plot correlation coefficient over length
plot(microns, MeanC)
hold on
title("Spatial Correlation Coefficient of Migration Speeds")
xlabel("Distance (um)")
ylabel("Correlation Coefficient")

% Find length where corr coefficient intercepts threshold
corrLength = interp1(MeanC(1:maxR/2.5),microns(1:maxR/2.5), threshold)

% plot threshold line and mark intercept
plot(microns, threshold*ones(size(microns)), ':k')
plot(corrLength, threshold, 'xr', 'MarkerSize',10)

end


