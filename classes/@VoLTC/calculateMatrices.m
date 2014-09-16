function obj = calculateMatrices(obj, debugPlots)
%tackle wavelength problem first...
%assume only 1 depth for now... to keep things simple

%formerly
%function [ obj.ACollection, obj.A1stCollection, obj.A2ndCollection ] = s3dVOLTCreateModel(lens, film, pSLocations)
% Create a volume of linear transformations for optics
%
%  [obj.ACollection obj.A1stCollection obj.A2ndCollection] = s3dVOLTCreateModel(lens, film, pSLocations)
%
% Creates a collection of linear transforms to summarize a lens. This
% particular function will only use 1 single depth, but multiple field
% heights.
%
% Inputs:
%   lens:  A lens class object
%   film:  A film class object that specifies the distance
%   pSLocations:
%     An n x 3 matrix that contains the locations of point sources for ray
%     tracing. A linear model will be computed for each one of these
%     locations.  These linear models will be placed in a collection.
%
% Returns:
%
% The 4-dimensions of each linear transformation convert the
% four-dimensional (column? row? vector) (xPos,yPos,ang1,ang2) into the
% corresponding set of values in theoutput light field. (Need to make sure
% these are x and y or y and x and what the angles represent)
%
% The matrices represent (field height, depth, wavelength)
%
% When we interogate the VOLT class for a linear transform at (x,y,z), we
% use length(x,y) for the field height, z for the depth, and angle(x,y) for
% the angle.  We begin with the linear transform at the right field height
% and depth, and we apply rotation matrices to the transform to get the
% proper angle corresponding to the (x,y,z) position.
%
%   obj.ACollection: a 4 x 4 x n matrix.  Each 4 x 4 layer contains the 4 x 4
% linear transform that summarizes light field transforms at that
% particular point.  (will become 4 x 4 x n x w)
%
%   obj.A1stCollection: The lens will be split into 2 layers: the first half
% leading to the aperture.  This matrix will be a 4 x 4 x n matrix that
% will be a collection of A matrices that summarizes the FIRST HALF of the
% lens only.
%
%   obj.A2ndCollection: a 4 x 4 x n matrix that contains a collection of A
%   matrices that summarizes the 2nd half of the lens.
%
% See also:
%    s3dVOLTRTOnePoint
%
% Example:
%
% AL Vistasoft, 2014
%
% TODO
%   We probably need a wavelength dimension
%   We may need a VOLT class
%   We could attach the interpolation function to the VOLT Class
%   A flag for turning on and off the graphs and speeding up
%   VOLT should get extended to

%% Argument checking here

if (ieNotDefined('debugPlots'))
    debugPlots = false;
end

%% Initialize complete matrices

depths = obj.get('depths');
fieldPositions = obj.get('fieldPositions');
wave = obj.get('wave');

obj.ACollection = zeros(4, 4, length(fieldPositions), length(depths), length(wave));
obj.A1stCollection = zeros(4, 4, length(fieldPositions), length(depths), length(wave));
obj.A2ndCollection = zeros(4, 4, length(fieldPositions), length(depths), length(wave));

%% Loop on points and make the various matrices

for depthIndex = 1:length(depths)
    pSLocations = obj.get('pSLocations', depthIndex);
    for pSIndex = 1:length(pSLocations)


        %% point sources (units are mm)
        pointSource = pSLocations(pSIndex, :);

        [ppsf, x, b, bMiddle] = s3dVOLTRTOnePoint(pointSource, obj.film, obj.lens);
        xFull = x;
        bFull = b;
        bMiddleFull = bMiddle;
        
        %% Do a speparate linear calculation for each wavelength.  
        % Split rays by wavelenghts and loop.
        
        for w = 1:length(wave)
            %%  We wonder about the full linear relationship
            %  b = Ax
            % To solve, we would compute
            % A = b\x

            % A = (x'\b')';
            % bEst = A * x;
            
            %only take data belonging to current wavelength
            
            waveIndex = ppsf.get('waveIndex');
            survivedWaveInd = waveIndex(ppsf.get('liveindices')); %remove nans to be on par with b and x
            
            %only use the x and b for the specific waveInd
            b = bFull(:, survivedWaveInd == w);  %this might be a cumbersome way to do this.  consider using bOrig
            x = xFull(:, survivedWaveInd == w);   %but see if this works first.  
            bMiddle = bMiddleFull(:, survivedWaveInd == w);
            
            A = b/x;
            bEst = A * x;

            % Scatter plot of positions
            for ii=1:4
                if(debugPlots)
                    vcNewGraphWin; plot(b(ii,:),bEst(ii,:),'o');
                    grid on;
                end
                
                meanAbsError = mean(abs(bEst(ii,:) - b(ii,:)));
                averageAmp = mean(abs(b(ii,:)));
                meanPercentError = meanAbsError/averageAmp * 100
            end

            obj.ACollection(:,:, pSIndex, depthIndex, w) = A;

            %% Calculate split A's: one for each half of the lens, divided by the middle aperture

            A1st = bMiddle/x;
            obj.A1stCollection(:,:, pSIndex, depthIndex, w) = A1st;
            bMiddleEst = A1st * x;

            A2nd = b/bMiddle;
            obj.A2ndCollection(:,:, pSIndex, depthIndex, w) = A2nd;

            %calculate final result
            bEst = A2nd * A1st * x;

            for ii=1:4
                ii
                if debugPlots
                    vcNewGraphWin; plot(b(ii,:),bEst(ii,:),'o');
                    grid on;
                end
                meanAbsError = mean(abs(bEst(ii,:) - b(ii,:)))
                averageAmpSplit = mean(abs(b(ii,:)));
                meanPercentErrorSplit = meanAbsError/averageAmpSplit * 100
            end

            close all;
        end
    end
end
obj.AMatricesUpdated = true;
end

