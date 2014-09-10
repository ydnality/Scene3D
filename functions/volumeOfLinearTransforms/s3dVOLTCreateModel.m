function [ AComplete, A1stComplete, A2ndComplete ] = s3dVOLTCreateModel(lens, film, pSLocations)
% Create a volume of linear transformations for optics
%
%  [AComplete A1stComplete A2ndComplete] = s3dVOLTCreateModel(lens, film, pSLocations)
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
%   AComplete: a 4 x 4 x n matrix.  Each 4 x 4 layer contains the 4 x 4
% linear transform that summarizes light field transforms at that
% particular point.  (will become 4 x 4 x n x w)
%
%   A1stComplete: The lens will be split into 2 layers: the first half
% leading to the aperture.  This matrix will be a 4 x 4 x n matrix that
% will be a collection of A matrices that summarizes the FIRST HALF of the
% lens only.
%
%   A2ndComplete: a 4 x 4 x n matrix that contains a collection of A
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



%% Initialize complete matrices
AComplete = zeros(4, 4, length(pSLocations));
A1stComplete = zeros(4, 4, length(pSLocations));
A2ndComplete = zeros(4, 4, length(pSLocations));

%% Loop on points and make the various matrices

for pSIndex = 1:length(pSLocations)
    
    
    %% point sources (units are mm)
    pointSource = pSLocations(pSIndex, :);
    
    [ppsf x b bMiddle] = s3dVOLTRTOnePoint(pointSource, film, lens);
    
    %%  We wonder about the full linear relationship
    %  b = Ax
    % To solve, we would compute
    % A = b\x
    
    % A = (x'\b')';
    % bEst = A * x;
    
    A = b/x;
    bEst = A * x;
    
    % Scatter plot of positions
    for ii=1:4
        vcNewGraphWin; plot(b(ii,:),bEst(ii,:),'o');
        grid on;
        
        meanAbsError = mean(abs(bEst(ii,:) - b(ii,:)));
        averageAmp = mean(abs(b(ii,:)));
        meanPercentError = meanAbsError/averageAmp * 100
    end
    
    AComplete(:,:,pSIndex) = A;
    
    %% Calculate split A's: one for each half of the lens, divided by the middle aperture
    
    A1st = bMiddle/x;
    A1stComplete(:,:, pSIndex) = A1st;
    bMiddleEst = A1st * x;
    
    A2nd = b/bMiddle;
    A2ndComplete(:,:, pSIndex) = A2nd;
    
    %calculate final result
    bEst = A2nd * A1st * x;
    
    for ii=1:4
        ii
        vcNewGraphWin; plot(b(ii,:),bEst(ii,:),'o');
        grid on;
        
        meanAbsError = mean(abs(bEst(ii,:) - b(ii,:)))
        averageAmpSplit = mean(abs(b(ii,:)));
        meanPercentErrorSplit = meanAbsError/averageAmpSplit * 100
    end
    
    close all;
end

end

