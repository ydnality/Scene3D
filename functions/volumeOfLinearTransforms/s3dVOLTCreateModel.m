function [ AComplete A1stComplete A2ndComplete ] = s3dVOLTCreateModel(lens, film, pSLocations)
%[ AComplete A1stComplete A2ndComplete ] = s3dVOLTCreateModel(lens, film, pSLocations)
%
% Creates a collection of linear transforms to summarize a lens.  
% This particular function will only use 1 single depth, but multiple field
% heights.
%
% pSLocations: a n x 3 matrix that contains the locations of point sources
% where raytracing will be performed, and a linear model will be computed
% for each one of these locations.  These linear models will be placed in a
% collection.  n represents the number of point sources.
%
% AComplete: a 4 x 4 x n matrix.  Each 4 x 4 layer contains the 4 x 4
% linear transform that summarizes light field transforms at that
% particular point.
%
% A1stComplete: The lens will be split into 2 layers: the first half
% leading to the aperture.  This matrix will be a 4 x 4 x n matrix that
% will be a collection of A matrices that summarizes the FIRST HALF of the
% lens only.
%
% A2ndComplete: a 4 x 4 x n matrix that contains a collection of A matrices
% that summarizes the 2nd half of the lens.


%initialize complete matrices
AComplete = zeros(4, 4, length(pSLocations));
A1stComplete = zeros(4, 4, length(pSLocations));
A2ndComplete = zeros(4, 4, length(pSLocations));

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

