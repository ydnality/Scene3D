function [ppsf, inLF, outLF, middleLF, inLFOrig, outLFOrig, ppsfCamera] = s3dVOLTRTOnePoint(pointSource, film, lens)
% This function performs the ray-trace for one point source, and
% returns light-field information in the form of the input, output,
% middle light-fields, the ppsf, and the ppsfCamera used for the
% calculation
%
% [ppsf x b bMiddle xOrig bOrig ppsfCamera] = s3dVOLTRTOnePoint(pointSource, film, lens)
%
%
% x: the light-field at the entrance of the lens (rays that don't make
% it through are removed).
% xOrig: the light-field at the entrance of the lens (rays that don't
% make it through are NOT removed).
% b: the light-field at the exit of the lens (rays that don't make it
% through are removed).
% bOrig: the light-field at the exit of the lens (rays that dont' make
% it through are NOT removed).
% bMiddle: the light-field that exists at the middle aperture location.
%
% TODO:  Point in the psCreate here ..
%
% Millimeters from last surface.  Always at least the lens thickness
% away.
%
%pointSourceDepth = 100;   % What is happening when 10,000?
%pointSourceDepth = max(pointSourceDepth,-(lens.get('totaloffset')+1));
%pointSourceFieldHeight = 0;
% pointSources = [ 0 0 -60];  %short distance test
%
% AL, VISTASOFT, 2014

%% ray trace and save ppsf - Not sure camera should have pointSources

% Use the multi element lens and film and a point source.  Combine into
% a camera that calculates the point spread function.
currentFilm = pbrtFilmC();
currentFilm.makeDeepCopy(film);
ppsfCamera = ppsfCameraC('lens', lens, 'film', currentFilm, 'pointSource', pointSource);

%%
nLines =  100;  % Draw the ray trace if nLines > 0
ppsf = ppsfCamera.estimatePPSF(nLines);

%% Record on film
ppsfCamera.recordOnFilm(nLines);

% Bring up the pointspread in an optics window
oi = ppsfCamera.oiCreate;
vcAddObject(oi); oiWindow;
plotOI(oi,'illuminance mesh log');

%% Calculate light field at the entrance pupil plane and exit pupil - estimate the linear transform

% do we need to put this in terms of angles? or are direction using
% cartesian coordinates good enough?
% put this into Ax = b form

% x will be a 4 x numSamples matrix containing the input lightfield
% b will be a 4 x numSamples matrix containing the output lightfield
% A will be the 4 x 4 least squares fit for this transformation

% Light field representation of the (x,y) positions in the entrance pupil
% of just those rays that actually make it to the film. There is a nan in
% the value when the ray does not make it.
%
cAEntranceXY = ppsf.aEntranceInt.XY';   % 2 x nSamples_in_aperture x nWave

% Eliminate nans
survivedRays = ~isnan(cAEntranceXY(1,:));
cAEntranceXYOrig = cAEntranceXY;
cAEntranceXY = cAEntranceXY(:, survivedRays);

% This is the effective aperture
vcNewGraphWin;
whichElement = 1;
r = lens.get('sdiameter',whichElement)/2; [inLF,y] = circlePoints([],r); plot(inLF,y,'.');
hold on; plot(cAEntranceXY(1,:),cAEntranceXY(2,:),'o'); axis equal
grid on
title(sprintf('Entrance points that make it to the exit'));

% Matrix of directions at entrance pupil.  This is 3 x nExitRays
% Write this: lf = ppsf.get('entrance lf')
%
% This direction is an (x,y,z) vector
entDirMatrix = ...
    [cAEntranceXY(1, :) - ppsf.pointSourceLocation(1);
    cAEntranceXY(2, :) - ppsf.pointSourceLocation(2);
    ppsf.aEntranceInt.Z * ones(size(cAEntranceXY(1,:))) - ppsf.pointSourceLocation(3)];
entDirMatrix = normvec(entDirMatrix, 'dim', 1);

%this is the original matrix WITH nans in it - we will take out the
%spurious rays later using the middle aperture.  We need these when we
%use the 2 linear models and middle aperture model.
entDirMatrixOrig = ...
    [cAEntranceXYOrig(1, :) - ppsf.pointSourceLocation(1);
    cAEntranceXYOrig(2, :) - ppsf.pointSourceLocation(2);
    ppsf.aEntranceInt.Z * ones(size(cAEntranceXYOrig(1,:))) - ppsf.pointSourceLocation(3)];
entDirMatrixOrig = normvec(entDirMatrixOrig, 'dim', 1);


% Here we have (x,y,z) positions in the entrace aperture.
% We also have the first two entries of the unit length vector direction of
% the ray.  Maybe we want the two angles of that ray.
% x = [cAEntranceXY(1,:);
%     cAEntranceXY(2,:);
%     ppsf.aEntranceInt.Z * ones(size(cAEntranceXY(1,:)));
%     entDirMatrix(1, :);
%     entDirMatrix(2,:)];

inLF = s3dVOLTLFFromPosDir(cAEntranceXY, entDirMatrix);

inLFOrig = s3dVOLTLFFromPosDir(cAEntranceXYOrig, entDirMatrixOrig);

%     inLF = [cAEntranceXY(1,:);
%         cAEntranceXY(2,:);
%         entDirMatrix(1, :);
%         entDirMatrix(2,:)];

%     inLFOrig = [cAEntranceXYOrig(1,:);
%         cAEntranceXYOrig(2,:);
%         entDirMatrixOrig(1, :);
%         entDirMatrixOrig(2,:)];

% Distribution of one of the angles
%     vcNewGraphWin;
%     hist(entDirMatrix(1,:),100);

%% All the ray directions from a common point
%     vcNewGraphWin;
%     nPoints = size(x,2);
%     if nPoints > 1000, s = randi(nPoints,[500,1]);
%     else              s = 1:nPoints;
%     end
%     z = zeros(1,length(s));
%     line([z;x(1,s)+x(3,s)],[z;x(2,s)+x(4,s)],[z;x(3,s)+entDirMatrix(3,s)])
%     view([20 58])
%     title(sprintf('%i distance\n%i samples\n%.1f aperture',pointSourceDepth,nSamples,apertureMiddleD))
%     set(gca,'xlim',[-10 10],'ylim',[-10,10],'zlim',[0 1.5]);

% The directions from the actual aperture position
% vcNewGraphWin;
% line([x(1,:);x(1,:)+x(4,:)],[x(2,:);x(2,:)+x(5,:)],[x(3,:);x(3,:)+entDirMatrix(3,:)])
% view([1 58])


%% Compute middle aperture lightfield

middleDir = ppsf.aMiddleDir';
middleDir = normvec(middleDir, 'p', 2, 'dim', 1);
middleDir = middleDir(:, survivedRays);

middleXY = ppsf.aMiddleInt.XY';
middleXY = middleXY(:, survivedRays);

middleLF = s3dVOLTLFFromPosDir(middleXY, middleDir);

%     middleLF =  [middleXY(1,:);
%         middleXY(2,:);
%         middleDir(1, :);
%         middleDir(2,:)];

%% Compute exit lightfield
cAExitXYOrig = ppsf.aExitInt.XY';
cAExitXY = cAExitXYOrig(:, survivedRays);

exitDirMatrixOrig = ppsf.aExitDir';
exitDirMatrix = exitDirMatrixOrig(:, survivedRays);

outLF = s3dVOLTLFFromPosDir(cAExitXY, exitDirMatrix);
outLFOrig = s3dVOLTLFFromPosDir(cAExitXYOrig, exitDirMatrixOrig);

%     outLF = [cAExitXY(1,:);
%         cAExitXY(2,:);
%         exitDirMatrix(1, :);
%         exitDirMatrix(2,:)];

%     outLFOrig = [cAExitXYOrig(1,:);
%         cAExitXYOrig(2,:);
%         exitDirMatrixOrig(1, :);
%         exitDirMatrixOrig(2,:)];

end
