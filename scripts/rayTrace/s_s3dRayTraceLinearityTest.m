%% ray-tracing for realistic lens - PPSF
%
% Lens element positions are all negative, and the final exit plan can be
% considered as z = 0
%
%  This uses Snell's law and a lens prescription to create a tray trace.
%  The script is too long, and we need to start writing functions so that
%  the length is shortened and the clarity increased.
%  We are only ray-tracing ideal point sources in order to extract out point
%  spread functions.
%
%  We are also experimenting with the plenoptic point spread function
%  (PPSF).  this is a very early experiment, where we are splitting the
%  calculation into 2 steps.  The first step traces the rays from a single
%  point in the scene towards the lens.  Ray-tracing is performed through
%  the lens, and out the back aperture.  At this point, the rays may be
%  saved as data.  Next, the rays are traced from the end of the lens to
%  the sensor (this process is reasonably efficient and doesn't take much
%  time).  Using this format, differrent sensor depths may be used to
%  access the PSF.
%
%  This specific script renders 2 PPSFs at a set field position , but with
%  different depths.  We attempt to interpolate a PPSF for a
%  position half-way between these 2, and check the results with the ground
%  truth.
%
% AL Vistalab, 2014
%%
s_initISET

%% film (sensor) properties
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = pbrtFilmObject('position', [0 0 60 ],'size', [10 10], 'wave', 400:10:700);

%% lens properties
% diffractionEnabled = false;
%turning on diffraction does NOT make sense yet since we have not modeled
%the transformation of uncertainty from the middle aperture to the end of the lens

%TODO: add diffraction into this somehow
%      Create a function lens.write() that makes a file from a lens
%      object.

%initialize and read multi-element lens from file
lensFileName = fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat');
nSamples = 101;
lens = lensMEObject('apertureSample', [nSamples nSamples], 'fileName', lensFileName);

%lens illustration - very small aperture in the middle
% lens.draw();

%% point sources (units are mm)

% Millimeters from last surface.  Always at least the lens thickness
% away.
pointSourceDepth = max(10000,-(lens.get('totaloffset')+1));
pointSources = [ 0 0 -pointSourceDepth];  %large distance test
pointSourceFieldHeight = 0;
% pointSources = [ 0 0 -60];  %short distance test


%% ray trace and save ppsf - Not sure camera should have pointSources

% Use the multi element lens and film and a point source.  Combine into
% a camera that calculates the point spread function.
ppsfCamera = ppsfCameraObject('lens', lens, 'film', film, 'pointSource', pointSources);
ppsf = ppsfCamera.estimatePPSF();

%% Record on film
ppsfCamera.recordOnFilm();

% Bring up the pointspread in an optics window
ppsfCamera.showFilm();
% There is now an ISET oi object we can review
% oi = vcGetObject('oi'); plotOI(oi,'illuminance mesh log');

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
cAEntranceXY = cAEntranceXY(:, survivedRays);
% vcNewGraphWin;
% r = lens.get('sradius',1); [x,y] = circlePoints([],r); plot(x,y,'.');
% hold on; plot(cAEntranceXY(1,:),cAEntranceXY(2,:),'o'); axis equal

% Matrix of directions at entrance pupil.  This is 3 x nExitRays
% Write this: lf = ppsf.get('entrance lf')
%
% This direction is an (x,y,z) vector
entDirMatrix = [cAEntranceXY(1, :) - ppsf.pointSourceLocation(1);
    cAEntranceXY(2, :) - ppsf.pointSourceLocation(2);
    ppsf.aEntranceInt.Z * ones(size(cAEntranceXY(1,:))) - ppsf.pointSourceLocation(3)];
entDirMatrix = normvec(entDirMatrix, 'dim', 1);

% Here we have (x,y,z) positions in the entrace aperture.
% We also have the first two entries of the unit length vector direction of
% the ray.  Maybe we want the two angles of that ray.
x = [cAEntranceXY(1,:);
    cAEntranceXY(2,:);
    ppsf.aEntranceInt.Z * ones(size(cAEntranceXY(1,:)));
    entDirMatrix(1, :);
    entDirMatrix(2,:)];
% We would like to have a look at this
% Here is a way to plot the rays at the entrance pupil.
% When the point is very far away, these will appear parallel.  WHen the
% point is close, they should be diverging.

% All the directions from a common point
vcNewGraphWin;
nPoints = size(x,2);
z = zeros(1,nPoints);
line([z;x(1,:)+x(4,:)],[z;x(2,:)+x(5,:)],[z;x(3,:)+entDirMatrix(3,:)])
view([1 58])
title(sprintf('%i distance %i samples',pointSourceDepth,nSamples))
% The directions from the actual aperture position
% vcNewGraphWin;
% line([x(1,:);x(1,:)+x(4,:)],[x(2,:);x(2,:)+x(5,:)],[x(3,:);x(3,:)+entDirMatrix(3,:)])
% view([1 58])


%%

% Compute exit lightfield
cAExitXY = ppsf.aExitInt.XY';
cAExitXY = cAExitXY(:, survivedRays);

exitDirMatrix = ppsf.aExitDir';
exitDirMatrix = exitDirMatrix(:, survivedRays);

b = [cAExitXY(1,:);
    cAExitXY(2,:);
    ppsf.aExitInt.Z * ones(size(cAExitXY(1,:)));
    exitDirMatrix(1, :);
    exitDirMatrix(2,:)];

% We wonder about the linear relationship
%  b = Ax
% To solve, we would compute
% b\x = A

%compute linear transformation
%SVD? pInverse?

%% Future development for modifying the rays.


% Make a second ppsf object
%     modifyRays = ppsfObject();
%
%     % Take the ppsfRays from the first object, copy the properties of the
%     % ppsfRays into real data, not just a pointer to the data.
%     modifyRays.makeDeepCopy(ppsfCamera.ppsfRays);
%
%     % Trace the rays lens to sensor
%     modifyRays.recordOnFilm(ppsfCamera.film);

%% Show the images

% vcNewGraphWin;
% imshow(film.image/ max(film.image(:)));

%         oi = oiCreate;
%         oi = initDefaultSpectrum(oi);
%         oi = oiSet(oi, 'wave', film.wave);
%         oi = oiSet(oi,'photons',film.image);
%
%
%         optics = oiGet(oi,'optics');
%         optics = opticsSet(optics,'focal length',lens.focalLength/1000);
%         optics = opticsSet(optics,'fnumber', lens.focalLength/(2*1));
%         oi = oiSet(oi,'optics',optics);
%         hfov = rad2deg(2*atan2(film.size(1)/2,lens.focalLength));
%         oi = oiSet(oi,'hfov', hfov);
%
%         temp = film.position;
%         filmDistance = temp(3);
%         oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
%         vcAddAndSelectObject(oi); oiWindow;


