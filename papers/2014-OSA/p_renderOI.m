%% p_renderOI
%
% Illustrates a point rendering for a simple lens 
%
% The point spread is calculated using lensMEObject and ppsfCamera object.
% The results are shown in the OI window.
% The method is illustrated
% Some ray traces are put up.
%
% Shortening of p_renderOiMatlabToolFull.m


%%
s_initISET

%% If modded pbrt is NOT installed on this system, run this command to
% load a scene file
% sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');
% scene = load(sceneFileName);
% scene = scene.scene;
% vcAddObject(scene); sceneWindow;

%% Render a PSF collection (see p_Figure1.m also)
%
% Rendering the collection will become a function at some point
%
% To render a high quality PSF, we first get an approximation at coarse
% scale. Then we calculate the centroid of that PSF,  zoom in, and
% calculate a "high quality" PSF using a much smaller film but higher spatial resolution.
%
% We need to keep track of spatial units on the film as we do this.
%
% The high-resolution PSFs will be added into a PSF Collection matrix,
% which contains the whole series of PSFs depending on field height,
% wavelength, depth.
%
% The PSF collection is used to produce the forward calculation
% rendered image.

% We will loop through the point positions.  Units are millimeters
pX = 0:-1000:-3000;
pY = 0;                % Assume radial symmetry, so only calculate X
pZ =[-70 -80 -90 -110];% Depth

% Assumed reference Z point.
% All other points will use the same field angle as this reference Z point
normalizingZ = -16000;           % mm
[X, Y, Z] = meshgrid(pX,pY,pZ);

%adjust for approximate difference in field position when Z changes
for i = 1:length(pZ)
    X(:,:,i) = X(:,:,i) *  pZ(i)/normalizingZ;
end
pointSources = [X(:), Y(:), Z(:)];

numDepths = length(pZ);
numFieldHeights = length(pX) * length(pY);
psfsPerDepth = size(pointSources, 1)/numDepths;
jitterFlag = true;   % Enable jitter for lens front element aperture samples
nLines = false;      % Number of lines to draw for debug illustrations.

%%  Declare film properties for PSF recording.

% wave = 500;
wave = 400:100:700;            % Wavelength
wList = 1:length(wave);
fX = 0; fY = 0; fZ = 135.5;    % mm.  Film position

% Film width and height for coarse calculation
fW = 80;  % mm
fH = 80;  % mm

% Film resolution (preview)
numPixelsW = 151;
numPixelsH = 151;

% Film resolution (final render)
numPixelsWHQ = 75;
numPixelsHHQ = 75;

% The film width for high quality is set manually - this
% should be automated in the future
newWidth = 10;    % mm

%% Describe the lens

lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
load(lensFile,'lens')
lens.apertureMiddleD = 5;

% Preview and high quality sampling on first lens aperture
nSamples = 25;           % On the first aperture. x,y, before cropping
nSamplesHQ = 801;        % Number of samples for the HQ render
fLength = 50;            % Todo: We should derive this using the lensmaker's equation
lens.apertureMiddleD = 5;

% lens.draw;

%% Pick a point, create its PSF 
% These psfs will be for different field heights, depths, and wavelengths
curPt = 1;
%  curInd = 1

%---initial low quality render
film = pbrtFilmObject('position', [fX fY fZ], 'size', [fW fH], ...
    'wave', wave, 'resolution', [numPixelsW numPixelsH length(wave)]);
psfCamera = psfCameraObject('lens', lens, 'film', film, 'pointsource', pointSources(curPt, :));
ds1 = psfCamera.get('spacing');

% What happens to each of the wavelengths?
oi = psfCamera.estimatePSF();
% To calculate and show, use this:
%   oi = psfCamera.showFilm;
%   oiGet(oi,'spatial resolution','mm')

% Figure out center pos by calculating the centroid of illuminance image
img = oiGet(oi,'illuminance');

% Force to unit area and flip up/down for a point spread
img = img./sum(img(:));
img = flipud(img);
% vcNewGraphWin; mesh(img);

% Calculate the weighted centroid/center-of-mass
xSample = linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1));
ySample = linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2));
[filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);

distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
centroidX = sum(sum(img .* filmDistanceX));
centroidY = sum(sum(img .* filmDistanceY));

sz = oiGet(oi,'size'); mid = round(sz(1)/2); 
u1 = plotOI(oi,'illuminance hline',[mid,mid]);

% Render image using new center position and width and higher resolution
smallFilm = pbrtFilmObject('position', [centroidX centroidY fZ], ...
    'size', [newWidth newWidth], ...
    'wave', wave, ...
    'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);

% Use more samples in the lens aperture to produce a high quality psf.
% NOTE:  Changing the number of samples also changes the oi size.
% This isn't good.  We need to change the sampling density without
% changing the size.
lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
psfCamera = psfCameraObject('lens', lens, 'film', smallFilm, 'pointsource', pointSources(curPt, :));
oi = psfCamera.estimatePSF(nLines, jitterFlag);
ds2 = psfCamera.get('spacing');

% Compare with coarse resolution
%   oi = psfCamera.showFilm;
%   oiGet(oi,'spatial resolution','mm');
sz = oiGet(oi,'size'); mid = round(sz(1)/2);
u2 = plotOI(oi,'illuminance hline',[mid,mid]);

vcNewGraphWin;
s1 = sum(u1.data(:))*ds1;
s2 = sum(u2.data(:))*ds2;
plot(u1.pos,u1.data/s1,'k-',u2.pos,u2.data/s2,'r-')

vcAddObject(oi); oiWindow;

%% Show the lens ray trace
psfCamera.estimatePSF(200);


%% Record on film
psfCamera.recordOnFilm();

% Show the point spread as an image
oi = psfCamera.oiCreate;
img = oiGet(oi,'rgb image');
vcNewGraphWin; image(img); axis image

% Bring up the pointspread in an optics window
psfCamera.showFilm();

% Plot the illuminance image
plotOI(oi,'illuminance mesh linear');

%% Plenoptic

ppsfCamera = ppsfCameraObject('lens', lens, 'film', film, 'pointSource', pointSources(curPt,:));

nLines =  100;  % Draw the ray trace if nLines > 0
ppsf = ppsfCamera.estimatePPSF(nLines);
ppsfCamera.recordOnFilm();
oi = ppsfCamera.showFilm();

%% END