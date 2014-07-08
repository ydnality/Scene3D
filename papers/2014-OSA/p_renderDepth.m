%% p_renderDepth
%
% Show focus as a function of point depth
%
% BW Copyright Vistasoft Team, 2014


%%
s_initISET

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
pX = 0;
pY = 0;                % Assume radial symmetry, so only calculate X
pZ =[-70:-5:-130];    % Depth range

% What is the normalizingZ thing inside of psCreate (BW)?
pointSources = psCreate(pX,pY,pZ);

nDepth = length(pZ);
nFH    = length(pX) * length(pY);

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

% lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
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

% ff = 1; dd = 1;
clear xLine
clear img
ff = 1;
for dd = 1:nDepth
    %---initial low quality render
    film = pbrtFilmObject('position', [fX fY fZ], ...
        'size', [fW fH], ...
        'wave', wave, ...
        'resolution', [numPixelsW numPixelsH length(wave)]);
    
    psfCamera = psfCameraObject('lens', lens, 'film', film, 'pointsource', pointSources{ff,dd});
    ds1 = psfCamera.get('spacing');
    
    % What happens to each of the wavelengths?
    oi = psfCamera.estimatePSF();
    sz = oiGet(oi,'size'); mid = round(sz(1)/2);
    [u1,g] = plotOI(oi,'illuminance hline',[mid,mid]); close(g)
    
    % Find the point spread centeroid
    centroid = psfCamera.get('image centroid');
    
    % Render image using new center position and width and higher resolution
    smallFilm = pbrtFilmObject('position', [fX fY fZ], ...
        'size', [newWidth newWidth], ...
        'wave', wave, ...
        'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
    
    % Use more samples in the lens aperture to produce a high quality psf.
    % NOTE:  Changing the number of samples also changes the oi size.
    % This isn't good.  We need to change the sampling density without
    % changing the size.
    lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
    psfCamera = psfCameraObject('lens', lens, 'film', smallFilm, 'pointsource', pointSources{ff,dd});
    oi = psfCamera.estimatePSF(nLines, jitterFlag);
    ds2 = psfCamera.get('spacing');
    
    % Compare with coarse resolution
    %   oi = psfCamera.showFilm;
    %   oiGet(oi,'spatial resolution','mm');
    sz = oiGet(oi,'size'); mid = round(sz(1)/2);
    [u2,g] = plotOI(oi,'illuminance hline',[mid,mid]); close(g);
    
    % Compare the coarse and fine plots through the center
    vcNewGraphWin;
    s1 = sum(u1.data(:))*ds1;
    s2 = sum(u2.data(:))*ds2;
    plot(u1.pos,u1.data/s1,'k-',u2.pos,u2.data/s2,'r-')
    title(sprintf('Point depth %.1f',pointSources{ff,dd}(3)))
    
    %%
    xLine(:,dd) = u2.data/s2;
    img(:,:,dd) = oiGet(oi,'illuminance');
end


%% Show the spread as a function of depth
vcNewGraphWin;
mesh(abs(pZ),u2.pos,xLine)
ylabel('position')
xlabel('depth')

%%
vcNewGraphWin([],'tall');
colormap(gray)
lst = 1:2:nDepth;
for dd=1:length(lst)
    subplot(length(lst),1,dd), imagesc(img(:,:,lst(dd))); axis image
end


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
%
% Some bug in here.
%
% ppsfCamera = ppsfCameraObject('lens', lens, 'film', film, 'pointSource', pointSources(ff,dd));
%
% nLines =  100;  % Draw the ray trace if nLines > 0
% ppsf = ppsfCamera.estimatePPSF(nLines);
% ppsfCamera.recordOnFilm();
% oi = ppsfCamera.showFilm();

%% END