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
%
% BW Vistasoft Team, Copyright July 2014


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

% Define point sources
pX = [0 10];
pY = [0];           
pZ =[-90:-20:-130];     % Depth range

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
fX = 0; fY = 0; fZ = 100;    % mm.  Film position should go with lens.

% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
% fX = 0; fY = 0; fZ = 135.5;    % mm.  Film position should go with lens.

% Maybe we need lensRead(fname,{param/val pairs}), such as
%  lensRead(lensFile,'wave',wave);
%
load(lensFile,'lens')
lens.set('wave', wave);

% Preview and high quality sampling on first lens aperture
nSamples = 25;           % On the first aperture. x,y, before cropping
nSamplesHQ = 801;        % Number of samples for the HQ render
fLength = 50;            % Todo: We should derive this using the lensmaker's equation
lens.apertureMiddleD = 10;

% lens.draw;

%% Pick a point, create its PSF 

ff = 1; dd = 1;
%---initial low quality render
for ff = 1:nFH
    for dd = 1:nDepth
        film = pbrtFilmObject('position', [fX fY fZ], ...
            'size', [fW fH], ...
            'wave', wave, ...
            'resolution', [numPixelsW numPixelsH length(wave)]);

        lens.apertureSample = ([nSamples nSamples]);
        psfCamera = psfCameraObject('lens', lens, ...
            'film', film, ...
            'pointsource', pointSources{ff,dd});
        ds1 = psfCamera.get('spacing');
        
        % This traces through to the film plane
        % toFilm = true; nLines = 200; psfCamera.draw(toFilm,nLines);
        
        % What happens to each of the wavelengths?
        % This traces to the back of the lens
        oi = psfCamera.estimatePSF();
        % vcAddObject(oi); oiWindow;
        
        % Find the point spread centroid
        centroid = psfCamera.get('image centroid');
        
        % To calculate and show, use this:
        %   oi = psfCamera.showFilm;
        %   oiGet(oi,'spatial resolution','mm')
        sz = oiGet(oi,'size'); mid = round(sz(1)/2);
        [u1,h] = plotOI(oi,'illuminance hline',[mid,mid]); close(h)
        
        % Render image using new center position and width and higher resolution
        smallFilm = pbrtFilmObject('position', [centroid.X centroid.Y fZ], ...
            'size', [newWidth newWidth], ...
            'wave', wave, ...
            'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
        
        % Use more samples in the lens aperture to produce a high quality psf.
        % NOTE:  Changing the number of samples also changes the oi size.
        % This isn't good.  We need to change the sampling density without
        % changing the size.
        lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
        psfCamera = psfCameraObject('lens', lens, ...
            'film', smallFilm, ...
            'pointsource', pointSources{ff,dd});
        oi = psfCamera.estimatePSF(nLines, jitterFlag);
        ds2 = psfCamera.get('spacing');
        % psfCamera.draw(true,500);
        
        % Compare with coarse resolution
        %   oi = psfCamera.showFilm;
        %   oiGet(oi,'spatial resolution','mm');
        sz = oiGet(oi,'size'); mid = round(sz(1)/2);
        [u2, h] = plotOI(oi,'illuminance hline',[mid,mid]); close(h)
        
        vcNewGraphWin;
        s1 = sum(u1.data(:))*ds1;
        s2 = sum(u2.data(:))*ds2;
        plot(u1.pos,u1.data/s1,'k-',u2.pos,u2.data/s2,'r-')
        legend({'Low res','High res'})
        
        vcAddObject(oi); oiWindow;
    end
end


%% Show the lens ray trace
psfCamera.draw(true,200);


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

ppsfCamera = ppsfCameraObject('lens', lens, 'film', film, 'pointSource', pointSources{1,1});

nLines =  100;  % Draw the ray trace if nLines > 0
ppsf = ppsfCamera.estimatePPSF(nLines);
ppsfCamera.recordOnFilm();
oi = ppsfCamera.showFilm();

%% END