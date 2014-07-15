%% p_renderWave
%
% Show focus as a function of point depth
%
% BW Copyright Vistasoft Team, 2014


%%
s_initISET

%% Render a PSF collection by wavelength
%
% Rendering the collection will become a function at some point
%
% The high-resolution PSFs will be added into a PSF Collection matrix,
% which contains the whole series of PSFs depending on field height,
% wavelength, depth.
%
% The PSF collection is used to produce the forward calculation
% rendered image.

% We will loop through the point positions.  Units are millimeters
pX = 0;
pY = 0;          % Assume radial symmetry, so only calculate X
pZ =-100;        % Depth range

% What is the normalizingZ thing inside of psCreate (BW)?
pointSources = psCreate(pX,pY,pZ);

nDepth = length(pZ);
nFH    = length(pX) * length(pY);

jitterFlag = true;   % Enable jitter for lens front element aperture samples
nLines = false;      % Number of lines to draw for debug illustrations.

%%  Declare film properties for PSF recording.

wave = 400:20:700;            % Wavelength
fX = 0; fY = 0; fZ = 100;    % mm.  Film position

% Film width and height for coarse calculation
fW = 80;  % mm
fH = 80;  % mm

% Film resolution (preview)
numPixelsW = 151;
numPixelsH = 151;

% Film resolution (final render)
numPixelsWHQ = 100;
numPixelsHHQ = 100;

% The film width for high quality is set manually - this
% should be automated in the future
newWidth = 10;    % mm

%% Describe the lens

lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');

% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
load(lensFile,'lens')
[p,lens.name] = fileparts(lensFile);

lens.apertureMiddleD = 5;
W = 400:100:700;
nW = [1.720,     1.6567,     1.6033,     1.5500];

% Preview and high quality sampling on first lens aperture
nSamples = 25;           % On the first aperture. x,y, before cropping
nSamplesHQ = 801;        % Number of samples for the HQ render
fLength = 50;            % Todo: We should derive this using the lensmaker's equation
lens.apertureMiddleD = 5;

% lens.draw;

%% Pick a point, create its PSF
% These psfs will be for different field heights, depths, and wavelengths

% ff = 1; dd = 1;
nWave = length(wave);
wbar = waitbar(0);
for ww = 1:nWave
    waitbar(ww/nWave,wbar,sprintf('Wave %d',wave(ww)));
    lens.set('wave',wave(ww));
    n = interp1(W,nW,wave(ww));
    lens.set('n all',n);
    
    % Render image using new center position and width and higher resolution
    smallFilm = pbrtFilmObject('position', [fX fY fZ], ...
        'size', [newWidth newWidth], ...
        'wave', wave(ww), ...
        'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
    
    % Use more samples in the lens aperture to produce a high quality psf.
    % NOTE:  Changing the number of samples also changes the oi size.
    % This isn't good.  We need to change the sampling density without
    % changing the size.
    lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
    psfCamera = psfCameraObject('lens', lens, ...
        'film', smallFilm, ...
        'pointsource', pointSources{1,1});
    
    oi = psfCamera.estimatePSF(nLines, jitterFlag);
    if ww == 1
        sz = size(oiGet(oi,'rgb image'));
        rgb = zeros(sz(1),sz(2),3,nWave);
    end
    rgb(:,:,:,ww) = oiGet(oi,'rgb image');
    % vcAddObject(oi); oiWindow;

end
close(wbar);

%% Make a video of the PSFs as a function of wavelength

vObj = VideoWriter('wavePSF.avi','Motion JPEG AVI');
vObj.FrameRate = 5;
open(vObj);

vcNewGraphWin;
v = uint8(ieScale(rgb,0,1)*255);
s = oiGet(oi,'sample spacing','mm');
x = ((1:sz(2)) - (sz(2)/2))*s(2);
y = ((1:sz(1)) - (sz(1)/2))*s(1);

for ii=1:nWave
    % v = uint8(ieScale(img(:,:,ii),0,1)*255);
    image(x,y,v(:,:,:,ii)); 
    text(x(5),y(5),sprintf('%.0f ',wave(ii)),'Color','w','FontSize',20);
    text(x(60),y(5),lens.name,'Color','w','FontSize',20);
    title(sprintf('%s',lens.name));
    grid on; set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
    xlabel('mm'); ylabel('mm')
    drawnow;
    
    f = getframe(gcf);
    writeVideo(vObj,f);
end

close(vObj);


%% Show the spread as a function of depth
% Not implemented here.  Get it from Depth if you want it.
% vcNewGraphWin;
% mesh(abs(pZ),u2.pos,xLine)
% ylabel('position')
% xlabel('depth')

%% Show the lens ray trace
% psfCamera.draw(200);


%% Plot the illuminance image
% plotOI(oi,'illuminance mesh linear');

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