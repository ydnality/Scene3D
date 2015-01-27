%% p_renderFH
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
pX = 0; %linspace(0,3.5,1);        
pY = linspace(-18,18,3);
pZ = -100;        % Depth range

% What is the normalizingZ thing inside of psCreate (BW)?
pointSources = psCreate(pX,pY,pZ);

nDepth = length(pZ);
nFH    = length(pX) * length(pY);

jitterFlag = true;   % Enable jitter for lens front element aperture samples
nLines = false;      % Number of lines to draw for debug illustrations.

%%  Declare film properties for PSF recording.

wave = 400:100:700;            % Wavelength
fX = 0; fY = 0; fZ = 100;    % mm.  Film position

% Film resolution (final render)
numPixelsWHQ = 300;
numPixelsHHQ = 300;

% The film width for high quality is set manually - this
% should be automated in the future
newWidth = 75;    % mm

%% Describe the lens

% lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
load(lensFile,'lens')
[p,lens.name] = fileparts(lensFile);
lens.set('wave',wave);

% Linearly space the index of refraction, just as an example.
W  = [400, 700];
nW = [1.720, 1.5500];

% Preview and high quality sampling on first lens aperture
nSamples = 25;           % On the first aperture. x,y, before cropping
nSamplesHQ = 801;        % Number of samples for the HQ render
fLength = 50;            % Todo: We should derive this using the lensmaker's equation
lens.apertureMiddleD = 10;

% lens.draw;

%% Pick a point, create its PSF
% These psfs will be for different field heights, depths, and wavelengths

% ff = 1; dd = 1;
wbar = waitbar(0);
% Could add a loop on depths for dd = 1:nDepth
dd = 1;
for ff = 1:nFH
    waitbar(ff/nFH,wbar,sprintf('Field height %d',ff));
    
    % Render image using new center position and width and higher resolution
    smallFilm = filmC('position', [fX fY fZ], ...
        'size', [newWidth newWidth], ...
        'wave', wave, ...
        'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
    
    % Use more samples in the lens aperture to produce a high quality psf.
    % NOTE:  Changing the number of samples also changes the oi size.
    % This isn't good.  We need to change the sampling density without
    % changing the size.
    lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
    psfCamera = psfCameraC('lens', lens, ...
        'film', smallFilm, ...
        'pointsource', pointSources{ff,dd});
    
    psfCamera.estimatePSF(nLines, jitterFlag);
    oi = psfCamera.oiCreate;
    % vcAddObject(oi); oiWindow;

    if ff == 1
        sz = size(oiGet(oi,'rgb image'));
        rgb = zeros(sz(1),sz(2),3,nFH);
        %         vcNewGraphWin; image(rgb(:,:,:,3))
        ps1 = psfCamera;
    elseif ff == round(nFH/2)
        ps2 = psfCamera;
    elseif ff == nFH
        ps3 = psfCamera;
    end
    
    rgb(:,:,:,ff) = oiGet(oi,'rgb image');
end
close(wbar);

%% Make a video of the PSFs as a function of wavelength

vObj = VideoWriter('fhPSF.avi','Motion JPEG AVI');
vObj.FrameRate = 4;
open(vObj);

vcNewGraphWin;
v = uint8(ieScale(rgb,0,1)*255);
s = oiGet(oi,'sample spacing','mm');
x = ((1:sz(2)) - (sz(2)/2))*s(2);
y = ((1:sz(1)) - (sz(1)/2))*s(1);

for ii=1:nFH
    % v = uint8(ieScale(img(:,:,ii),0,1)*255);
    image(x,y,v(:,:,:,ii)); axis image
    text(x(60),y(5),lens.name,'Color','w','FontSize',20);
    grid on; set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
    xlabel('mm'); ylabel('mm')
    drawnow;
    
    f = getframe(gcf);
    writeVideo(vObj,f);
end

close(vObj);

%%
toFilm = true; nLines = 75;
ps1.draw(toFilm,nLines);
set(gca,'ylim',[-30 30],'xlim',[-120 100])
ps2.draw(toFilm,nLines);
set(gca,'ylim',[-30 30],'xlim',[-120 100])
ps3.draw(toFilm,nLines);
set(gca,'ylim',[-30 30],'xlim',[-120 100])

%%
oi = ps1.oiCreate;
vcAddObject(oi); oiWindow;


oi = ps2.oiCreate;
vcAddObject(oi); oiWindow;

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