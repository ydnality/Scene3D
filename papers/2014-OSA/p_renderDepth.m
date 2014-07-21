%% p_renderDepth
%
% Show focus as a function of point depth
%
% BW Copyright Vistasoft Team, 2014


%%
s_initISET

%% Render a PSF collection
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
% What is the normalizingZ thing inside of psCreate (BW)?
pX = 0;
pY = 0;                % Assume radial symmetry, so only calculate X
pZ =[-70:-10:-100];    % Depth range
pointSources = psCreate(pX,pY,pZ);

nDepth = length(pZ);
nFH    = length(pX) * length(pY);

jitterFlag = true;   % Enable jitter for lens front element aperture samples
nLines = false;      % Number of lines to draw for debug illustrations.

%%  Declare film properties for PSF recording.

wave = 500;
% wave = 400:100:700;            % Wavelength
wList = 1:length(wave);
fX = 0; fY = 0; fZ = 135.5;      % mm.  Film position

% Film width and height for coarse calculation
% fW = 80;  % mm
% fH = 80;  % mm

% Film resolution (preview)
% numPixelsW = 151;
% numPixelsH = 151;

% Film resolution (final render)
numPixelsWHQ = 100;
numPixelsHHQ = 100;

% The film width for high quality is set manually - this
% should be automated in the future
newWidth = 10;    % mm

%% Describe the lens

lensFile = fullfile(s3dRootPath, 'data', 'lens', 'dgauss.50mm.mat');
% lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.mat');
load(lensFile,'lens')
lens.apertureMiddleD = 5;

[p,lens.name] = fileparts(lensFile);
lens.set('wave',wave);
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
wbar = waitbar(0,'Depth renderings');
for dd = 1:nDepth
    
    % Render image using new center position and width and higher resolution
    smallFilm = pbrtFilmC('position', [0 0 fZ], ...
        'size', [newWidth newWidth], ...
        'wave', wave, ...
        'resolution', [numPixelsWHQ numPixelsHHQ length(wave)]);
    
    % Use more samples in the lens aperture to produce a high quality psf.
    % NOTE:  Changing the number of samples also changes the oi size.
    % This isn't good.  We need to change the sampling density without
    % changing the size.
    waitbar(dd/nDepth,wbar,'Estimating PSFs');
    lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
    
    % Create the camera, calculate the PSF and create an OI
    psfCamera = psfCameraC('lens', lens, 'film', smallFilm, 'pointsource', pointSources{ff,dd});
    psfCamera.estimatePSF(nLines, jitterFlag);
    oi = psfCamera.oiCreate;
    % vcAddObject(oi); oiWindow;
    
    % Save the point spread data for plotting
    if dd == 1
        sz = oiGet(oi,'size'); 
        mid = round(sz(1)/2);
        img = zeros(sz(1),sz(2),nDepth);
        xLine = zeros(sz(1),nDepth);
        ps1 = psfCamera; 
    elseif dd == round(nDepth/2)
        ps2 = psfCamera; 
    elseif dd == nDepth
        ps3 = psfCamera; 
    end
    
    [uData,g] = plotOI(oi,'illuminance hline',[mid,mid]); close(g);
    if dd == 1,
        pos = uData.pos/1000; % Return is um, converts to mm
        xLine = zeros(length(pos),nDepth);
    end
    img(:,:,dd) = oiGet(oi,'illuminance');
    xLine(:,dd) = uData.data;
    
end
close(wbar);


%%  Make a video of the different depth point spreads

vObj = VideoWriter('depthPSF.avi','Motion JPEG AVI');
vObj.FrameRate = 5;
open(vObj);

vcNewGraphWin; colormap(gray)
v = uint8(ieScale(img,0,1)*255);
[r, c] = size(v(:,:,1));
x = pos;
y = pos;

for ii=1:nDepth
    % v = uint8(ieScale(img(:,:,ii),0,1)*255);
    image(x,y,v(:,:,ii)); brighten(0.5); axis image;
    text(x(5),y(5),sprintf('%.0f ',-pZ(ii)),'Color','w','FontSize',20);
    text(x(60),y(5),lens.name,'Color','w','FontSize',20);
    title(sprintf('%s',lens.name));
    grid on; set(gca,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
    xlabel('mm'); ylabel('mm')
    drawnow;
    
    f = getframe(gcf);
    writeVideo(vObj,f);
end

close(vObj);

%% Plot the spread as a function of depth
vcNewGraphWin;
mesh(abs(pZ),pos,xLine)
xlabel('depth (mm)')
ylabel('position (mm)')
zlabel('Illuminance (lux)')

%% Show ray traces for three depths

ps1.draw(true,200); 
set(gca,'xlim',[-100 150]);
ps2.draw(true,200); 
set(gca,'xlim',[-100 150]);
ps3.draw(true,200); 
set(gca,'xlim',[-100 150]);


%% Record on film
psfCamera.recordOnFilm();

% Show the point spread as an image
oi = psfCamera.oiCreate;

% Plot the illuminance image
plotOI(oi,'illuminance mesh linear');

% Bring up the oiWindow;
vcAddObject(oi); oiWindow;


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