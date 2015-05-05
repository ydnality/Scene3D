%% p_renderOiMatlabToolFull
%
% 2014 OSA Conference: Full Forward calculation
%
% This script takes a scene image with depth map, and applies the series of
% PSF's to it, depending on the position, wavelength, and depth.
%
%
% AL Copyright Vistasoft Team 2014

%%
ieInit

%% The scene file, rendered using a pinhole model
sceneFileName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'pinholeSceneFile.mat');
load(sceneFileName,'scene');
vcAddAndSelectObject(scene); sceneWindow;

%% Basic point source properties (distances)

% Make a function that provides points along a constant ray whose centroid
% should be at different distances along the same field height.
wave = 400:100:700;            % Wavelength

% We will loop through the point positions
pX = 0:-1500:-4500;
pY = 0;
pZ =[-70 -80 -90 -100 -110];  % millimeters
normalizingZ = -16000;        % mm assumed reference Z point.
% All other points will use the same field angle as this reference Z point

% I don't understand the normalization here.  This probably has to do with
% ray angles.  Ask AL (BW).
% Actually, the pZ(ii) below doesn't make any sense.
pointSources = cell(nFH,nDepth);
for ii=1:nFH
    for dd = 1:nDepth
        pointSources{ii,dd} = [pX(ii) *  (pZ(ii)/normalizingZ), pY, pZ(dd)];
    end
end

nWave = length(wave);
nFH = length(pX) * length(pY);
nDepth = length(pZ);
nPoints = nFH*nDepth;

jitterFlag   = true;   %enable jitter for lens front element aperture samples
nLines       = false;  %number of lines to draw for debug illustrations.

% We should plot the points in 3 Space
% vcNewGraphWin
% for ff=1:nFH
%     for dd=1:nDepth
%         pt = pointSources{ff,dd};
%         plot3(pt(1),pt(2),pt(3),'o');
%         hold on
%     end
% end
% xlabel('X'); ylabel('Y'); zlabel('Z')

%%  Declare film properties for PSF recording (NOT for the forward calculation)

fX = 0; fY = 0; fZ = 135.5;    % mm.  Film position

% Film width and height
fW = 120;  % mm
fH = 120;  % mm

% Film resolution (preview, large film size)
lowRes = 151;   % Number of sample points across the film surface

% Film resolution (final render, small film size)
highRes = 75;   % Number of samples across the whole film

%for now - the width of the high quality sensor is set manually - this should be
%somewhat automated in the future
newWidth = 10;    %mm

%% Create the lens from a lens file
% see TwoElLens for how this lens was created
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens');
load(lensFile,'lens')
lens.apertureMiddleD = 5;

% Preview and high quality sampling on first lens aperture
nSamples = 25;           % On the first aperture. x,y, before cropping
nSamplesHQ = 801;        % Number of samples for the HQ render

% lens.draw;

%% Make the low resolution camera
% We will send this in to the PSFArray routine.

film = filmC('position', [fX fY fZ], ...
    'size', [fW fH], ...
    'wave', wave, ...
    'resolution', [lowRes lowRes length(wave)]);

lens.apertureSample = ([nSamples nSamples]);
psfCamera = psfCameraC('lens', lens, ...
    'film', film, ...
    'pointsource', pointSources{1,1});

%% Should replace cell below

PSF = psfCamera.PSFArray(pointSources);

vcNewGraphWin;
colormap(gray)
% Could loop through these
% ww = 3; 
% for ff=1:4
%     for dd=1:5
%         s = PSF(:,:,ww,dd,ff);
%         imagesc(s);
%         pause(0.2);
%     end
% end

%% Original calculation
% 
% wbar = waitbar(0,sprintf('Creating %i point spreads ...',numel(pointSources)));
% oiList = cell(nFH,nDepth);
% 
% % ii = 3; dd = 2;
% for ii = 1:nFH
%     waitbar(ii/nFH,wbar);
%     
%     for dd = 1:nDepth
%         
%         %--- Initial low quality render
%         film = filmC('position', [fX fY fZ], ...
%             'size', [fW fH], ...
%             'wave', wave, ...
%             'resolution', [lowRes lowRes length(wave)]);
%         
%         lens.apertureSample = ([nSamples nSamples]);
%         psfCamera = psfCameraC('lens', lens, ...
%             'film', film, ...
%             'pointsource', pointSources{ii,dd});
%         
%         % From here could could be psfCamera.get('image centroid')
%         
%         % What happens to each of the wavelengths?
%         % psfCamera.estimatePSF();
%         % psfCamera.oiCreate;
%         
%         % Replaced by oiGet centroid, I think
%         
% %         % Figure out center pos by calculating the centroid of illuminance image
% %         img = oiGet(oi,'illuminance');
% %         
% %         % Force to unit area and flip up/down for a point spread
% %         img = img./sum(img(:));
% %         img = flipud(img);
% %         % vcNewGraphWin; mesh(img);
% %         
% %         % Calculate the weighted centroid/center-of-mass
% %         % The units here are millimeters, I think (BW)
% %         xSample = linspace(-film.size(1)/2, film.size(1)/2, film.resolution(1));
% %         ySample = linspace(-film.size(2)/2, film.size(2)/2, film.resolution(2));
% %         [filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);
% %         
% %         distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
% %         centroidX = sum(sum(img .* filmDistanceX));
% %         centroidY = sum(sum(img .* filmDistanceY));
%         
%         % The centroid calculation produces the same answer
%         % disp([ii dd])
% %         fprintf('field height %i\nDepth %i',ii,dd);
% %         disp([centroidX, centroidY])
%         
%         c = psfCamera.centroid();
%         centroidX = c.X; centroidY = c.Y;
%         fprintf('field height %i\nDepth %i',ii,dd);
%         disp([centroidX, centroidY])
% 
%         % sz = oiGet(oi,'size'); mid = round(sz(1)/2);
%         
%         % Render image using new center position and width and higher resolution
%         film = filmC('position', [centroidX centroidY fZ], ...
%             'size', [newWidth newWidth], ...
%             'wave', wave, ...
%             'resolution', [highRes highRes length(wave)]);
%         
%         % Use more samples in the lens aperture to produce a high quality psf.
%         % NOTE:  Changing the number of samples also changes the oi size.
%         % This isn't good.  We need to change the sampling density without
%         % changing the size.
%         lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
%         psfCamera = psfCameraC('lens', lens, ...
%             'film', film, ...
%             'pointsource', pointSources{ii,dd});
%         psfCamera.estimatePSF(nLines, jitterFlag);
%         oiList{ii,dd} = psfCamera.oiCreate();
%         
%         % vcAddObject(oiList{1,1}); oiWindow;
%     end
% end
% 
% delete(wbar)
% 
% % Form PSF matrix
% PSF = zeros(highRes, highRes, nWave, nDepth, nFH);
% for ww = 1:nWave
%     for dd = 1:nDepth
%         for ii = 1:nFH
%             PSF(:,:,ww,dd,ii) = oiGet(oiList{ii,dd}, 'photons',wave(ww));
%         end
%     end
% end

%%
% Key data to know for interpolation later.
% This should become the ray trace structure in ISET.

% We store the field heights in terms of angle
PSFStructure.fHAngle = atan(pX/normalizingZ) * 180/pi;

% The depth is in mm
PSFStructure.depth   = -pZ(:);

% Wavelength, PSF array, and so forth
PSFStructure.wave    = wave;
PSFStructure.PSF     = PSF;   % The spatial samples are in???
PSFStructure.film    = film;

% Make some plots of this PSF structure data.
% What we really want is just this:  To make the PSFStructure here
% compatibility with the ray trace point spread functions in ISET, and
% then to be able to run oiCompute in raytrace mode with the depth
% image.
%
% Doing this involves eliminating the code Andy copied from ISET, but
% upgrading the ISET code to 3D.

%% This section appears to determine the size of the output image


% % Dimensions of scene data
% numRows = sceneGet(scene, 'rows');
% numCols = sceneGet(scene,'cols');
% 
% % Size of PSF film pixel
% smallPixelSize = 75/size(PSFStructure.film.image,1);     %mm
% 
% % Size of scene sample pixel
% largePixelSize = 70/sqrt(2)/numRows;  %mm
% 
% %should match the PSFStructure film Z to be consistent
% filmDistance = PSFStructure.film.position(3);
% 
% %70/sqrt(2) represents row size/diagonal  TODO: find a way to automate this
% %This size should correspond to the size sensor you wish to use for the
% %forward calculation.
% 
% % Amount that PSF film needs to be scaled to be equivalent to scene sample
% %size
% scaleFactor = smallPixelSize/largePixelSize;

%%  Copy the scene data into an oi structure?
renderWave = 400:100:700;

% Create an oi
oi = oiCreate;
oi = initDefaultSpectrum(oi);

% Use scene data as photons
scene   = sceneSet(scene,'wave', renderWave);
dM      = sceneGet(scene, 'depth map');
photons = sceneGet(scene, 'photons');

% Apply this to the oi
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi, 'photons', photons);

% This should be set so that the resolution of the PSF samples match the
% resolution of the oi
oi = oiSet(oi,'optics focal length',abs(fZ)*1.1*1e-3);
oiGet(oi,'sample spacing','mm')
% Adjust to make the spacing 100 um, or ...

% Pad the oi 
padAmount = oiGet(oi,'size')*0.1;
oi = oiPad(oi, padAmount);
vcAddObject(oi); oiWindow;

% Obtain unblurred image, to be used later
unBlurredPhotons = oiGet(oi, 'photons');

% This is where the output will go
photonSum = zeros(size(oiGet(oi,'photons')));

%% Apply proper PSF for every pixel - this will become a function

% For each pixel,
% - Find field height, angle, wavelength, and depth
% - Find the corresponding PSF
% - Multiply the PSF by the pixel radiance
% - Add the result to the full image

numRows = sceneGet(scene, 'rows');
numCols = sceneGet(scene,'cols');
filmDistance = PSFStructure.film.position(3);

% Size of scene sample pixel
% Not sure where the 70 comes from.
% Not sure about this whole thing here.
largePixelSize = 70/sqrt(2)/numRows;  %mm
smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm
scaleFactor    = smallPixelSize/largePixelSize;

[X,Y] = meshgrid(1:numCols,1:numRows);
X = X - (numCols/2);
Y = Y - (numRows/2);
hyp   = sqrt( Y.^2 + X.^2)*largePixelSize;
% vcNewGraphWin; mesh(hyp);

% I think the X-axis is supposed to be zero deg.
% This isn't quite there yet.  But it is close.
ang = atan2d(Y, X);
% vcNewGraphWin; mesh(ang);

%% Precompute the PSFStructure to match the spatial sampling of the OI

% The PSF is sampled at 100 microns, we think.  
% The OI has a different sample.  Resample the PSFs to match the OI spatial
% sampling.

% We could also replace all of the PSFs with a functional form, say a
% Gaussian with a couple of parameters.  Then produce the Gaussians with
% the appropriate sampling as we step through the field heights.  The
% parameters that are interpolated are the Gaussian parameters, not the PSF
% values.

%% Set up the dM, angle, and field height maps for each pixel
fhMap
dMap
aMap

%% Loop through wavelength
wBar = waitbar(0,'Slow forward rendering');
for waveInd = 1:length(renderWave)
    % Should be able to stop looping on wavelength
    
    curWave = renderWave(waveInd);
    fprintf('Wavelength %.0f\n',curWave);
    
    % Loop through pixel rows and cols
    for ii = 1:numRows
        waitbar(ii/numRows,wBar);
        for jj = 1:numCols
            
            % Calculate field height (in degrees) and angle(relative to
            % center)
            % hypotenuse = sqrt((ii - numRows/2)^2 + (jj - numCols/2)^2) * largePixelSize;
            hypotenuse = hyp(ii,jj);
            curFieldHeightAngle = atan(hypotenuse/filmDistance) * 180/pi;
            
            % Pixel depth
            depth = dM(ii, jj);
            
            % Calculate the PSF for the current field height, angle, depth,
            % wavelength, given the PSf structure
            currentPSF = s3dPSFLookUp(curFieldHeightAngle, depth, curWave, PSFStructure);
            
            %rotatePSF to correct orientation
            scaledCPSF = imrotate(currentPSF, ang(ii,jj), 'bilinear', 'crop' );
            
            % Scale PSF to the right size
            scaledCPSF = imresize(scaledCPSF, scaleFactor);
            scaledCPSF = scaledCPSF./sum(scaledCPSF(:)); %normalize
            
            scaledPSFNumRows = size(scaledCPSF, 1);
            scaledPSFNumCols = size(scaledCPSF, 2);
            
            % Define center and starting window
            centerPos = ceil(size(scaledCPSF, 1)/2);
            
            % Put it in the photonSum variable
            startingRow = ii + padAmount - centerPos;
            endingRow   = startingRow + scaledPSFNumRows - 1;
            
            startingCol = jj + padAmount - centerPos;
            endingCol   = startingCol + scaledPSFNumCols - 1;
            
            % Add to the final image
            photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) = ...
                photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) + ...
                unBlurredPhotons(ii + padAmount,jj + padAmount, waveInd)*scaledCPSF;
        end
    end
end
delete(wBar);

% Assign sum as oi values
oi = oiSet(oi, 'wave', renderWave);
oi = oiSet(oi, 'photons',photonSum);
vcAddObject(oi); oiWindow;


% %% Apply proper PSF for every pixel - this will become a function
% 
% % For each pixel,
% % - Find field height, angle, wavelength, and depth
% % - Find the corresponding PSF
% % - Multiply the PSF by the pixel radiance
% % - Add the result to the full image
% 
% numRows = sceneGet(scene, 'rows');
% numCols = sceneGet(scene,'cols');
% filmDistance = PSFStructure.film.position(3);
% 
% % Size of scene sample pixel
% % Not sure where the 70 comes from.
% % Not sure about this whole thing here.
% largePixelSize = 70/sqrt(2)/numRows;  %mm
% smallPixelSize = PSFStructure.film.size(1)/size(PSFStructure.film.image,1);     %mm
% scaleFactor    = smallPixelSize/largePixelSize;
% 
% [X,Y] = meshgrid(1:numCols,1:numRows);
% X = X - (numCols/2);
% Y = Y - (numRows/2);
% hyp   = sqrt( Y.^2 + X.^2)*largePixelSize;
% % vcNewGraphWin; mesh(hyp);
% 
% % I think the X-axis is supposed to be zero deg.
% % This isn't quite there yet.  But it is close.
% ang = atan2d(Y, X);
% % vcNewGraphWin; mesh(ang);
% 
% % Loop through wavelength
% wBar = waitbar(0,'Slow forward rendering');
% for waveInd = 1:length(renderWave)
%     
%     curWave = renderWave(waveInd);
%     fprintf('Wavelength %.0f\n',curWave);
%     
%     % Loop through pixel rows and cols
%     for ii = 1:numRows
%         waitbar(ii/numRows,wBar);
%         for jj = 1:numCols
%             
%             % Calculate field height (in degrees) and angle(relative to
%             % center)
%             % hypotenuse = sqrt((ii - numRows/2)^2 + (jj - numCols/2)^2) * largePixelSize;
%             hypotenuse = hyp(ii,jj);
%             curFieldHeightAngle = atan(hypotenuse/filmDistance) * 180/pi;
%             
%             % Pixel angle
%             angle = ang(ii,jj);
%             
%             % Pixel depth
%             depth = dM(ii, jj);
%             
%             % Calculate the PSF for the current field height, angle, depth,
%             % wavelength, given the PSf structure
%             currentPSF = s3dPSFLookUp(curFieldHeightAngle, depth, curWave, PSFStructure);
%             
%             %rotatePSF to correct orientation
%             scaledCPSF = imrotate(currentPSF, angle, 'bilinear', 'crop' );
%             
%             % Scale PSF to the right size
%             scaledCPSF = imresize(scaledCPSF, scaleFactor);
%             scaledCPSF = scaledCPSF./sum(scaledCPSF(:)); %normalize
%             
%             scaledPSFNumRows = size(scaledCPSF, 1);
%             scaledPSFNumCols = size(scaledCPSF, 2);
%             
%             % Define center and starting window
%             centerPos = ceil(size(scaledCPSF, 1)/2);
%             
%             startingRow = ii + padAmount - centerPos;
%             endingRow   = startingRow + scaledPSFNumRows - 1;
%             
%             startingCol = jj + padAmount - centerPos;
%             endingCol   = startingCol + scaledPSFNumCols - 1;
%             
%             % Add to the final image
%             photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) = ...
%                 photonSum(startingRow:endingRow,startingCol:endingCol, waveInd) + ...
%                 unBlurredPhotons(ii + padAmount,jj + padAmount, waveInd)*scaledCPSF;
%         end
%     end
% end
% delete(wBar);
% 
% % Assign sum as oi values
% oi = oiSet(oi, 'wave', renderWave);
% oi = oiSet(oi, 'photons',photonSum);
% vcAddObject(oi); oiWindow;

%%
scene = sceneCreate;
scene = sceneSet(scene,'wave',oiGet(oi,'wave'));

p = oiGet(oi,'photons');

scene = sceneSet(scene,'photons',p);
vcAddObject(scene);
sceneWindow;

% % debugging
% vcNewGraphWin;
% imshow(photonSum(:,:,1)./ max(photonSum(:)))
% vcNewGraphWin;
% imshow(unBlurredPhotons(:,:,1)./ max(unBlurredPhotons(:)))

%% Debugging s3dLookUp: Key data: PSF and PSFFieldHeightSamples

% wavelength = 700;
% fieldHeightAngle = 0;
% depth = 4000;
% PSFStructure
%
% %TODO: keep track of differences in size maybe?  For now, assume app PSFs
% %are plotted on the same size sensor
%
% currentPSF = s3dPSFLookUp(fieldHeightAngle, depth, wavelength,PSFStructure);
% vcNewGraphWin; imagesc(currentPSF); title(['PSF Depth: ' int2str(depth)]);


%% For comparison this is the backwards calculation

% Reduce to 400:100:700
% Note: there is no way for a scene to be restructured for wavelengths so
% we need to make this into an Oi and do it this way.  Perhaps this should
% be a new feature of ISE to allow for scenes to change their wave samples
% and interpolate.

scene = sceneCreate;
% scene = sceneSet(scene, 'photons', oiGet(oi,' photons'));
scene = sceneSet(scene, 'photons', oiGet(opticalimage,' photons'));
scene = sceneSet(scene, 'wave', 400:100:700);
% If PBRT is installed, we can run this to compare with the backward
% calculation.

% sceneName = fullfile(s3dRootPath, 'data', 'pbrtScenes', 'metronome', 'mainDefocused2El.pbrt');
% oi = s3dRenderOI(sceneName, .050, 'indObj');
%
% %reduce to 400:100:700
% scene = sceneCreate;
% scene = sceneSet(scene, 'photons', oiGet(oi,' photons'));
% scene = sceneSet(scene, 'wave', 400:100:700);
%
% oi = oiSet(oi, 'wave', 400:100:700);
% oi = oiSet(oi, 'photons', sceneGet(scene, 'photons'));
% vcAddObject(oi); oiWindow;

%%
