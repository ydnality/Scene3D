function PSFArray = PSFArray(psfCamera,pointSources,varargin)
% Loop on wavelength and depth to create a cell array of PSFs
%
%  If we put this in the psfCamera object, then we could call
%
%   PSFArray = psfCamera.PSFArray(pointSources)
%
%  That gives us an array of PSFs that we can apply in the forward (slow)
%  rendering tool, p_renderOiMatlabToolFull.m, for example.
%
%   pointSources = makeMe
%   PSFArray = psfCamera.PSFArray(pointSources);
%
% The psfCamera has a film and lens.  
% We need to adjust the resolution somehow.
% We would use the varargin to set the high or low spatial resolution of
% the PSF data.
%
% AL/BW Vistasoft 2015

% These psfs will be for different field heights, depths, and
% wavelengths
wbar = waitbar(0,sprintf('Creating %i point spreads ...',numel(pointSources)));
oiList = cell(nFH,nDepth);
for ii = 1:nFH
    waitbar(ii/nFH,wbar);
    
    for dd = 1:nDepth
        
        %--- Initial low quality render
        film = filmC('position', [fX fY fZ], ...
            'size', [fW fH], ...
            'wave', wave, ...
            'resolution', [lowRes lowRes length(wave)]);
        
        lens.apertureSample = ([nSamples nSamples]);
        psfCamera = psfCameraC('lens', lens, ...
            'film', film, ...
            'pointsource', pointSources{ii,dd});
        
        % From here could could be psfCamera.get('image centroid')
        
        % What happens to each of the wavelengths?
        psfCamera.estimatePSF();
        oi = psfCamera.oiCreate;
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
        
        % to here
        
        sz = oiGet(oi,'size'); mid = round(sz(1)/2);
        
        % Render image using new center position and width and higher resolution
        film = filmC('position', [centroidX centroidY fZ], ...
            'size', [newWidth newWidth], ...
            'wave', wave, ...
            'resolution', [highRes highRes length(wave)]);
        
        % Use more samples in the lens aperture to produce a high quality psf.
        % NOTE:  Changing the number of samples also changes the oi size.
        % This isn't good.  We need to change the sampling density without
        % changing the size.
        lens.apertureSample = ([nSamplesHQ nSamplesHQ]);
        psfCamera = psfCameraC('lens', lens, ...
            'film', film, ...
            'pointsource', pointSources{ii,dd});
        psfCamera.estimatePSF(nLines, jitterFlag);
        oiList{ii,dd} = psfCamera.oiCreate();
        
        % vcAddObject(oiList{1,1}); oiWindow;
    end
end

delete(wbar)

%% Compute PSF collection matrix.
% This will serve as a lookup table for later parts of the script

% Form PSF matrix.
% See if we can move this up into the loop above.
PSF = zeros(highRes, highRes, nWave, nDepth, nFH);
for ww = 1:nWave
    for dd = 1:nDepth
        for ii = 1:nFH
            PSF(:,:,ww,dd,ii) = oiGet(oiList{ii,dd}, 'photons',wave(ww));
        end
    end
end

%% Key data to know for interpolation later.
% This should become the ray trace structure in ISET.

% PSFFieldHeightSamples = atan(pX/normalizingZ) * 180/pi;
% PSFDepthSamples = -pZ(:);
PSFArray.fHAngle = atan(pX/normalizingZ) * 180/pi;
PSFArray.depth   = -pZ(:);
PSFArray.wave    = wave;
PSFArray.PSF     = PSF;
PSFArray.film    = film;

% Make some plots of this PSF structure data.


% What we really want is just this:  To make the PSFArray here
% compatibility with the ray trace point spread functions in ISET, and
% then to be able to run oiCompute in raytrace mode with the depth
% image.
%
% Doing this involves eliminating the code Andy copied from ISET, but
% upgrading the ISET code to 3D.