function blurScene(obj, overSampleRate)
% Renders the scene in a VoLTCamera according to the included VoLT model,
% scene, and film.  For each point in the scene, a light field is
% calculated from that point to the front most aperture element of the
% lens.  Next, the VoLT will apply a linear transform on this input light
% field to transform it to output light field, which is at the back-most
% element of the lens.  Lastly, the light field is projected onto the film.
%
% overSampleRate: in order to reduce aliasing artifacts, we must sample the
% scene at a higher resolution than the film.  The overSampleRate is a 1x2
% vector that determines this resolution.  For example, an overSampleRate
% of [2 3] will sample the scene twice as much in the x direction, and
% three times as much in the y direction. (Default: [2 2])
%
% Example: VoLTCameraObject = VoLTCameraC(); 
%          VoLTCameraObject.applyOTF([2 2]);
%

if (ieNotDefined('overSampleRate'))
    overSampleRate = [2 2];
end

%tic

%profile on;
filmResolution = obj.film.resolution;
filmPosition = obj.film.position(3);

% Resize scene and all other related vars to the desired size(we want to oversample the film)
unBlurredPhotons = sceneGet(obj.scene, 'photons');
unBlurredPhotons = imresize(unBlurredPhotons, [obj.film.resolution(1) * overSampleRate(1), obj.film.resolution(2) * overSampleRate(2)]);
resizedDepth = imresize(sceneGet(obj.scene, 'depth map'), [obj.film.resolution(1) * overSampleRate(1), obj.film.resolution(2) * overSampleRate(2)]);
blurredPhotons = zeros(filmResolution(1), filmResolution(2), filmResolution(3));

% Declare some useful vars
oiSize = size(unBlurredPhotons); 
center = oiSize./2;
sceneHFOV = sceneGet(obj.scene, 'hfov');
adjustedMiddleApertureRadius = 1;  %this is the size of the middle aperture

wave = obj.film.wave;
% set p the parallel pool

% originally 3580 sec
% with parfor: 376 sec
if (matlabpool('size') > 0)
    matlabpool close;
end

%TODO:need to detect # of cores first
matlabpool open 8;

parfor i = 1:oiSize(2)
    i
    film = pbrtFilmC('position', [0 0 filmPosition ], ...
    'size', [40 40], 'resolution', filmResolution,  ...
    'wave', wave);
    for j = 1:oiSize(1)
        film.clear();
        
        %figure out rotation of pixel
        %this will be the rotation from the positive y = 0 line in a
        %counter-clockwise fashion
        x = -(j - center(2));  %consider vectorizing for speed
        y = i - center(1);
        thetaRad = atan2(x,y);
        thetaDeg = thetaRad/pi * 180;
        
        %figure out field height
        fieldHeight = sqrt((x)^2 + (y)^2); %make into function
        %this is only the field height with respect to the pixels on the
        %sensor
        
        %convert this to the PSF location in 3 space... somehow... using
        %the depth map and some geometry
        wantedPSLocation = [0 0 0];
        
        %wantedPSLocation(2) = fieldHeight/2;  %works as a placeholder
        %wantedPSLocation(2) = fieldHeight/oiSize(2) * resizedDepth(i,j)/filmDistance;
        
        %this gives the current angle with respect to optical axis, when
        %using the radially symmetric field height (in radians)
        currentAngle = fieldHeight/(oiSize(2)/2) * (sceneHFOV/2) * (pi/180);   % figure out FOV stuff... do we use scene FOV or oi FOV?
        %currentDepth = 110; %assumed to be 103 for now for simplicity
        currentDepth = resizedDepth(i,j);
        %wantedPSLocation(3) = -currentDepth;   %old - not completely true
        %wantedPSLocation(2) = tan(currentAngle) * currentDepth;
        wantedPSLocation(2) = sin(currentAngle) * currentDepth;
        wantedPSLocation(3) = -sqrt(currentDepth^2 - wantedPSLocation(2)^2); 
        %wantedPsLocation(3) = -resizedDepth(i,j);
        %wantedPSLocation = [0 15 -103]; %some testing with PS locations
        
        % --- Interpolate PSF for current point in scene ---
        %first get the linear transform
        LTObject = obj.VoLTObject.interpolateAllWaves(wantedPSLocation);
        
        % calculate the lightfield from this particular point source
        [inputLF]  = s3dLightFieldEntrance(wantedPSLocation, obj.lens);   %traces rays to entrance only. 
        
        % Make an LT (linear transform) object and apply the LT on the inputLF
        outputLFObject = LTObject.applyOnLF(inputLF, adjustedMiddleApertureRadius);
        
        % Apply linear rotation transform on LF
        %thetaRad = theta/180 * pi;
        rotationMatrix = [cos(thetaRad)    sin(thetaRad)      0           0;
            -sin(thetaRad)   cos(thetaRad)      0           0;
            0             0               cos(thetaRad)  sin(thetaRad)
            0             0               -sin(thetaRad) cos(thetaRad)];
        
        rotationMatrixFull = repmat(rotationMatrix, [1 1 length(wave)]);
        RotationObject = LTC('AInterp', rotationMatrixFull, 'wave', wave);
        rotatedLFObject = RotationObject.applyOnLF(outputLFObject, adjustedMiddleApertureRadius);
        
        % Visualize PSF and phase space
        oiI = rotatedLFObject.createOI(obj.lens,film);
        psfPhotons = oiGet(oiI, 'photons');
        
        %convert the oi of the PSF and multiply by spectral radiance of the
        %scene, and add together, in order to blur the scene
        
        if(isnan(psfPhotons(:)))
            warning('nan photons');   
            %problem: when oiI is all 0's... for whatever reason, the min
            %and max is not set correctly and results in nans.  Investigate
            %this in the future;
            psfPhotons = zeros(size(psfPhotons));   %this is a hack for now
        end
        
        %weigh each channel of PSf according to the weight of the original
        %Unblurred Image
        illuminanceWeight = repmat(unBlurredPhotons(i,j, :), [size(psfPhotons, 1) size(psfPhotons, 2)]);
        
        %add blurred PSF to the existing sum
        blurredPhotons = blurredPhotons + psfPhotons .* illuminanceWeight; 
    end
end

% Close the matlab pool
if (matlabpool('size') > 0)
    matlabpool close;
end

% Create an oi and assign blurred photons to oi
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi, 'wave', wave);
oi = oiSet(oi, 'cphotons', blurredPhotons);
oi = oiSet(oi,'hfov', sceneHFOV);
vcAddObject(oi); oiWindow;

%profile viewer;
%profile clear; 
%toc

end