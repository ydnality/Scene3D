%% ray-tracing for realistic lens - PPSF
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
% AL Vistalab, 2014
%%
s_initISET


%% first PPSF at the center
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 0 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmObject([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%% ray trace and save ppsf
s_s3dRayTracePPSFInterpHelper;
% save ppsf
firstRays = ppsfObject();
firstRays.makeDeepCopy(modifyRays);


%% 2nd PPSF at the center to the side
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 10 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmObject([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%%
s_s3dRayTracePPSFInterpHelper;
% save ppsf
secondRays = ppsfObject();
secondRays.makeDeepCopy(modifyRays);




%% 3rd PPSF in between
%% point sources
pointSourceDepth = 20000;
pointSourceFieldHeight = 0;
pointSources = [ 5 0 -pointSourceDepth];  %large distance test
% pointSources = [ 0 0 -60];  %short distance test

%% film properties -
film = pbrtFilmObject([0 0 36.4],[1 1], 400:10:700, [(400:10:700)' (1:31)'], []);   %large distance

%% lens properties
diffractionEnabled = false;   
%turning on diffraction does NOT make sense yet since we have not modeled 
%the transformation of uncertainty from the middle aperture to the end of the lens

%initialize to default
%lensRealisticObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center, diffractionEnabled, wave)
lens = lensRealisticObject([],[],[],[], 8, 50, [], diffractionEnabled);

%read lens from file
lens.readLensFile(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'))

%intialize ray=tracing
randJitter = false;
lens.calculateApertureSample([301 301], randJitter);

%lens illustration
lens.drawLens();
%%
s_s3dRayTracePPSFInterpHelper;
% save ppsf
middleRays = ppsfObject();
middleRays.makeDeepCopy(modifyRays);


%% Interpolate between 1st and 2nd PPSF to try to produce 3rd one
