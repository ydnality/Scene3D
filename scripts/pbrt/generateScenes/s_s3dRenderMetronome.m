%% Render the Metronome scene with a standard camera or wide field
%
%  With some parameter settings below, we also can render the scene with a
%  very wide field of view and with a curved sensor. 
%  
%  Read the comments to find the key place where those parameters are set.
%
%  There is also a lightfield camera version of this rendering.
%
% AL VISTASOFT, 2014

%%
ieInit

%% Rendering parameters
param.wideField    = false;
param.sensorCurved = false;
param.addLight     = true;
param.pinhole      = false;

% Governs the number of samples per pixel. A minimum is 32.
% For the scarlet one we ran it with 12*2048.
nSamples = 2*2048;

%%
tic

%initialization
clear curPbrt;
curPbrt = pbrtObject();

% Specify properties files (materials and geometry)
matFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-mat.pbrt');
geoFile = fullfile(dataPath, 'twoFlashDepth', 'indObject', 'pbrt', 'default-geom.pbrt');

% light properties
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]*1e5);
lightFrom = [-56.914787 -105.385544 35.0148];  % Position of source
lightTo =   [-56.487434 -104.481461 34.8  ];   % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  
    
%% Add spot light and infinite light sources

if param.addLight
    % Make the spot light
    spotLight = pbrtLightSpotObject();
    spotLight.setName('spot');
    spotLight.setSpectrum(pbrtSpectrumObject('rgb I', [1000 1000 1000]*1e5));
    spotLight.setAngle(180);
    spotLight.setDeltaAngle(180);
    spotLight.setFrom([-142.3855 -286.2024  13.0082]);
    spotLight.setTo([ -141.9582 -285.2984   13.0082]);
    curPbrt.addLightSource(spotLight);
    
    %infinite light (for diffuse lighting)
    infiniteLight = pbrtLightInfiniteObject();
    curPbrt.addLightSource(infiniteLight);
    
end

%% Standard camera properties
position = ...
    [-56.914787 -105.385544 35.0148;    % From
     -56.487434 -104.481461 34.8 ;      % To
      0 0 1];
% position = [from; to; 0 0 1];
curPbrt.camera.setPosition(position);

%% specify the wide field camera lens and aperture
    
if param.pinhole
    lens         = pbrtLensPinholeObject();
    filmDistance = 140;
    filmDiag     = 50.9117;
    curPbrt.camera.setLens(pbrtLensPinholeObject(filmDistance, filmDiag));
    % curPbrt.camera.setResolution(300, 300);
else
    curPbrt.camera.setLens(pbrtLensRealisticObject());
    curPbrt.camera.lens.filmDistance = 90; % 90 70;  % 133.33;
    curPbrt.camera.lens.filmDiag     = 70;  % 70
    % curPbrt.camera.lens.specFile = '2ElLens.dat';
    curPbrt.camera.lens.specFile      = 'dgauss.50mm.dat';
    curPbrt.camera.lens.apertureDiameter = 3;  % 16 for narrow DOF, 1 for deep (mm)
    curPbrt.camera.lens.curveRadius = 0;       % Experimental
end
curPbrt.camera.setResolution(450, 300);


%% Wide field camera properties
% 
% Use this to over-write the standard lens, above.
% Simplify the logic, though the code might be OK.
%
if param.wideField
    % Move the camera back 80 millimeters from the standard position
    v1 = [-56.914787 -105.385544 13.014802];
    v2 = [-56.487434 -104.481461 13.014835;];
    camDir = v2 - v1;
    camDir = camDir./norm(camDir, 2);
    camPos = [v1 + camDir * 80 ;
        v2 + camDir* 80 ;
        -0.000013 -0.000031 1.000000;];
    
    curPbrt.camera.setPosition(camPos);
    curPbrt.camera.setLens(pbrtLensRealisticObject());
    
    curPbrt.camera.lens.filmDistance = 15;  %90; %70; %133.33;
    curPbrt.camera.lens.filmDiag = 43.75; %; 70;
    curPbrt.camera.lens.specFile = '2ElLens13.5mm.dat';
    curPbrt.camera.lens.apertureDiameter = .5; % in mm
    
    curPbrt.camera.setResolution(450, 300);  
end

%% AL experiments with curved sensor

% The sensor is built inside of PBRT
if param.sensorCurved
    % Experimental work with the curved sensor (even though it says lens)
    curPbrt.camera.lens.curveRadius = -13.5;  % This is really the SENSOR curve radius.
    % If we set this to 0, then the sensor is flat.
end

%% add old parts, put in new ones
% Candidate for a function
curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Render the oi

dockerFlag = true;
oiName = 'metronome';
oi = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

toc

vcAddAndSelectObject(oi);
oiWindow;

%% Save  the  variables that enable running the previous cell
%
% This could be saved on scarlet.

save('metronome','oi','curPbrt','param');

%% END
