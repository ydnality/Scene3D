%% Render a slanted bar target
%
% Goal is to emphasize the chromatic aberration measurement
%
% As things stand, we need too many samples.  We really want to have the
% rays be confined to the aperture.  Or we want to blur a little bit in a
% wavelength independent way.
% You can see the chromatic aberration here, but the picture isn't yet
% quite what we want.
%
% AL Vistasoft, 2015

%%
ieInit
clear curPbrt;

%% render scene with PBRT using pbrtObjects (front flash)

% Specify properties files (materials and geometry)
matFile = fullfile(dataPath, 'pbrtScenes', 'slantedBar', 'slantedBar-mat.pbrt');
geoFile = fullfile(dataPath, 'pbrtScenes', 'slantedBar', 'slantedBar-geom.pbrt');

% This is the camera position, taken from the teacup.pbrt file.  I think
% that must be the 'LookAt' parameter.  Notice that they are the same as
% the light positions, above.
position = [
    0 0 8001.0
    0 0 8000.0
    0.000000 1.00000 0.0000005];

% light properties - DOn't understand all of these, but in particular where
% the direction came from.
spectrum = pbrtSpectrumObject('rgb I', [1000 1000 1000]);

% Put the light at the camera
lightFrom = position(1,:);
lightTo =   position(2,:);
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???
lightSource = pbrtLightSpotObject('light', spectrum, coneAngle, coneDeltaAngle, lightFrom, lightTo);  

%% Set up the pbrt material structures and properties
curPbrt = pbrtObject();
curPbrt.camera.setPosition(position);

curPbrt.removeMaterial();
curPbrt.addMaterial(matFile);
curPbrt.removeGeometry();
curPbrt.addGeometry(geoFile);
curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
nSamples = 8192;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Set up film and lens - from slantedBar.pbrt
filmDist = 36.4;   
filmDiag = 0.125;   

numPinholesW = [];   % Set if you want a light field
numPinholesH = [];

% Film sample resolution
rows = 128;
cols = 128;

% Lens propertiess
specFile = 'dgauss.50mm.dat'; 
apertureDiameter    = 2;         
diffraction         = false;
chromaticAberration = true;

%assign pinhole position to PBRT, and figure out correct cropWindow
lens = pbrtLensRealisticObject(filmDist, filmDiag, ...
    specFile, apertureDiameter, ...
    diffraction, chromaticAberration, [], [], [], ...
    numPinholesW, numPinholesH);
curPbrt.camera.setLens(lens);

curPbrt.camera.setResolution(rows, cols);

%% render
focalLength = 0.050;   % In meters
oiName = [];
dockerFlag = true;
oi = s3dRenderOIAndDepthMap(curPbrt,focalLength,oiName,dockerFlag);
vcAddObject(oi); oiWindow;
plotOI(oi,'irradiance hline',[1 round(rows/2)])

%% To save depends on whether light field or note
if isempty(numPinholesW) || isempty(numPinholesH)
    save('slantedBar','oi','focalLength','curPbrt');
else
    save('slantedBarLF','oi','numPinholesW','numPinholesH','focalLength','curPbrt');
end

%% END
% 
% %% PBRT will run the PBRT script
% chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
% 
% % list of all chromatic aberration renderings - uncomment and run the one
% % you wish to run
% 
% % slanted bar rendering
% unix([fullfile(pbrtHome, '/src/bin/pbrt') 'chromaticAberration.pbrt']);
% 
% %% ISET will read the PBRT output
% % scene = sceneSet(scene,'fov', 8);
% oi = pbrt2oi('output_d.dat');
% % oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
% vcAddAndSelectObject(oi);
% oiWindow;
% 
% m = oiGet(oi, 'mean illuminance')
% unix('cd ..');
