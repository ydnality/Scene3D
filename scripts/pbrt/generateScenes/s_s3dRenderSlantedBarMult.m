%% Render a slanted bar target
%
% Renders 3 targets with slanted bars.  Each at a different depth.  This
% will be used for MTF evaluation for a LF camera in the future.  For now
% it is good to demonstrate defocus abilities.
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
    0 0 0.0
    0 0 -1
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


matRGB= [1 1 1];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('color Kd', matRGB));
curPbrt.addMaterial(newMaterial);


curPbrt.removeGeometry();
%curPbrt.addGeometry(geoFile);
backDropDepth = -400;
backDropTransform = ...
    [100 0 0 0;
    0 100 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
backDrop = pbrtGeometryObject('backdrop', 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

foregroundDepth = -200;
%add a foreground target
foregroundTransform = ...
    [8 0 0 0;
    0 8 0 0 ;
    0 0 1 0;
    0 0 foregroundDepth  1];
frontSquare = pbrtGeometryObject('slantedBar1', 'Material', [], [], foregroundTransform);
curPbrt.addGeometry(frontSquare);

foregroundDepth2 = -300;
%add a foreground target
foregroundTransform2 = ...
    [8 0 0 0;
    0 8 0 0 ;
    0 0 1 0;
    30 0 foregroundDepth2  1];
frontSquare = pbrtGeometryObject('slantedBar2', 'Material', [], [], foregroundTransform2);
curPbrt.addGeometry(frontSquare);

foregroundDepth0 = -100;
%add a foreground target
foregroundTransform0 = ...
    [5 0 0 0;
    0 5 0 0 ;
    0 0 1 0;
    -10 0 foregroundDepth0  1];
frontSquare = pbrtGeometryObject('slantedBar0', 'Material', [], [], foregroundTransform0);
curPbrt.addGeometry(frontSquare);


curPbrt.removeLight();
curPbrt.addLightSource(lightSource);

% set sampler
curPbrt.sampler.removeProperty();
nSamples = 256;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Set up film and lens - from slantedBar.pbrt
filmDist = 68; %36.4;   
filmDiag = 35;   

numPinholesW = 80;   % Set if you want a light field
numPinholesH = 80;

% Film sample resolution
rows = 720;   %these must be a multiple of numPinholesW and numPinholesH
cols = 720;

% Lens propertiess
%specFile = 'dgauss.50mm.dat'; 
specFile = '2ElLens.dat';
apertureDiameter    = 12;         
diffraction         = false;
chromaticAberration = false;

%assign pinhole position to PBRT, and figure out correct cropWindow
lens = pbrtLensRealisticObject(filmDist, filmDiag, ...
    specFile, apertureDiameter, ...
    diffraction, chromaticAberration, [], [], [], ...
    numPinholesW, numPinholesH);
curPbrt.camera.setLens(lens);

curPbrt.camera.setResolution(rows, cols);

%% render
oiName = 'multSlantedBar';
dockerFlag = true;
oi = s3dRenderOIAndDepthMap(curPbrt,oiName,dockerFlag);

vcAddObject(oi); oiWindow;
% plotOI(oi,'irradiance hline',[1 round(rows/2)])

%% To save depends on whether light field or note
if isempty(numPinholesW) || isempty(numPinholesH)
    save('slantedBarMult','oi','curPbrt');
else
    save('slantedBarMultLF80720','oi','numPinholesW','numPinholesH','curPbrt');
end

%% END

