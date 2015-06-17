% Runs PBRT and imports the result into an ISET oi for the 'cone' scene. 
%
% Creates the light field OI used to illustrate transverse chromatic
% aberration in AL's dissertation.  The image has a black line at the
% bottom that can be zoomed to show the different wavelength dispersion.
%
% NOTES: 
% The number of samples needs to be large to make a reasonable image.  The
% default is set to low resolution (see 'quick' below to increase
% resolution).
%
% This routine could be made more efficient using MP's code based on the
% effective aperture.  In general, if we could integrate the bbm code with
% what we are doing, that would be good.
%
% AL Copyright Vistasoft Team 2015

%%
ieInit;

%% Call a PBRT wrapper to make the scene cones.pbrt

%initialization
clear curPbrt;
curPbrt = pbrtObject();

includeFile = fullfile(dataPath, 'pbrtScenes', 'cones', 'cones.pbrt');
if ~exist(includeFile,'file'), error('Missing key include file'); end

curPbrt.addInclude(includeFile);
curPbrt.camera.addTransform(pbrtTransformObject('Scale', [1000 1000 1000]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate',  [-3 1 0 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Rotate', [52 0 1 0]));
curPbrt.camera.addTransform(pbrtTransformObject('Translate', [-2.3 -.05 .5]));

%% Make a lens and add to pbrt object
quick    = false; 
if quick, nSamples = 2048;             % This is too few, but it gives the idea
else      nSamples = 2048*10;
end

filmDist = 40;           % This brings the front part of the image into reasonable focus
filmDiag = 10;
specFile = 'dgauss.50mm.dat';
apertureDiameter    = 3;  % mm
diffraction         = false;
chromaticAberration = true;

% set sampler
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

% Create and set lenss
lens = pbrtLensRealisticObject(filmDist, filmDiag, specFile, apertureDiameter, diffraction, chromaticAberration, []);
curPbrt.camera.setLens(lens);

%% make an area light and add to pbrt object
light = pbrtAreaLightObject('area', pbrtSpectrumObject('color L', [1000 1000 1000]));
areaLightShape = pbrtShapeObject('disk', 'radius', 8);
light.addShape(areaLightShape);
light.addTransform(pbrtTransformObject('Translate', [ 0 9.9 0]));
light.addTransform(pbrtTransformObject('Rotate', [90 1 0 0]));

% light.removeProperty();
light.addProperty(pbrtPropertyObject('integer nsamples', 4));

% Remove the default and add the one that was just made
curPbrt.removeLight();
curPbrt.addLightSource(light);

% An optional alternative light source
% curPbrt.removeLight();
% curPbrt.addLightSource(pbrtLightInfiniteObject('infiniteLight', 16, [], [], []));
% 

% Remove the defaults
curPbrt.removeMaterial();
curPbrt.removeGeometry();

%% Run it

oiName     = 'coneArray';
dockerFlag = true;
frontOi    = s3dRenderOIAndDepthMap(curPbrt, oiName, dockerFlag);

% Visualize it
vcAddObject(frontOi); oiWindow;

%% END
