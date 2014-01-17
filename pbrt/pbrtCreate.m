function pbrt = pbrtCreate
% Create a pbrt object that we will use for writing a pbrt file
%
% We use blender to create graphics files that PBRT can render.
%
% This function creates a pbrt structure that includes a variety of
% parameters (top level) for telling pbrt how to render a specific set of
% objects.
%
% The parameters that can be set are
%
%   LIST TO GO HERE
%
%
% (AL) VISTASOFT TEAM Copyright 2013

%% Initialize the pbrt structure
pbrt.name = 'default';
pbrt.type = 'pbrt';

% This is some blender to pbrt coordinate frames (AL)
pbrt.scale = [-1 1 1];

% Placement of the camera.  Defined by a vector that looks in a certain
% direction, and then tells you which way is up.
pbrt.cameraPosition = [
    4.5 -80 7 % Starting up
    4.5 -79 7 % Ending up 
    0 0 1];   % Which way is up

% Example lens
pbrt.lensFile = 'idealLens-50mm.pbrt';

% Typical sensor
pbrt.film.name = 'image';
pbrt.film.xresolution = 200;
pbrt.film.yresolution = 200;

% 
% Sampler
pbrt.sampler.type = 'lowdiscrepancy';
pbrt.sampler.pixelsamples = 512;

% SurfaceIntegrator
pbrt.surfaceIntegrator.type  = 'directlighting';
pbrt.surfaceIntegrator.maxdepth = 0;

% Renderer
pbrt.renderer = 'sample';

% World Begin

%  Attribute Begin
%   Light source
pbrt.lightSource = cell(1,1);

%default white light at camera position
whiteLight.type = 'spot';
whiteLight.spectrum.type = 'rgb I';
whiteLight.spectrum.value = [1000 1000 1000];
whiteLight.coneangle = 180;
whiteLight.conedeltaangle = 180;
whiteLight.from = [4.5 -90 8.5];
whiteLight.to = [4.5 -89 8.5];

%green spotlight
greenLight.type = 'spot';
greenLight.spectrum.type = 'rgb I';
greenLight.spectrum.value = [0 1000 0];
greenLight.coneangle = 180;
greenLight.conedeltaangle = 180;
greenLight.from = [4.5 -90 8.5];
greenLight.to = [4.5 -89 8.5];

%red spotlight
redLight.type = 'spot';
redLight.spectrum.type = 'rgb I';
redLight.spectrum.value = [1000 0 0];
redLight.coneangle = 180;
redLight.conedeltaangle = 180;
redLight.from = [4.5 -90 8.5];
redLight.to = [4.5 -89 8.5];

pbrt.lightSource{1} = greenLight;
%  Attribute End

%  Material file
pbrt.materialsFile = 'depthTargetDepths-mat.pbrt';
%  Geometry file
pbrt.geometryFile = 'depthTargetDepths-geom.pbrt';
% WorldEnd


end
