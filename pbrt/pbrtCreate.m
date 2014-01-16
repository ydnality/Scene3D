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

% SurfaceIntegrator
% Sampler
% Renderer

% World Begin
%  Attribute Begin
%   Light source
%  Attribute End
%  Material file
%  Geometry file
% WorldEnd


end
