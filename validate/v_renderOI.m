%% v_renderOI
%
%  Validate the docker rendering code for a low resolution slanted bar scene
%
% AL/BW

%% Quickly load the pbrt object of a slanted bar

pbrt = pbrtC('slanted bar');

oi = s3dRenderOI(pbrt,'slantedbar',true);

% Run some assert this and that on the oi ...
%

%% End