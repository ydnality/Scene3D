%% Validate lens ray tracing of a single point
%
%  Read in a lens file
%  Create a point source
%  Create a film
%  Visualize Ray trace the point through the lens to the film
%  Create an optical image of the ray trace
%
% AL/BW VISTASOFT 2014

% We could also do this for a couple of film distances and point distances
%
ieInit

%% Make a point far away.  A little off center and 100 mm from the back surface

point = psCreate(0,2,-1000);

%% Read a lens file and create a lens

lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');

nSamples = 251;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

%% Create a film (sensor)

% wavelength samples
wave = lens.get('wave');

% position - relative to center of final lens surface
%   Positive is where the image is formed, negative is where the objects
%   are
% size - 'mm' 

% In focus is about 100 mm
film = filmC('position', [0 0 40], ...
    'size', [10 10], ...
    'wave', wave);

%% Ray trace the point to the film

camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
nLines = 100;
jitter = true;
camera.estimatePSF(nLines,jitter);

%% Show the point spread in the optical image window

oi = camera.oiCreate;
vcAddObject(oi); oiWindow;


%% END
