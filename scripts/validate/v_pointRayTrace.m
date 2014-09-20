%% Validate creation of a lens and single point ray trace
%
%  Read in a lens file
%  Create a point source
%  Create a film
%  Ray trace the point through the lens to the film
%
% AL/BW VISTASOFT 2014

%%
s_initISET

%% Make a volume of point.  A little boring in this case

point = psCreate(0,1,-120);

%% Read a lens file and create a lens

lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

%% Create a film (sensor)

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
wave = lens.get('wave');

film = pbrtFilmC('position', [0 0 100 ], ...
    'size', [10 10], ...
    'wave', wave);


%% Ray trace the point to the film

camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
camera.estimatePSF;
oi = camera.oiCreate;

vcAddObject(oi); oiWindow;

%% END
