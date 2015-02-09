%% Test Spherical sensor capability
%
%  Read in a lens file
%  Create a point source
%  Create a film
%  Ray trace the point through the lens to the film
%
%  Repeat with spherical sensor and observe the difference
% AL/BW VISTASOFT 2014

%%
ieInit

%% Make a volume of point.  A little boring in this case

%point = psCreate(0,15,-120);

point = psCreate(0,20,-120);
%point = psCreate(sqrt(112.5), sqrt(112.5),-120);
%point = psCreate(15,0,-120);
%point = psCreate(-sqrt(112.5), sqrt(112.5),-120);

%% Read a lens file and create a lens

lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');

nSamples = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
    'fileName', lensFileName, ...
    'apertureMiddleD', apertureMiddleD);

%% make lens have chromatic aberrations

wave = lens.get('wave');
numSurfaces = lens.get('nsurfaces');
offset = .05;
for i = 1:numSurfaces
   curN = lens.surfaceArray(i).get('n');
   
   if (curN(1)~=0 && curN(1)~=1)
       nVector = linspace(curN(1) - offset, curN(1) + offset, length(wave));
       lens.surfaceArray(i).set('n', nVector);
   end
end


%% Create a film (planar sensor)

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
wave = lens.get('wave');

film = filmC('position', [0 0 95 ], ...
    'size', [50 50], ...
    'wave', wave);


%% Ray trace the point to the film

camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
camera.estimatePSF(20, true);
oi = camera.oiCreate;
oi = oiSet(oi,'name','Planar sensor');
vcAddObject(oi); oiWindow;


%% Create a film (spherical sensor)

% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
wave = lens.get('wave');

film = filmSphericalC('position', [0 0 95 ], ...
    'size', [50 50], ...
    'wave', wave, ...
    'radius', -25);


%% Ray trace the point to the film

camera = psfCameraC('lens',lens,'film',film,'point source',point{1});

% Sequence of events for estimating the PSF, 
camera.estimatePSF(20, true);
oi = camera.oiCreate;
oi = oiSet(oi,'name','Spherical sensor');
vcAddObject(oi); oiWindow;

%% END
