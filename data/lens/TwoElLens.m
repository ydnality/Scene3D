%% TwoElLens.m
%
% Specification of the 2ElLens File data
%
% This is a simple multi element lens properties
%
% This is one way to make the lens.  BW created a .dat file in the proper
% format to make this lens just by calling 
% This goes from the light through the lens to the film
wave = 400:100:700;
zPos   = [-3 -1.5 0];   % Z intercept positions of lens surfaces
radius   = [67 0 -67];    % Radius of curvature, 0 means aperture
aperture = [10 10 10];       % Circular apertures, these are the radii in mm

% Index of refraction to the right of each surface
%(ray.wavelength - 550) * -.04/(300) + curEl.n;
firstN = (wave - 550) * -.04/(300) + 1.65; %linearly changes the 1.65 material
n = [firstN' zeros(length(wave), 1) ones(length(wave),1)]; %index of refraction (wavelength x element)

nSamples = 25;           % On the first aperture. x,y, before cropping
diffractionEnabled = false;    %disable diffraction for this example
idx = find(radius==0);  % This is the middle of the lens aperture size mm
fLength = 50;           % Todo: We should derive this using the lensmaker's equation
% For multiple lenses, we add up the power using something from the web

% Populate lens surface array ussing given property arrays
lensSurfaceArray = lensSurfaceObject();
for i = 1:length(zPos)
    lensSurfaceArray(i) = lensSurfaceObject('sRadius', radius(i), 'apertureD', aperture(i), 'wave',wave,'zPos', zPos(i), 'n', n(:, i));
end

% Declare lens
lens = lensMEObject('surfaceArray', lensSurfaceArray, ...
    'focalLength', fLength, ...
    'diffractionEnabled', diffractionEnabled,...
    'wave', wave, ...
    'aperturesample', [nSamples nSamples]);
lens.apertureMiddleD = 5;
% lens.draw;

lensFile = fullfile(s3dRootPath,'data','lens','2ElLens.mat');
save(lensFile,'lens')

%% The PBRT lens data are stored as 2ElLens.dat 
%  But this has no way to include the additional index of refraction 

% But, we haven't got a way to represent the variations in n with
% wavelength 
lensFile = fullfile(s3dRootPath, 'data', 'lens', '2ElLens.dat');
wave = 400:10:700;
fLength = 50;
nSamples = 25;           % On the first aperture. x,y, before cropping
diffractionEnabled = false;    % disable diffraction for this example
lens = lensMEObject('fileName', lensFile, ...
    'focalLength', fLength, ...
    'diffractionEnabled', diffractionEnabled, ...
    'wave', wave, ...
    'aperturesample', [nSamples nSamples]);
lens.apertureMiddleD = 5;
lens.draw;

%% END
