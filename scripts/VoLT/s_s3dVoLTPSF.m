%% Create a PSF using the VoLT process
%
%  Create a point
%  Create a lens
%  Create a film
%  Calculate point spread on the film
%
% This illustrates that the aperture can be changed quickly in the
% simulation, without recomputing the VOLT class
%
% AL Vistalab, 2014
%%
s_initISET

%% Specify point,lens and film.

%desired pSLocation for interpolation
pSource = [0 0.0001 103];  %note that it is now in [dummyvalue phi (degrees)  depth] format!  (depth MUST BE POSITIVE!!!)

%  Lens 
lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');
nSamples     = 151;
apertureMiddleD = 8;   % mm
lens = lensC('apertureSample', [nSamples nSamples], ...
             'fileName', lensFileName, ...
             'apertureMiddleD', apertureMiddleD);

wave = lens.get('wave');

% Film
% position - relative to center of final lens surface
% size - 'mm'
% wavelength samples
film = filmC('position', [0 0 100 ], ...
    'size', [2 2], ...
    'wave', wave);


%% Compute VOLT model for this point

%VoLTObject = VoLTC('lens', lens, ...
%    'film', film, ...
%    'field positions', pSource(1:2), 'depths', pSource(3), ...
%    'wave', wave);

VoLTObject = VoLTC('lens', lens, ...
    'film', film, ...
    'field positions', [0 1], 'depths', pSource(3), ...
    'wave', wave);

% Calculate the critical ABCD matrices


debug = false;
VoLTObject.calculateMatrices(debug);


%% Calculate the PSF, using the VoLT matrices and the aperture

% Calculate the input light field.  See notes in file  ....
%s3dlightField still uses the old format!

% sampledLocations = VoLTObject.getPSLocations();
% pSource = sampledLocations( 1, :);

pSource = [ 0 0 -103];   %***TODO:  We need to figure out which direction the polar goes in with 0 phi.  There is a discrepency between s3dLightField, and VoltObject.LT.  They use different coordinates.
%[~,~,inputLF]  = s3dLightField(pSource, lens);
[inputLF]  = s3dLightFieldEntrance(pSource, lens);


%polar pSource
polarPSource = [0 0 103];
LT = VoLTObject.LT(polarPSource);

% Apply the LT on the inputLF given a radius
adjustedMiddleApertureRadius = 8;
outputLFObject = LT.applyOnLF(inputLF, adjustedMiddleApertureRadius);

%% Visualize PSF and phase space
%why is this broken??
oiI = outputLFObject.createOI(lens,film);
oiI = oiSet(oiI,'name','Light Field');
vcAddObject(oiI); oiWindow;

uI = plotOI(oiI,'illuminance hline',[1 135]);
title(sprintf(oiGet(oiI,'name')));


%%
