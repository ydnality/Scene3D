%% s_s3dTransverseCA
%
% Show how changing the aperture position with respect to the lens causes a
% shift in the magnification (transverse chromatic aberration) with respect
% to wavelength
%
% This is built on the notes DHB left us about 18 months ago
%
% AL/BW Vistasoft Team, Copyright 2014


%% Make a lens with an aperture in the middle.
lensFileName = fullfile(s3dRootPath,'data', 'lens', 'dgauss.50mm.dat');
apertureD = 8;
lens = lensC('fileName', lensFileName,'aperture middle d',apertureD);
lens.draw;

% Set the n values for all the lenses to range from 1.3 to 1.7
% Should the aperture (n = 6) have all ones?
nSurfaces = lens.get('n surfaces');
for ii=1:(nSurfaces-1)
    lens.surfaceArray(ii).n = linspace(1.3,1.7,lens.get('nwave'));
end

film = pbrtFilmC;

% Make a point source
ps = [20 0 200];

ppsfCamera = ppsfCameraC('lens',lens,'film',film,'point source',ps);

%% 
ppsfCamera.bbmCreate;

% Still doesn't work.  Ask MP what's going on
ppsfCamera.bbmGetValue('focal')

ppsfCamera.bbmGetValue('effective focal length')



%% End