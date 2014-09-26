function oi = createOI(obj, lens, film)
% Create OI from a LF and film 
%
%  oi = LF.createOI(obj, lens, film)
%
% Example:
%
% See Also: ppsfCameraC, rayC
%
% AL Vistasoft Copyright 2014

%% Parameters
% Convert from light field format to ray format (origin/direction)
r         = obj.get('ray');

%rayOrigin = r.get('origin');
%rayDir    = r.get('direction');

% Set up the wavelength information
wave      = obj.get('wave');
%waveIndex = obj.get('waveIndex');

% We should be worried about removing the linear transforms at this stage,
% too.  Say
%
%   obj2 = obj.removeNans
%
%remove nans - aperture NOT specified.  just get rid of nan's

%keep = r.get('liveIndices');
%rayOrigin = rayOrigin(:,keep);   %todo: make this cleaner.  perhaps the ray object can return the non nan's automatically
%rayDir    = rayDir(:,keep);
%waveIndex = waveIndex(keep);

rayOrigin = r.get('origin', 'survivedRaysOnly', true);
rayDir = r.get('direction', 'survivedRaysOnly', true);
waveIndex = obj.get('waveIndex', 'survivedRaysOnly', true);

% Each ray has a wavelength from wave, and which one is stored in waveIndex
calculatedRays = rayC('origin', rayOrigin', ...
    'direction', rayDir',...
    'wave', wave,...
    'waveIndex', waveIndex);
% calculatedRays.plotPhaseSpace();
calculatedRays.recordOnFilm(film);

ppsfCamera = ppsfCameraC('film',film,'lens',lens);
ppsfCamera.ppsfRays = calculatedRays;

oi = ppsfCamera.oiCreate;
% vcAddObject(oi); oiWindow

end

