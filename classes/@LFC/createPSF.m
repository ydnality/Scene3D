function [oi] = createPSF(obj, ppsfCamera)
% CREATE PSF from a plenoptic psf camera
%
%  [oi] = createPSF(obj, ppsfCamera)
%
% A ppsfCamera, use the light field rays and reproject them to construct a
% PSF for visualization.  Also plot the phase space at the sensor.
% TODO: this still needs to be cleaned up.  Maybe there's a better way to
% accomplish what this function accomplishes
%
% AL Vistasoft Copyright 2014

%find the z position of film 
% zPos = ppsfCamera.film.position(3);

%% Parameters

% Get the light field matrix
LF = obj.get('LF');

% Convert from light field format to ray format (origin/direction)
% We do this so we can put the rays into a ppsf object to create the PSF.
%
% Could become
%   [rayOrigin, rayDir] = lf2ray(LF);
%
% or
%
%  [rayO,rayD] = LF.get('ray')
%
r = LF.get('ray');
rayOrigin = r.rayOrigin;
rayDir    = r.rayDir;

% rayOrigin = zeros(3, size(LF, 2));
% rayDir = rayOrigin;
% 
% rayOrigin(1,:) = LF(1,:);
% rayOrigin(2,:) = LF(2,:);
% rayOrigin(3,:) = 0;        % Probably should be a variable name
% 
% rayDir(1,:) = LF(3,:);
% rayDir(2,:) = LF(4,:);
% rayDir(3,:) = 1 - rayDir(1,:).^2 + rayDir(2,:).^2;

%TODO: error checking: if wave is the same for ppsfCameraC and LFC

% Set up the wavelength information
wave = obj.get('wave');
waveIndex = obj.get('waveIndex');

% We should be worried about removing the linear transforms at this stage,
% too.  Say
%
%   obj2 = obj.removeNans
%
%remove nans - aperture NOT specified.  just get rid of nan's
waveIndex = waveIndex(~isnan(waveIndex));  

% Each ray has a wavelength from wave, and which one is stored in waveIndex
calculatedRays = rayC('origin', rayOrigin', 'direction', rayDir',...
    'wave', wave, 'waveIndex', waveIndex);
calculatedRays.plotPhaseSpace();

film = filmC('position', [0 0 100 ], ...
    'size', [10 10], ...
    'wave', 400:50:700);

nLines = 100; % Number of lines used in creating the ray images
calculatedRays.recordOnFilm(film, nLines);

ppsfCamera.film = film;

oi = ppsfCamera.oiCreate;

end

