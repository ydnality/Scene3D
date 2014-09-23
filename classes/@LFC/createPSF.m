function [oi] = createPSF(obj, ppsfCamera)
%CREATEPSF 

% [oi] = createPSF(obj, ppsfCamera, withinAperture, waveIndexIn)
% given a bEstInterp vector (containing the lightfield at the sensor), and
% a ppsfCamera, use the light field rays and reproject them to construct a
% PSF for visualization.  Also plot the phase space at the sensor.
% TODO: this still needs to be cleaned up.  Maybe there's a better way to
% accomplish what this function accomplishes

    %find the z position of rays
    zPos = ppsfCamera.film.position(3);

    LF = obj.get('LF');
    rayOrigin = zeros(3, size(LF, 2));
    rayDir = rayOrigin;

    rayOrigin(1,:) = LF(1,:);
    rayOrigin(2,:) = LF(2,:);
    rayOrigin(3,:) = 0;

    rayDir(1,:) = LF(3,:);
    rayDir(2,:) = LF(4,:);
    rayDir(3,:) = 1 - rayDir(1,:).^2 + rayDir(2,:).^2;

    %TODO: error checking: if wave is the same for ppsfCameraC and LFC
    wave = obj.get('wave');  
    waveIndex = obj.get('waveIndex');
    
    %TODO: error checking.  Don't rely on luck for these vectors to be the
    %same length.
    waveIndex = waveIndex(~isnan(waveIndex));  %remove nans - aperture NOT specified.  just get rid of nan's
    
%     if(ieNotDefined('withinAperture'))
%        
%     else
%         waveIndex = waveIndex(withinAperture);  %remove nans if aperture is specified  %TODO: see if this is necessary..
%     end


    calculatedRays = rayC('origin', rayOrigin', 'direction', rayDir', 'wave', wave, 'waveIndex', waveIndex);
    calculatedRays.plotPhaseSpace();

    newFilm = pbrtFilmC('position', [0 0 100 ], ...
        'size', [10 10], ...
        'wave', 400:50:700);

    nLines = 100;
    calculatedRays.recordOnFilm(newFilm, nLines); 
    ppsfCamera.film = newFilm; 
    oi = ppsfCamera.oiCreate;
end

