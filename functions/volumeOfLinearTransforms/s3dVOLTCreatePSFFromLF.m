function ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp, withinAperture, waveIndexIn)
% ppsfCamera = s3dVOLTCreatePSFFromLF(ppsfCamera, bEstInterp)
% given a bEstInterp vector (containing the lightfield at the sensor), and
% a ppsfCamera, use the light field rays and reproject them to construct a
% PSF for visualization.  Also plot the phase space at the sensor.
% TODO: consider putting this in a class later (either ppsfCamera, ppsf, or
% something else)
% TODO: this still needs to be cleaned up.  Maybe there's a better way to
% accomplish what this function accomplishes
%
% waveIndexIn: sometimes a new waveIndex must be specified, depending how
% bEstInterp was constructed.  In multiple wavelength mode, we will
% concatenate rays from different wavelenghts.  This is kept track of using
% a custom waveIndex list.  This is when this would be used
%
    %find the z position of rays
    zPos = ppsfCamera.film.position(3);

    rayOrigin = zeros(3, size(bEstInterp, 2));
    rayDir = rayOrigin;

    rayOrigin(1,:) = bEstInterp(1,:);
    rayOrigin(2,:) = bEstInterp(2,:);
    rayOrigin(3,:) = 0;

    rayDir(1,:) = bEstInterp(3,:);
    rayDir(2,:) = bEstInterp(4,:);
    rayDir(3,:) = 1 - rayDir(1,:).^2 + rayDir(2,:).^2;

    wave = ppsfCamera.ppsfRays.get('wave');  %check to see if this is right... used to be ppsf.get
    waveIndex = ppsfCamera.ppsfRays.get('waveIndex');
    
    
    if (ieNotDefined('waveIndexIn'))
        if(ieNotDefined('withinAperture'))
            waveIndex = waveIndex(~isnan(waveIndex));  %remove nans - aperture NOT specified.  just get rid of nan's
        else
            waveIndex = waveIndex(withinAperture);  %remove nans if aperture is specified  %TODO: see if this is necessary..
        end
    else
        waveIndex = waveIndexIn;
    end
    
    
    calculatedRays = rayC('origin', rayOrigin', 'direction', rayDir', 'wave', wave, 'waveIndex', waveIndex);
    calculatedRays.plotPhaseSpace();

    newFilm = pbrtFilmC('position', [0 0 100 ], ...
        'size', [10 10], ...
        'wave', 400:50:700);

    nLines = 100;
    calculatedRays.recordOnFilm(newFilm, nLines); 

    ppsfCamera.film = newFilm; 
end