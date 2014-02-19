classdef filmObject <  handle
    % Create a film object
    %
    %   lens = filmObject(filmDistance,filmDiag);  % Units are mm
    %
    % Presently we only represent spherical lenses and apertures.
    %
    % These are defined by a series of surfaces. We code the offset to each
    % surface radius to the center of the spherical lens.  Positive means
    % to the right (or left???). Aperture parameters (a single number is a
    % diameter in mm). index of refraction (n) for the material to the left
    % of the surface.
    % 
    % pinhole cameras have no aperture and the pinhole lens will inherit
    % this superclas. This will be a superclass that will be inherited by
    % other classes in the future
    %
    % We aim to be consistent with the PBRT lens files, and maybe the Zemax
    % as far possible ?
    %
    % This could become a camera, or we could make a camera object that has
    % a lens and film.
    % 
    % Example:
    %   lensObject
    %   lensObject(30,250)
    %
    % AL Vistasoft Copyright 2014
    
    properties
        filmDistance;   % Probably will go away
        filmDiag;       % Probably will go away, or we will call this camera
    end
    
    methods
        
        %default constructor
        function obj = filmObject(inFilmDistance, inFilmDiag)
                        
            if (ieNotDefined('inFilmDistance')), obj.filmDistance = 140;
            else                                 obj.filmDistance = inFilmDistance;
            end
            
            if (ieNotDefined('inFilmDiag')),     obj.filmDiag = 43.267;
            else                                 obj.filmDiag = inFilmDiag;
            end
            
            
        end
    end
    
end