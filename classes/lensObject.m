classdef lensObject <  handle
    % Create a lens object
    %
    %   lens = lensObject(filmDistance,filmDiag);  % Units are mm
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
        offset;
        radius;
        aperture;
        n;
    end
    
    methods
        
        %default constructor
        function obj = lensObject(offset, radius, aperture, n)
            
            % Units are mm
            if (ieNotDefined('offset')), obj.offset = 0;
            else                         obj.offset = offset;
            end
            
            % Units are mm
            if (ieNotDefined('radius')), obj.radius = 10;
            else                         obj.radius = radius;
            end
            
            % Units are mm
            if (ieNotDefined('aperture')), obj.aperture = 10;
            else                           obj.aperture = aperture;
            end
            
            % Units are mm
            if (ieNotDefined('n')), obj.n = 1;
            else                    obj.n = n;
            end
            
            
        end
    end
    
end