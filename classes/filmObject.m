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
        position;
        size;  %in mm
        wave;
        waveConversion;
        resolution;  %decide whether this includes the 3rd dimension or not
        image;
    end
    
    methods
        
        %default constructor
        function obj = filmObject(position, size,  wave, waveConversion, resolution)
                        
            if (ieNotDefined('position')), obj.position = [0 0 100];
            else                           obj.position = position;
            end
            
            if (ieNotDefined('size')),     obj.size = [48 48];
            else                           obj.size = size;
            end
            
            if (ieNotDefined('wave')),     obj.wave = [400 550 700];  % in nm;
            else                           obj.wave = wave;
            end      
            
            %this field might go away soon
            if (ieNotDefined('waveConversion')),     obj.waveConversion = [400 1; 550 2; 700 3];  % in nm;
            else                           obj.waveConversion = waveConversion;
            end  
            
            if (ieNotDefined('resolution')),     obj.resolution = [200 200 length(obj.wave)];
            else                           obj.size = resolution;
            end
            %TODO: error checking.  make sure all dimensions are good
            
            obj.image = zeros(obj.resolution);
        end
        
        
    end
    
end