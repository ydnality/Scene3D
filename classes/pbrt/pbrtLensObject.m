classdef pbrtLensObject < handle
    % this is the simplest lens - no aperture is required because pinhole
    % cameras have no aperture and the pinhole lens will inherit this
    % superclas. This will be a superclass that will be inherited by other
    % classes in the future
    
    properties
        filmDistance;     % Film to lens distance in millimeters
        filmDiag;         % Film diagonal in millimeters
    end
    
    methods
        
        %default constructor
        function obj = pbrtLensObject(inFilmDistance, inFilmDiag)
            % Could switch to 'param',val with varargin here, some day.
            if (ieNotDefined('inFilmDistance'))
                % Example lens
                obj.filmDistance = 140;
            else
                obj.filmDistance = inFilmDistance;
            end
            
            if (ieNotDefined('inFilmDiag'))
                % Example lens
                obj.filmDiag = 43.267;
            else
                obj.filmDiag = inFilmDiag;
            end
            
        end
        
    end
end