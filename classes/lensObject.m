% this is the simplest lens - no aperture is required because pinhole cameras
% have no aperture and the pinhole lens will inherit this superclas.  
% This will be a superclass that will be inherited by other classes in the future
classdef lensObject <  handle
    properties 
        filmDistance;
        filmDiag;
    end
    methods
        
        %default constructor
        function obj = lensObject(inFilmDistance, inFilmDiag)
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