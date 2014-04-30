classdef ppsfObject < rayObject
    % Create a ray object
    %
    % ray = rayObject(origin,direction,wavelength)
    %
    % Example:
    %   rayObject
    %
    % AL Vistasoft Copyright 2014
    
    properties
%         origin;
%         direction;
%         wavelength;
%         waveIndex;
          pointSourceDepth;          %depth of the originating point source
          pointSourceFieldHeight;    %field height of the originating point source
          outsideApertureLocation;   %where rays intersect further most lens surface
          apertureLocation;    %where rays intersect the actual lens aperture (in the middle usually)
    end
    
    methods
        function obj = ppsfObject(origin, direction, wavelength, pointSourceDepth, pointSourceFieldHeight)
            %default constructor
            
            % rayObject properties
            if (ieNotDefined('origin')),   obj.origin = [0,0,0];
            else                           obj.origin = origin;
            end
            
            if (ieNotDefined('direction')), obj.direction = [0,0,1]; 
            else                            obj.direction = direction;
            end
            
            if ieNotDefined('wavelength'),  obj.wavelength = 550;
            else                            obj.wavelength = wavelength;
            end
            
            %new properties
            if ieNotDefined('pointSourceDepth'),  obj.pointSourceDepth = 100;
            else                            obj.pointSourceDepth = pointSourceDepth;
            end
            
            if ieNotDefined('pointSourceFieldHeight'),  obj.pointSourceFieldHeight = 0;
            else                            obj.pointSourceFieldHeight = pointSourceFieldHeight;
            end
        end
        
        function setApertureLocation(obj, apertureLocation)
            %setApertureLocation(obj, apertureLocation)
            %
            %sets the intersection of rays with actual lens aperture
            %location
            %TODO: error handling
            
           obj.apertureLocation = apertureLocation;
           return;
        end
    end
    
end