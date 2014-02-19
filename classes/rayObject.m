classdef rayObject <  handle
    % Create a ray object
    %
    % ray = rayObject(origin,direction,wavelength)
    %
    % Example:
    %   rayObject
    %
    % AL Vistasoft Copyright 2014
    
    properties
        origin;
        direction;
        wavelength;
    end
    
    methods
        
        %default constructor
        function obj = rayObject(origin, direction, wavelength)
            if (ieNotDefined('origin')),   obj.origin = [0,0,0];
            else                           obj.origin = origin;
            end
            
            if (ieNotDefined('direction')), obj.direction = [0,0,1]; 
            else                            obj.direction = direction;
            end
            
            if ieNotDefined('wavelength'),  obj.wavelength = 550;
            else                            obj.wavelength = wavelength;
            end
        end
        
    end
    
end