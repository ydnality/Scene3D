classdef rayTracerObject <  handle
    % Create a ray object
    %
    % ray = rayObject(origin,direction,wavelength)
    %
    % Example:
    %   rayObject
    %
    % AL Vistasoft Copyright 2014
    
    properties
        %think of some things to put in here
    end
    
    methods
        
        %default constructor
        function obj = rayTracerObject()
            
        end
        
        function rays = traceLightToLens(obj, curPointSource, lens)
            rays.origin = repmat(curPointSource, [size(lens.apertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
            rays.direction = [(lens.apertureSample.X(:) -  rays.origin(:,1)) (lens.apertureSample.Y(:) -  rays.origin(:,2)) (lensCenterPosition(3) - rays.origin (:,3)) .* ones(size(lens.apertureSample.Y(:)))];
            rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction

        end
        
        
        function traceThroughLens(obj, lens)
            
        end
        
        function traceToFilm(obj, film)
            
        end
    end
    
end