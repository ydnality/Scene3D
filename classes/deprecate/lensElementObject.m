classdef lensElementObject <  handle
    % Create a lens element object
    %
    %   lens = lensElementObject(offset, radius, aperture, n);  % Units are mm
    %
    % Presently we only represent spherical lenses and apertures.
    %
    % These are defined by a series of surfaces. We code the offset to each
    % surface radius to the center of the spherical lens.  Positive means
    % the sphere center is to the right. Aperture parameters (a single number is a
    % diameter in mm). index of refraction (n) for the material to the left
    % of the surface.
    % 
    % 
    % Example:
    %   lensElementObject
    %   lensElementObject(30,250)
    %
    % AL Vistasoft Copyright 2014
    
    properties  %eventually this will be set access = private
        offset;
        radius;
        aperture;
        n;
        sphereCenter;
        zIntercept;
    end
    
    methods
        
        %default constructor
        function obj = lensElementObject(offset, radius, aperture, n)
            
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

            %initial value for sphere center
            obj.sphereCenter = 10;
            
            obj.zIntercept = 0;
        end

    end
    
end