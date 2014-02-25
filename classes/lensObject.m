classdef lensObject <  handle
    % Create a lens object
    %
    %   lens = lensObject( aperture, focalLength);  % Units are mm
    %
    % Presently we only represent spherical lenses and apertures.
    %
    % This is meant to be the superclass for additional lens objects.  This
    % is not meant to be a stand-alone object.  This object contains some
    % basic properties common to almost all lenses.  Also - it contains the
    % calculateApertureSamples function which outputs samples on the
    % aperture to aid in ray-tracing. 
    %
    % TODO: determine which aperture is it best to sample?
    % 
    % AL Vistasoft Copyright 2014
    
    properties
        apertureRadius;
        apertureSample;
        focalLength;
        centerPosition;
    end
    
    methods
        
        %default constructor
        function obj = lensObject( aperture, focalLength, center)
            
            % Units are mm
            if (ieNotDefined('aperture')), obj.apertureRadius = 3;
            else                           obj.apertureRadius = aperture;
            end            
            
            % Units are mm
            if (ieNotDefined('focalLength')), obj.focalLength = 50;
            else                           obj.focalLength = focalLength;
            end          
            
             % Units are mm
            %TODO: error checking
            if (ieNotDefined('center')), obj.centerPosition = [0 0 0];
            else                           obj.centerPosition = center;
            end 
            
            obj.calculateApertureSample();
        end
        

        %creates a uniform circular sampling patern on the aperture.
        %Sample resolution: the resolution of the rectangular grid before
        %it is cropped by the aperture shape
        function obj = calculateApertureSample(obj, sampleResolution)

            % Units are mm
            if (ieNotDefined('sampleResolution')), sampleResolution = [3 3];
            end      
            
            %loop through aperture positions and uniformly sample the aperture
            %everything is done in vector form for speed
            [rectApertureSample.X, rectApertureSample.Y] = meshgrid(linspace(-1, 1, sampleResolution(1)),linspace(-1, 1, sampleResolution(2))); %adjust this if needed - this determines the number of samples per light source
            
            %assume a circular aperture, and make a mask that is 1 when the pixel
            %is within a circle of radius 1
            apertureMask = (rectApertureSample.X.^2 + rectApertureSample.Y.^2) <= 1;
            scaledApertureSample.X = rectApertureSample.X .* obj.apertureRadius;
            scaledApertureSample.Y = rectApertureSample.Y .* obj.apertureRadius;
            
            %remove cropped sections of aperture
            croppedApertureSample.X =  scaledApertureSample.X;
            croppedApertureSample.X(apertureMask == 0) = [];
            croppedApertureSample.X = croppedApertureSample.X';
            croppedApertureSample.Y =  scaledApertureSample.Y;
            croppedApertureSample.Y(apertureMask == 0) = [];
            croppedApertureSample.Y = croppedApertureSample.Y';
            
            obj.apertureSample = croppedApertureSample;
        end
        
        
        %traces rays from a point source to a sampling function on the lens
        function obj = rayTraceSourceToLens(obj, curPointSource, rays)
            rays.origin = repmat(curPointSource, [size(obj.apertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
            rays.direction = [(obj.apertureSample.X(:) -  rays.origin(:,1)) (obj.apertureSample.Y(:) -  rays.origin(:,2)) (obj.centerPosition(3) - rays.origin (:,3)) .* ones(size(obj.apertureSample.Y(:)))];
            rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction
        end
        
%         function obj =  rayTraceThroughLens(obj, rays)  - TODO: make this
%         into an abstract function for all lenses
        
        
    end
    
end