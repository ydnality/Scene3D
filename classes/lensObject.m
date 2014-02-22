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
        elementArray;
        totalOffset;
        numEls;
        apertureRadius;
        apertureSample;
        centerPosition;
    end
    
    methods
        
        %default constructor
        function obj = lensObject(elOffset, elRadius, elAperture, elN, aperture, center)
            
            % Units are mm
            if (ieNotDefined('elOffset')), elOffset = 0;
            end
            
            % Units are mm
            if (ieNotDefined('elRadius')), elRadius = 10;
            end
            
            % Units are mm
            if (ieNotDefined('elAperture')), elAperture = 10;
            end
            
            % Units are mm
            if (ieNotDefined('elN')), elN = 1;
            end
           
            % Units are mm
            if (ieNotDefined('aperture')), obj.apertureRadius = 3;
            else                           obj.apertureRadius = aperture;
            end            
            
            
            % Units are mm
            %TODO: error checking
            if (ieNotDefined('center')), obj.centerPosition = [0 0 -1.5];
            else                           obj.centerPosition = center;
            end 
            
            obj.numEls = length(elOffset); % we must update numEls each time we add a lens element
            
            %error checking
            if (obj.numEls~= length(elRadius) || obj.numEls~= length(elAperture) || obj.numEls ~= length(elN))
                error('input vectors must all be of the same lengths');
            end
            
            obj.elementArray = lensElementObject();
            for i = 1:length(elOffset)
                obj.elementArray(i) = lensElementObject(elOffset(i), elRadius(i), elAperture(i), elN(i));
            end
            
            obj.computeCenters();
            obj.calculateApertureSample();
        end
        
        %Calculates the total offset of the lens by adding all existing
        %offsets
        %make this private later
        function computeTotalOffset(obj)
            obj.totalOffset  = 0;
            for i = 1:obj.numEls
                obj.totalOffset = obj.totalOffset + obj.elementArray(i).offset;
            end
        end

        %computes the spherical centers of each element
        function computeCenters(obj)

            obj.computeTotalOffset();
            prevSurfaceZ = -obj.totalOffset;
            for i = length(obj.elementArray):-1:1
                obj.elementArray(i).zIntercept = prevSurfaceZ + obj.elementArray(i).offset;  %TODO: these will change later with sets
                obj.elementArray(i).sphereCenter = [0 0 obj.elementArray(i).zIntercept+ obj.elementArray(i).radius];
                prevSurfaceZ = obj.elementArray(i).zIntercept;
            end
        end
        
        function obj = calculateApertureSample(obj)
           
            %loop through aperture positions and uniformly sample the aperture
            %everything is done in vector form for speed
            [rectApertureSample.X, rectApertureSample.Y] = meshgrid(linspace(-1, 1, 3),linspace(-1, 1, 3)); %adjust this if needed - this determines the number of samples per light source
            
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
    end
    
end