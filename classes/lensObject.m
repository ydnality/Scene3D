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
        diffractionEnabled;
    end
    
    methods
        
        %default constructor
        function obj = lensObject( aperture, focalLength, center, diffractionEnabled)
            
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
            
             % Units are mm
            %TODO: error checking
            if (ieNotDefined('diffractionEnabled')), obj.diffractionEnabled = false;
            else                           obj.diffractionEnabled = diffractionEnabled;
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

        
        %Performs the Heisenburg Uncertainty Ray Bending method on the rays, given
        %a circular aperture radius, and lens intersection position
        %This function accepts both vector forms of inputs, or individual inputs
        function obj = rayTraceHURB(obj, rays, lensIntersectPosition, curApertureRadius)
            
            %potentially vectorize later for speed
            for i = 1:size(rays.direction, 1)
                %calculate the distance of the intersect point to the center of the lens
                curLensIntersectPosition = lensIntersectPosition(i, :);
                curRay.origin = rays.origin(i, :);
                curRay.direction = rays.direction(i, :);
                curRay.wavelength = rays.wavelength(i);

                ipLength = norm(curLensIntersectPosition);

                %calculate directionS and orthogonal directionL
                directionS = [curLensIntersectPosition(1) curLensIntersectPosition(2) 0];
                directionL = [-curLensIntersectPosition(2) curLensIntersectPosition(1) 0];
                directionS = directionS./norm(directionS);
                directionL = directionL./norm(directionL);

                pointToEdgeS = curApertureRadius - ipLength;   %this is 'a' from paper  //pointToEdgeS stands for point to edge short
                pointToEdgeL = sqrt((curApertureRadius* curApertureRadius) - ipLength * ipLength);  %pointToEdgeS stands for point to edge long

                lambda = curRay.wavelength * 1e-9;  %this converts lambda to meters
                sigmaS = atan(1/(2 * pointToEdgeS *.001 * 2 * pi/lambda));  %the .001 converts mm to m
                sigmaL = atan(1/(2 * pointToEdgeL * .001 * 2 * pi/lambda));

                %this function regenerates a 2D gaussian sample and
                %returns it randOut
                %gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    %experiment for now
                [randOut] = randn(1,2) .* [sigmaS sigmaL];

                %calculate component of these vectors based on 2 random degrees
                %assign noise in the s and l directions according to data at these pointers
                noiseS = randOut(1);
                noiseL = randOut(2);

                %project the original ray (in world coordinates) onto a new set of basis vectors in the s and l directions
                projS = (curRay.direction(1) * directionS(1) + curRay.direction(2) * directionS(2))/sqrt(directionS(1) * directionS(1) + directionS(2) * directionS(2));
                projL = (curRay.direction(1) * directionL(1) + curRay.direction(2) * directionL(2))/sqrt(directionL(1) * directionL(1) + directionL(2) * directionL(2));
                thetaA = atan(projS/curRay.direction(3));   %azimuth - this corresponds to sigmaS
                thetaE = atan(projL/sqrt(projS*projS + curRay.direction(3)* curRay.direction(3)));   %elevation - this corresponds to sigmaL

                %add uncertainty
                thetaA = thetaA + noiseS;
                thetaE = thetaE + noiseL;

                %convert angles back into cartesian coordinates, but in s,l space
                newprojL = sin(thetaE);
                smallH = cos(thetaE);   %smallH corresponds to the projection of the ray onto the s-z plane
                newprojS = smallH * sin(thetaA);
                curRay.direction(3) = smallH * cos(thetaA);

                %convert from s-l space back to x-y space
                curRay.direction(1) = (directionS(1) * newprojS + directionL(1) * newprojL)/sqrt(directionS(1) * directionS(1) + directionL(1) * directionL(1));
                curRay.direction(2) = (directionS(2) * newprojS + directionL(2) * newprojL)/sqrt(directionS(2) * directionS(2) + directionL(2) * directionL(2));
                curRay.direction = curRay.direction./norm(curRay.direction);

                %reassign ray
                rays.origin(i,:) = curRay.origin;
                rays.direction(i, :) = curRay.direction;
                rays.wavelength(i,:) = curRay.wavelength;
            end
        end
    end
    
end