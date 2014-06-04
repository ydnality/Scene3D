classdef lensObject <  handle
    % Create a lens object
    %
    %   lens = lensObject( parameter, value, ....);
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
    % Is it possible to write this in the form
    %  lensObject(lensType, relevant parameter list ... )
    %
    % Examples:
    %   lens = lensObject; lens.apertureSample = [51,51];
    %   fGrid = lens.fullGrid;
    %   aMask = lens.apertureMask; vcNewGraphWin; imagesc(aMask)
    %   aGrid = lens.apertureGrid; vcNewGraphWin; plot(aGrid.X(:),aGrid.Y(:),'o')
    %
    %   pointSource = [0 0 -50];
    %   rays = lens.rtSourceToEntrance(pointSource);
    %   ro =  rays.origin;
    %   rd =  rays.origin + rays.direction;
    %   vcNewGraphWin; hold on;
    %   for ii=1:10:size(rd,1), line([ro(ii,1),rd(ii,1)],[ro(ii,2),rd(ii,2)],[ro(ii,3),rd(ii,3)]); end
    %
    % AL Vistasoft Copyright 2014
    
    properties
        name = 'default';
        type = 'lens';
        apertureRadius = 1;        % mm
        
        % apertureGrid = 
        % [X,Y] =
        % meshgrid([-apertureRadius:apertureSample(1):apertureRadius],
        % ...);  
        % Then eliminate outside the aperture
        apertureSample = [11 11];  % Number of spatial samples in the aperture.  Use odd number
        focalLength = 50;          % mm
        centerPosition = [0 0 0];  % almost always on axis
        diffractionEnabled = false;% On/off
        wave = 400:50:700;         % nm
    end
    
    methods
        
        % %%%%% Lens object constructor %%%%%%%
        function obj = lensObject(varargin)

            for ii=1:2:length(varargin)
                switch varargin{ii}
                    case 'apertureRadius'
                        % Units are mm
                        obj.apertureRadius = varargin{ii+1};

                    case 'focalLength'
                        obj.focalLength = varargin{ii+1};
                        
                    case 'centerPosition'
                        obj.centerPosition = varargin{ii+1};
                        
                    case 'diffractionEnabled'
                        obj.diffractionEnabled = varargin{ii+1};
                        
                    case 'wave'
                        obj.wave = varargin{ii+1};
                        
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
        end
        
        function res = get(obj,paramName,varargin)
            % Get various derived lens properties though this call
            res = lensGet(obj,paramName,varargin);
        end
        
        
        function apertureMask = apertureMask(obj)
            % Creates a sampling grid on the circular aperture.
            %
            % We first build the samples on a unit radius disk, and then
            % scale to the aperture radius.
            %
            % The returned pattern is placed in apertureSamples.
            % These are a set of (x,y) positions.
            %
            % sampleResolution is a 2d vector. It contains the number of
            % samples in the rectangular grid before the samples are
            % cropped by the aperture shape
            %
            aGrid = obj.fullGrid;
            
            % We assume a circular aperture. This is a mask that is 1 when
            % the pixel is within a circle of radius
            apertureMask = (aGrid.X.^2 + aGrid.Y.^2) <= obj.apertureRadius^2;
            % vcNewGraphWin;  mesh(double(apertureMask))
            
        end
        
        function aGrid = fullGrid(obj,randJitter)
            
            if (ieNotDefined('randJitter')), randJitter = false; end
            
            % First make the rectangular samples.
            xSamples = linspace(-obj.apertureRadius, obj.apertureRadius, obj.apertureSample(1));
            ySamples = linspace(-obj.apertureRadius, obj.apertureRadius, obj.apertureSample(2));
            [X, Y] = meshgrid(xSamples,ySamples);
            
            %Add a random jitter.  Uniform distribution.  Scales to plus or
            %minus half the sample width
            if(randJitter)
                X = X + (rand(size(X)) - .5) * (xSamples(2) - xSamples(1));
                Y = Y + (rand(size(Y)) - .5) * (ySamples(2) - ySamples(2));
            end
            aGrid.X = X; aGrid.Y = Y;
            
        end
        
        function aGrid = apertureGrid(obj,randJitter)
            % Find the grid, mask it and return only the (X,Y) positions
            % inside the masked region.  This is the usual set of positions
            % that you want.
            if ieNotDefined('randJitter'), randJitter = false; end
            
            aGrid = fullGrid(obj,randJitter);
            aMask = apertureMask(obj);
            aGrid.X = aGrid.X(aMask);
            aGrid.Y = aGrid.Y(aMask);
        end
        
        function rays = rtSourceToEntrance(obj, pointSource)
            % Ray trace from a point to the aperture grid 
            % 
            % See also: rtEntranceToExit in meLens
            
            rays = rayObject;
            
            % We create rays from the point source to each aperture grid point
            aGrid = obj.apertureGrid;
            nPts  = numel(aGrid.X(:));
            
            %the new origin is the position of the current point source
            rays.origin = repmat(pointSource, [nPts, 1, 1] );
            
            aGrid.Z = repmat(obj.centerPosition(3),[nPts,1]);
            ePoints = [aGrid.X(:),aGrid.Y(:),aGrid.Z(:)];
            
            rays.direction = rayDirection(rays.origin,ePoints);
                        
            % rays.direction = [(aGrid.X(:) -  rays.origin(:,1)) (aGrid.Y(:) -  rays.origin(:,2)) (obj.centerPosition(3) - rays.origin (:,3)) .* ones(size(obj.apertureSample.Y(:)))];
            % rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction
        end
                      
        
        function obj = rtHURB(obj, rays, lensIntersectPosition, curApertureRadius)
            %Performs the Heisenburg Uncertainty Ray Bending method on the
            %rays, given a circular aperture radius, and lens intersection
            %position This function accepts both vector forms of inputs, or
            %individual inputs
            
            %potentially vectorize later for speed
%             for i = 1:size(rays.direction, 1)
                %calculate the distance of the intersect point to the center of the lens
                
%                 curLensIntersectPosition = lensIntersectPosition(i, :);
                
                %we don't care about the z coordinate so we remove it
%                 curLensIntersectPosition = curLensIntersectPosition(1:2);
                
%                 curRay.origin = rays.origin(i, :);
%                 curRay.direction = rays.direction(i, :);
%                 curRay.wavelength = rays.wavelength(i);
                
                ipLength = sqrt(sum(dot(lensIntersectPosition(:, (1:2)), lensIntersectPosition(:, (1:2)), 2), 2));
                
                %calculate directionS and orthogonal directionL
                
                directionS = [lensIntersectPosition(:, 1) lensIntersectPosition(:,2) zeros(length(lensIntersectPosition), 1)];
                directionL = [-lensIntersectPosition(:,2) lensIntersectPosition(:,1) zeros(length(lensIntersectPosition), 1)];
                
                
                
                normS = repmat(sqrt(sum(dot(directionS, directionS, 2), 2)), [1 3]);
                normL = repmat(sqrt(sum(dot(directionL, directionL, 2), 2)), [1 3]);
                divideByZero = sum(normS,2) ==0;
                
                
                directionS(~divideByZero, :) = directionS(~divideByZero, :)./normS(~divideByZero, :);
                directionL(~divideByZero, :) = directionL(~divideByZero, :)./normL(~divideByZero, :);
                directionS(divideByZero, :) = [ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 2)];
                directionL(divideByZero, :) = [zeros(sum(divideByZero == 1), 1) ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 1)];
                
%                 if (norm(directionS)~= 0)
%                     directionS = directionS./norm(directionS);
%                     directionL = directionL./norm(directionL);
%                 else
%                     %this is the case that the ray hits the absolute center
%                     %Then, there is a special case to avoid division by 0
%                     directionS = [1 0 0];
%                     directionL = [0 1 0];
%                 end
                
                pointToEdgeS = curApertureRadius - ipLength;   %this is 'a' from paper  //pointToEdgeS stands for point to edge short
                pointToEdgeL = sqrt((curApertureRadius* curApertureRadius) - ipLength .* ipLength);  %pointToEdgeS stands for point to edge long
                
                lambda = rays.wavelength * 1e-9;  %this converts lambda to meters
                sigmaS = atan(1./(2 * pointToEdgeS *.001 * 2 * pi./lambda));  %the .001 converts mm to m
                sigmaL = atan(1./(2 * pointToEdgeL * .001 * 2 * pi./lambda));
                
                %this function regenerates a 2D gaussian sample and
                %returns it randOut
                %gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    %experiment for now
                [randOut] = randn(length(sigmaS),2) .* [sigmaS sigmaL];
                
                %calculate component of these vectors based on 2 random degrees
                %assign noise in the s and l directions according to data at these pointers
                noiseS = randOut(:,1);
                noiseL = randOut(:,2);
                
                %project the original ray (in world coordinates) onto a new set of basis vectors in the s and l directions
                projS = (rays.direction(: , 1) .* directionS(: ,1) + rays.direction(: , 2) .* directionS(:,2))./sqrt(directionS(:,1) .* directionS(:,1) + directionS(:,2) .* directionS(:,2));
                projL = (rays.direction(: , 1) .* directionL(:, 1) + rays.direction(: , 2 ) .* directionL(:,2))./sqrt(directionL(:,1) .* directionL(:,1) + directionL(:,2) .* directionL(:,2));
                thetaA = atan(projS./rays.direction(: , 3));   %azimuth - this corresponds to sigmaS
                thetaE = atan(projL./sqrt(projS.*projS + rays.direction(:, 3).* rays.direction(: , 3)));   %elevation - this corresponds to sigmaL
                
                %add uncertainty
                thetaA = thetaA + noiseS;
                thetaE = thetaE + noiseL;
                
                %convert angles back into cartesian coordinates, but in s,l space
                newprojL = sin(thetaE);
                smallH = cos(thetaE);   %smallH corresponds to the projection of the ray onto the s-z plane
                newprojS = smallH .* sin(thetaA);
                rays.direction(:, 3) = smallH .* cos(thetaA);
                
                %convert from s-l space back to x-y space
                rays.direction(:, 1) = (directionS(:, 1) .* newprojS + directionL(:, 1) .* newprojL)./sqrt(directionS(:,1) .* directionS(:,1) + directionL(:,1) .* directionL(:,1));
                rays.direction(:, 2) = (directionS(:, 2) .* newprojS + directionL(:, 2) .* newprojL)./sqrt(directionS(:,2) .* directionS(:,2) + directionL(:,2) .* directionL(:,2));
                normDirection = repmat(sqrt(sum(dot(rays.direction, rays.direction, 2),2)), [1 3]);
                rays.direction = rays.direction./normDirection;
                
                %reassign ray
%                 rays.origin(i,:) = curRay.origin;
%                 rays.direction(i, :) = curRay.direction;
%                 rays.wavelength(i,:) = curRay.wavelength;
%             end
        end
        
    end
    
end