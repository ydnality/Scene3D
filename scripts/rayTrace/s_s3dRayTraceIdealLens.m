%% ray-tracing for an ideal lens


diffractionEnabled = true;

% point sources

% declare point sources in world space.  The camera is usually at [0 0 0],
% and pointing towards -z.  We are using a right-handed coordinate system.
% pointSources = [0 0 -100;
%                 0 50 -100;
%                 0 -50 -100;
%                 50 0 -100;
%                 50 50 -100;
%                 50 -50 -100;
%                 -50 0 -100;
%                 -50 50 -100;
%                 -50 -50 -100;];

% [XGrid YGrid] = meshgrid(-4000:1000:4000,-4000:1000:4000);
% pointSources = [XGrid(:) YGrid(:) ones(size(XGrid(:))) * -20000];
pointSources = [0 0 -20000];
 
% wavelength sampling
wave = 400:10:700;  % in nm
wavelengthConversion = [wave' (1:31)'];  %[475 3; 550 2; 600 1];

% sensor properties
sensor.size = [1 1]; 
%sensor.size = [48 48]; %in mm
sensor.position = [0 0 50]; 
sensor.resolution = [200 200 length(wave)];
sensor.image = zeros(sensor.resolution);
sensor.focalLength = 50; %in mm

apertureRadius =.1; % in mm
%loop through aperture positions and uniformly sample the aperture
%everything is done in vector form for speed
[apertureSample.X, apertureSample.Y] = meshgrid(linspace(-1, 1, 1000),linspace(-1, 1, 2000)); %adjust this if needed



% this position should NOT change
lensCenterPosition = [0 0 0];

%loop through all point sources
for curInd = 1:size(pointSources, 1);
    
    curPointSource = pointSources(curInd, :);
    
    %trace ray from point source to lens center, to image.  This helps
    %determine the point of focus
    centerRay.origin = curPointSource;
    centerRay.direction = lensCenterPosition - centerRay.origin;
    centerRay.direction = centerRay.direction./norm(centerRay.direction);
    
    %calculate the in-focus position using thin lens equation
    inFocusDistance = 1/(1/sensor.focalLength - -1/curPointSource(3));
    inFocusT = (inFocusDistance - centerRay.origin(3))/centerRay.direction(3);
    inFocusPosition = centerRay.origin + inFocusT .* centerRay.direction;
    

    
    %assume a circular aperture, and make a mask that is 1 when the pixel
    %is within a circle of radius 1
    apertureMask = (apertureSample.X.^2 + apertureSample.Y.^2) <= 1;
    
    scaledApertureSample.X = apertureSample.X .* apertureRadius; 
    scaledApertureSample.Y = apertureSample.Y .* apertureRadius; 
    
    %remove cropped sections of aperture
    croppedApertureSample.X =  scaledApertureSample.X;
    croppedApertureSample.X(apertureMask == 0) = [];
    croppedApertureSample.X = croppedApertureSample.X';
    croppedApertureSample.Y =  scaledApertureSample.Y;
    croppedApertureSample.Y(apertureMask == 0) = [];
    croppedApertureSample.Y = croppedApertureSample.Y';
    
    %calculate the origin and direction of the rays
    rays.origin = repmat(curPointSource, [size(croppedApertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
    rays.direction = [(croppedApertureSample.X(:) -  rays.origin(:,1)) (croppedApertureSample.Y(:) -  rays.origin(:,2)) (lensCenterPosition(3) - rays.origin (:,3)) .* ones(size(croppedApertureSample.Y(:)))];
    rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction
       
    
    %add different wavelengths
    %first duplicate the existing entries, and create one for each
    %wavelength
    subLength = size(rays.origin, 1);
    rays.origin = repmat(rays.origin, [length(wave) 1]);
    rays.direction = repmat(rays.direction, [length(wave) 1]);
    %creates a vector representing wavelengths... for example: [400 400 400... 410 410 410... ..... 700]
    rays.wavelength = vectorize((wave' * ones(1, subLength))');  
    

    %when intersecting ideal lens, change the direction to intersect the
    %inFocusPosition, and update the origin
    lensIntersectT = (lensCenterPosition(3) - rays.origin(:,3))./ rays.direction(:,3);
    lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
    %calculate new direction
    newRays.origin = lensIntersectPosition;
    newRays.direction = repmat(inFocusPosition , [size(newRays.origin,1) 1 ]) - newRays.origin;
    newRays.wavelength = rays.wavelength;
    
    
    % -- diffraction -- (make this into a function)
    if (diffractionEnabled)
        %vectorize this later for speed
        for i = 1:size(rays.direction, 1)
            %calculate the distance of the intersect point to the center of the lens
            curLensIntersectPosition = lensIntersectPosition(i, :);
            curRay.origin = newRays.origin(i, :);
            curRay.direction = newRays.direction(i, :);
            curRay.wavelength = newRays.wavelength(i);
            
            ipLength = norm(curLensIntersectPosition);
            
            %calculate directionS and orthogonal directionL
            directionS = [curLensIntersectPosition(1) curLensIntersectPosition(2) 0];
            directionL = [-curLensIntersectPosition(2) curLensIntersectPosition(1) 0];
            directionS = directionS./norm(directionS);
            directionL = directionL./norm(directionL);
            
            pointToEdgeS = apertureRadius- ipLength;   %this is 'a' from paper  //pointToEdgeS stands for point to edge short
            pointToEdgeL = sqrt((apertureRadius* apertureRadius) - ipLength * ipLength);  %pointToEdgeS stands for point to edge long
            
            lambda = curRay.wavelength * 1e-9;  %this converts lambda to meters
            sigmaS = atan(1/(2 * pointToEdgeS *.001 * 2 * pi/lambda));  %the .001 converts mm to m
            sigmaL = atan(1/(2 * pointToEdgeL * .001 * 2 * pi/lambda));
            
            %this function regenerates a 2D gaussian sample and
            %returns it randOut
            %gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    %experiment for now
            [randOut] = randn(1,2) .* [sigmaS sigmaL];   
            
            %calculate component of these vectors based on 2 random degrees
            %assign noise in the s and l directions according to data at these pointers
            %             noiseS = (float)(*noiseSPointer);
            %             noiseL = (float)(*noiseLPointer);
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
            newRays.origin(i,:) = curRay.origin;
            newRays.direction(i, :) = curRay.direction;
            newRays.wavelength(i,:) = curRay.wavelength;
        end
    end
    % -- end diffraction --
    
    
    %calculate intersection point at sensor
    intersectZ = repmat(sensor.position(3), [size(newRays.origin, 1) 1]);
    intersectT = (intersectZ - newRays.origin(:, 3))./newRays.direction(:, 3);
    intersectPosition = newRays.origin + newRays.direction .* repmat(intersectT, [1 3]);
    
     %imagePixel is the pixel that will gain a photon due to the traced ray
    imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)]; 
    imagePixel.position = real(imagePixel.position); %add error handling for this
    imagePixel.position = round(imagePixel.position * sensor.resolution(1)/sensor.size(1) + ...
        repmat( sensor.resolution(1:2)./2, [size(imagePixel.position,1) 1]));   %
    %scale the position to a sensor position
    imagePixel.position(imagePixel.position < 1) = 1; %make sure pixel is in range
    imagePixel.position = min(imagePixel.position, repmat(sensor.resolution(1:2), [size(imagePixel.position,1) 1]));
    imagePixel.wavelength = newRays.wavelength; 
    
    
    
    %add a value to the intersection position
    for i = 1:size(rays.origin , 1)
        wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(wavelengthConversion == imagePixel.wavelength(i))];  %pixel to update
        sensor.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) =  sensor.image(wantedPixel(1), wantedPixel(2), wantedPixel(3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
        %illustrations for debugging
%         line([newRays.origin(i, 3) intersectPosition(i, 3)] ,  [newRays.origin(i, 2);  intersectPosition(i, 2)]);
    end
end

% figure; imshow(sensor.image/ max(sensor.image(:)));


%assign as optical image
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',sensor.image);
vcAddAndSelectObject(oi); oiWindow;


