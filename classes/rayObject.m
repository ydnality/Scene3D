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
        
        
        function obj = rayObject(origin, direction, wavelength)
            %default constructor
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
        
        
        function obj = traceSourceToLens(obj, curPointSource, lens)
            %traces rays from a point source to a sampling function on the lens
            obj.origin = repmat(curPointSource, [size(lens.apertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
            obj.direction = [(lens.apertureSample.X(:) -  obj.origin(:,1)) (lens.apertureSample.Y(:) -  obj.origin(:,2)) (lens.centerPosition(3) - obj.origin (:,3)) .* ones(size(lens.apertureSample.Y(:)))];
            obj.direction = obj.direction./repmat( sqrt(obj.direction(:, 1).^2 + obj.direction(:, 2).^2 + obj.direction(:,3).^2), [1 3]); %normalize direction
        end
        

        function obj =  traceThroughLens(obj, lens)
            %performs ray-trace of the lens, given an input bundle or rays
            %outputs the rays that have been refracted by the lens
            %TODO: consdier moving this to the lens
            prevN = 1;  %assume that we start off in air
            
            %initialize newRays to be the old ray.  We will update it later.
%             newRays = rays;
            
            prevSurfaceZ = -lens.totalOffset;
            
            for lensEl = lens.numEls:-1:1
                curEl = lens.elementArray(lensEl);
%                 curEl.center = [0 0 prevSurfaceZ + curEl.offset + curEl.radius];
                
                %illustrations for debug
                zPlot = linspace(curEl.sphereCenter(3) - curEl.radius, curEl.sphereCenter(3) + curEl.radius, 10000);
                yPlot = sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                yPlotN = -sqrt(curEl.radius^2 - (zPlot - curEl.sphereCenter(3)) .^2);
                arcZone = 5;
                %TODO:find a better way to plot the arcs later - this one is prone to potential problem
                withinRange = and(and((yPlot < curEl.aperture),(zPlot < prevSurfaceZ + curEl.offset + arcZone)), (zPlot > prevSurfaceZ + curEl.offset - arcZone));
                line(zPlot(withinRange), yPlot(withinRange));
                line(zPlot(withinRange), yPlotN(withinRange));
                
                %vectorize this operation later
                for i = 1:size(obj.origin, 1)
                    %get the current ray
                    ray.direction = obj.direction(i,:);   %TODO: replace with real ray object
                    ray.origin = obj.origin(i,:);
                    ray.wavelength = obj.wavelength(i);
                    
                    %calculate intersection with spherical lens element
                    radicand = dot(ray.direction, ray.origin - curEl.sphereCenter)^2 - ...
                        ( dot(ray.origin -curEl.sphereCenter, ray.origin -curEl.sphereCenter)) + curEl.radius^2;
                    if (curEl.radius < 0)
                        intersectT = (-dot(ray.direction, ray.origin - curEl.sphereCenter) + sqrt(radicand));
                    else
                        intersectT = (-dot(ray.direction, ray.origin - curEl.sphereCenter) - sqrt(radicand));
                    end
                    
                    %make sure that T is > 0
                    if (intersectT < 0)
                        disp('Warning: intersectT less than 0.  Something went wrong here...');
                    end
                    
                    intersectPosition = ray.origin + intersectT * ray.direction;
                    
                    %illustrations for debugging
%                     lensIllustration(max(round(intersectPosition(2) * 100 + 150),1), max(-round(intersectPosition(3) * 1000), 1)) = 1;  %show a lens illustration
                    hold on;
                    line([ray.origin(3) intersectPosition(3) ], [ray.origin(2) intersectPosition(2)] ,'Color','b','LineWidth',1);
                    
                    normalVec = intersectPosition - curEl.sphereCenter;  %does the polarity of this vector matter? YES
                    normalVec = normalVec./norm(normalVec);
                    if (curEl.radius < 0)  %which is the correct sign convention? This is correct
                        normalVec = -normalVec;
                    end
                    
                    %modify the index of refraction depending on wavelength
                    %TODO: have this be one of the input parameters (N vs. wavelength)
                    if (curEl.n ~= 1)
                        curN = (ray.wavelength - 550) * -.04/(300) + curEl.n;
                    else
                        curN = 1;
                    end
                    

                    ratio = prevN/curN;    %snell's law index of refraction
                    
                    %Vector form of Snell's Law
                    c = -dot(normalVec, ray.direction);
                    newVec = ratio *ray.direction + (ratio*c -sqrt(1 - ratio^2 * (1 - c^2)))  * normalVec;
                    newVec = newVec./norm(newVec); %normalize
                    
                    %update the direction of the ray
                    obj.origin(i, : ) = intersectPosition;
                    obj.direction(i, : ) = newVec;
                end
                prevN = curN;
                
                prevSurfaceZ = prevSurfaceZ + curEl.offset;
            end
        end
        
        
        function obj = recordOnFilm(obj, film)
        %records the ray on the film
          
            %check if no rays - only record if there are any
            if(~isempty(obj.origin))
                
                %calculate intersection point at sensor
                intersectZ = repmat(film.position(3), [size(obj.origin, 1) 1]);
                intersectT = (intersectZ - obj.origin(:, 3))./obj.direction(:, 3);
                intersectPosition = obj.origin + obj.direction .* repmat(intersectT, [1 3]);
                
                %imagePixel is the pixel that will gain a photon due to the traced ray
                imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
                imagePixel.position = real(imagePixel.position); %add error handling for this
                imagePixel.position = round(imagePixel.position * film.resolution(1)/film.size(1) + ...
                    repmat(-film.position(1:2) + film.resolution(1:2)./2, [size(imagePixel.position,1) 1]));   %
               
                %scale the position to a sensor position
                %             imagePixel.position(imagePixel.position < 1) = 1; %make sure pixel is in range
                %             imagePixel.position = real(min(imagePixel.position, repmat(film.resolution(1:2), [size(imagePixel.position,1) 1])));
                
                imagePixel.wavelength = obj.wavelength;
                
                
                %attempted vectorized version - runs out of memory
%                 
%                 %todo - somehow fix this
%                 convertChannel = uint8((imagePixel.wavelength - 400)/10 + 1);
%                 
%                 wantedPixel = [imagePixel.position(:, 1) imagePixel.position(:,2) convertChannel];  %pixel to update
%                 
%                 recordablePixels =and(and(and(wantedPixel(:, 1) >= 1,  wantedPixel(:,1) <= film.resolution(1)), (wantedPixel(:, 2) > 1)), wantedPixel(:, 2) <= film.resolution(2));
%                 
%                 %remove the nonrecordablePixels
%                 wantedPixel = wantedPixel(recordablePixels, :);
%                 
%                 
%                 film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) =  film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;

%                 for i = 1:10:size(obj.origin , 1)
%                     %illustrations for debugging (out of bounds rays will
%                     %still be displayed)
%                     line(real([obj.origin(i, 3) intersectPosition(i, 3)]) ,  real([obj.origin(i, 2);  intersectPosition(i, 2)]));
%                 end

                
                %non-vectorized
                %add a value to the intersection position
                for i = 1:size(obj.origin , 1)
                    %                 wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(film.waveConversion == imagePixel.wavelength(i))];  %pixel to update
                    wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(film.waveConversion == imagePixel.wavelength(i))];  %pixel to update
                    yPixel = film.resolution(1) + 1 - wantedPixel(2);
                    xPixel = wantedPixel(1);
                    
                    %check bounds - if out of bounds, do not display on film
                    if (xPixel >= 1 && xPixel <= film.resolution(1) && yPixel > 1 && yPixel <= film.resolution(2))
                        film.image(yPixel, xPixel, wantedPixel(3)) =  film.image(yPixel, xPixel, wantedPixel(3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
                    end
                    
                    %illustrations for debugging (out of bounds rays will
                    %still be displayed)
                    line(real([obj.origin(i, 3) intersectPosition(i, 3)]) ,  real([obj.origin(i, 2);  intersectPosition(i, 2)]));
                end
            end
        end
        
        %replicates the ray bundle to cover a series of wavelengths
        function obj = expandWavelengths(obj, wave)
            subLength = size(obj.origin, 1);
            obj.origin = repmat(obj.origin, [length(wave) 1]);
            obj.direction = repmat(obj.direction, [length(wave) 1]);
            tmp = (wave' * ones(1, subLength))'; %creates a vector representing wavelengths... for example: [400 400 400... 410 410 410... ..... 700]
            obj.wavelength = tmp(:);
        end
    end
end