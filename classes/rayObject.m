classdef rayObject <  clonableHandleObject
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
        waveIndex;
    end
    
    methods
          
        function obj = rayObject(varargin)
            % Constructor

           for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'origin'
                        % Define, please
                        obj.origin = varargin{ii+1};
                    case 'direction'
                        % Must be a 2 element vector that represents ???
                        obj.direction = varargin{ii+1};  
                    case 'wavelength'
                        % In nanometers.  Can it be a vector? Or just a
                        % scalar?
                        obj.wavelength = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
           end
        end
        
        % Get parameters about the rays
        function val = get(obj,param,varargin)
            p = ieParamFormat(param);
            switch p
                case 'nrays'
                    val = size(obj.origin,1);
                otherwise
                    error('Unknown parameter %s\n',p);
            end
        end
        
        
        function obj = traceSourceToLens(obj, curPointSource, lens)
            % Deprecate?
            %traces rays from a point source to a sampling function on the lens
            % Is this old code that was moved to the lens object?
            obj.origin = repmat(curPointSource, [size(lens.apertureSample.Y(:), 1) 1] );   %the new origin will just be the position of the current light source
            obj.direction = [(lens.apertureSample.X(:) -  obj.origin(:,1)) (lens.apertureSample.Y(:) -  obj.origin(:,2)) (lens.centerPosition(3) - obj.origin (:,3)) .* ones(size(lens.apertureSample.Y(:)))];
            obj.direction = obj.direction./repmat( sqrt(obj.direction(:, 1).^2 + obj.direction(:, 2).^2 + obj.direction(:,3).^2), [1 3]); %normalize direction
        end
        

        function obj =  traceThroughLens(obj, lens)
            % Deprecate?
            % This seems like a duplicate of the lens function.
            % If so, maybe we should get rid of it.  Or that one.
            %
            % Performs ray-trace of the lens, given an input bundle or rays
            % outputs the rays that have been refracted by the lens
            % TODO: consdier moving this to the lens - II think it was.
            %       So, delete? (BW)
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
                    %                     hold on;
                    %                     line([ray.origin(3) intersectPosition(3) ], [ray.origin(2) intersectPosition(2)] ,'Color','b','LineWidth',1);
                    %
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
        
        function obj = projectOnPlane(obj, planeLocation)
        % ray.projectOnPlane(planeLocation)
        %
        % planeLocation: the z coordinate of a plane that is parallel to the
        %                x-y plane.
        % intersects the rays with a plane at a specified location.  
        % This is meant to make it easier to analyze the lightfield function.
        
        %remove dead rays
        liveIndices = ~isnan(obj.wavelength);
        
        %             planeCoordinate = [ 0 0 planeLocation];
        intersectZ = repmat(planeLocation, [size(obj.wavelength(liveIndices, 1), 1) 1]);
        intersectT = (intersectZ - obj.origin(liveIndices, 3))./obj.direction(liveIndices, 3);
        intersectPosition = obj.origin(liveIndices, :) + obj.direction(liveIndices, :) .* repmat(intersectT, [1 3]);
        
        obj.origin(liveIndices, :) = intersectPosition;
        
        end
        
        function obj = recordOnFilm(obj, film)
        %records the ray on the film surface
        
            %make a clone of the rays
            liveRays = rayObject();
            liveRays.makeDeepCopy(obj);
            
            %remove dead rays
            deadIndices = isnan(obj.waveIndex);      
            
%             props = properties(liveRays);
%             for i = 1:length(props)
%                 liveRays.(props{i})(deadIndices) = [];
                



                %***TODO - perhaps remove all NANs..... for all data
                %members
                liveRays.origin(deadIndices, : ) = [];
                liveRays.direction(deadIndices, : ) = [];
                liveRays.wavelength(deadIndices) = [];
                liveRays.waveIndex(deadIndices) = [];                
                
                %                 liveRays.apertureSamples.X(deadIndices) = [];   %THIS
                %                 NEEDS TO BE FIXED!!! the difference between rayobject
                %                 andd ppsfObject
                %                 liveRays.apertureSamples.Y(deadIndices) = [];
                %                 liveRays.apertureLocation(deadIndices, :) = [];
                %             end
                
            
            % Record the real rays - if there are any
            if(~isempty(liveRays.origin))
                
                %calculate intersection point at sensor
                intersectZ = repmat(film.position(3), [size(liveRays.origin, 1) 1]);
                intersectT = (intersectZ - liveRays.origin(:, 3))./liveRays.direction(:, 3);
                intersectPosition = liveRays.origin + liveRays.direction .* repmat(intersectT, [1 3]);
                
                %imagePixel is the pixel that will gain a photon due to the traced ray
                imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
                imagePixel.position = real(imagePixel.position); %add error handling for this
                imagePixel.position = round(imagePixel.position * film.resolution(2)/film.size(2) + ...
                    repmat(-film.position(2:-1:1)*film.resolution(2)/film.size(2)  + (film.resolution(2:-1:1) + 1)./2, [size(imagePixel.position,1) 1]));   %
               
                %scale the position to a sensor position
                %             imagePixel.position(imagePixel.position < 1) = 1; %make sure pixel is in range
                %             imagePixel.position = real(min(imagePixel.position, repmat(film.resolution(1:2), [size(imagePixel.position,1) 1])));
                
                imagePixel.wavelength = liveRays.wavelength;
                
                
                %attempted vectorized version - runs out of memory
                %
                %                 %todo - somehow fix this - it is not general enough
                %                 convertChannel = uint8((imagePixel.wavelength - 400)/10 + 1);  %THIS LINE IS PROBLEMATIC!!!!!
                
                
                convertChannel = liveRays.waveIndex;
                wantedPixel = [imagePixel.position(:, 1) imagePixel.position(:,2) convertChannel];  %pixel to update
                
                recordablePixels =and(and(and(wantedPixel(:, 1) >= 1,  wantedPixel(:,1) <= film.resolution(1)), (wantedPixel(:, 2) > 1)), wantedPixel(:, 2) <= film.resolution(2));
                
                %remove the nonrecordablePixels
                wantedPixel = wantedPixel(recordablePixels, :);
                
                %correct for y coordinates
                wantedPixel(:, 1) =  film.resolution(1) + 1 - wantedPixel( :, 1);
                
                %make a histogram of wantedPixel in anticipation of adding
                %to film
                %                 [count bins] = hist(single(wantedPixel));
                %                  [count bins] = hist(single(wantedPixel), unique(single(wantedPixel), 'rows'));
                
                uniqueEntries =  unique(single(wantedPixel), 'rows');
                
                % Serializes the unique entries
                % I had an error here once where the wavelength in the
                % uniqueEntries was 7 but there were only 4 wavelength
                % dimensions in the film image.  Hmm. (BW).
                serialUniqueIndex = sub2ind(size(film.image), uniqueEntries(:,1), uniqueEntries(:,2), uniqueEntries(:,3));
                
                serialUniqueIndex = sort(serialUniqueIndex);
                
                serialWantedPixel = sub2ind(size(film.image), single(wantedPixel(:,1)), single(wantedPixel(:,2)), single(wantedPixel(:,3)));
                
                
                [countEntries] = hist(serialWantedPixel, serialUniqueIndex);

                %old slower code
                %                 countEntries = zeros(length(uniqueEntries), 1);
                %count the entries - see if there's a better way to do
                %this...
                %                 for ii = 1:length(uniqueEntries)
                %                     countEntries(ii) = sum(sum(wantedPixel == repmat(uniqueEntries(ii,:), [length(wantedPixel) 1]), 2) == 3);
                %                 end
                
                
                %serialize the film, then the indices, then add by
                %countEntries
                serializeFilm = film.image(:);
                serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + countEntries';
                
                film.image = reshape(serializeFilm, size(film.image));
                
                %                 film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) =  film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
                %this line doesn't work - wrong indexing, and cannot count
                %mroe than 1
                %                 film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) =  film.image(wantedPixel(:, 1), wantedPixel(:, 2), wantedPixel(:, 3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
                
                %                 for i = 1:10:size(obj.origin , 1)
                %                     %illustrations for debugging (out of bounds rays will
                %                     %still be displayed)
                %                     line(real([obj.origin(i, 3) intersectPosition(i, 3)]) ,  real([obj.origin(i, 2);  intersectPosition(i, 2)]));
                %                 end
                
                
                %                 %non-vectorized
                %                 %add a value to the intersection position
                %                 for i = 1:size(obj.origin , 1)
                %                     %                 wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(film.waveConversion == imagePixel.wavelength(i))];  %pixel to update
                %                     wantedPixel = [imagePixel.position(i,1) imagePixel.position(i,2) find(film.waveConversion == imagePixel.wavelength(i))];  %pixel to update
                %                     yPixel = film.resolution(1)+1 - wantedPixel(1);
                %                     xPixel = wantedPixel(2);
                %
                %                     %check bounds - if out of bounds, do not display on film
                %                     if (xPixel >= 1 && xPixel <= film.resolution(1) && yPixel > 1 && yPixel <= film.resolution(2))
                %                         film.image(yPixel, xPixel, wantedPixel(3)) =  film.image(yPixel, xPixel, wantedPixel(3)) + 1;  %sensor.image(imagePixel(:,1), imagePixel(:,2)) + 1;
                %                     end
                %
                %                     %illustrations for debugging (out of bounds rays will
                %                     %still be displayed)
                %                     line(real([obj.origin(i, 3) intersectPosition(i, 3)]) ,  real([obj.origin(i, 2);  intersectPosition(i, 2)]));
                %                 end
            end
        end
        
        % Replicates the ray bundle to cover a series of wavelengths
        % Why do we need to do this?  (BW)
        function obj = expandWavelengths(obj, wave, waveIndex)
            
            if ieNotDefined('waveIndex')
                waveIndex = 1:length(wave);
            end
            
            % This is the number of rays.
            subLength = size(obj.origin, 1);
            
            % Not sure why we need to repmat the origin and direction
            obj.origin = repmat(obj.origin, [length(wave) 1]);
            obj.direction = repmat(obj.direction, [length(wave) 1]);
            
            % Creates a vector representing wavelengths... for example:
            % [400 400 400... 410 410 410... ..... 700] 
            tmp = (wave' * ones(1, subLength))'; 
            obj.wavelength = tmp(:);
            
            %assign the indices of this wavelength expansion 
            %  TODO: maybe make this simpler and cleaner somehow ...
            tmp = (waveIndex' * ones(1, subLength))';
            obj.waveIndex = tmp(:);
        end
    end
end