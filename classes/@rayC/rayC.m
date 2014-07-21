classdef rayC <  clonableHandleObject
    % Create a ray object
    %
    % ray = rayObject(origin,direction,wavelength)
    %
    % Example:
    %   r = rayC
    %
    % Programming:
    %   A lot of these seem to overlap with other routines.  Delete one of
    %   them, either here or the lens version or ... (BW)
    %
    % AL Vistasoft Copyright 2014
    
    properties
        origin;
        direction;
        wavelength;
        waveIndex;
    end
    
    methods (Access = public)
          
        function obj = rayC(varargin)
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
        
        function obj = projectOnPlane(obj, planeLocation)
        % ray.projectOnPlane(planeLocation)
        %
        % planeLocation: the z coordinate of a plane that is parallel to the
        %                x-y plane.
        % intersects the rays with a plane at a specified location.  
        % This is meant to make it easier to analyze the lightfield function.
        
        %remove dead rays
        liveIndices = ~isnan(obj.wavelength);
        
        intersectZ = repmat(planeLocation, [size(obj.wavelength(liveIndices, 1), 1) 1]);
        intersectT = (intersectZ - obj.origin(liveIndices, 3))./obj.direction(liveIndices, 3);
        intersectPosition = obj.origin(liveIndices, :) + obj.direction(liveIndices, :) .* repmat(intersectT, [1 3]);
        
        obj.origin(liveIndices, :) = intersectPosition;
        
        end
        
        function obj = recordOnFilm(obj, film)
        %records the ray on the film surface
        
            %make a clone of the rays
            liveRays = rayC();
            liveRays.makeDeepCopy(obj);
            
            %remove dead rays
            deadIndices = isnan(obj.waveIndex);      
            
            liveRays.origin(deadIndices, : ) = [];
            liveRays.direction(deadIndices, : ) = [];
            liveRays.wavelength(deadIndices) = [];
            liveRays.waveIndex(deadIndices) = [];
            
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
                              
                imagePixel.wavelength = liveRays.wavelength;
                
                convertChannel = liveRays.waveIndex;
                wantedPixel = [imagePixel.position(:, 1) imagePixel.position(:,2) convertChannel];  %pixel to update
                
                recordablePixels =and(and(and(wantedPixel(:, 1) >= 1,  wantedPixel(:,1) <= film.resolution(1)), (wantedPixel(:, 2) > 1)), wantedPixel(:, 2) <= film.resolution(2));
                
                %remove the nonrecordablePixels
                wantedPixel = wantedPixel(recordablePixels, :);
                
                %correct for y coordinates
                wantedPixel(:, 1) =  film.resolution(1) + 1 - wantedPixel( :, 1);
                
                %make a histogram of wantedPixel in anticipation of adding
                %to film
                %  [count bins] = hist(single(wantedPixel));
                %  [count bins] = hist(single(wantedPixel), unique(single(wantedPixel), 'rows')); 
                uniqueEntries =  unique(single(wantedPixel), 'rows');
                
                % Serializes the unique entries
                % I had an error here once where the wavelength in the
                % uniqueEntries was 7 but there were only 4 wavelength
                % dimensions in the film image.  Hmm. (BW).
                serialUniqueIndex = sub2ind(size(film.image), uniqueEntries(:,1), uniqueEntries(:,2), uniqueEntries(:,3));
                
                serialUniqueIndex = sort(serialUniqueIndex);
                
                serialWantedPixel = sub2ind(size(film.image), single(wantedPixel(:,1)), single(wantedPixel(:,2)), single(wantedPixel(:,3)));
                
                
                [countEntries] = hist(serialWantedPixel, serialUniqueIndex);

                %serialize the film, then the indices, then add by countEntries
                serializeFilm = film.image(:);
                serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + countEntries';
                
                film.image = reshape(serializeFilm, size(film.image));
                
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