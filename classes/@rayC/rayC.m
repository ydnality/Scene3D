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
        %wavelength;  %no longer used.  use obj.get('wavelength')
        waveIndex;   
        wave;  %this is the same wave in the film, which gives a range of wavelengths used
        drawSamples;   
        plotHandle = [];
    end
    
    methods (Access = private)
        function angles = toSphericalAngles(obj)
            %angles = toSphericalAngles(obj)
            %
            %Conversion of cartesian to spherical coordinates for
            %direction.
            %We will use official spherical coordinates as defined here: 
            %http://en.wikipedia.org/wiki/Spherical_coordinate_system
            %The first column will be the azimuth angle, and the second
            %column will be the polar angle.  All angles are in degrees.
            
            angles = zeros(size(obj.direction,1), 2);
            y = obj.direction(:, 2);
            x = obj.direction(:, 1);
            z = obj.direction(:, 3);
            azimuth = atan(y./x) * 180/pi;
            polar = acos(z./sqrt(x.^2 + y.^2 + z.^2));
            
            angles(:,1) = azimuth;
            angles(:,2) = polar;
        end
        
        function angles = toProjectedAngles(obj)
            %angles = toProjectedAngles(obj)
            
            %Conversion of cartesian to projected angles
            %
            %Let theta_x denote the angle cast by the vector, when
            %projected onto the z-x plane.
            %Similar, let theta_y denote the angle cast by the vector, when
            %projected onto the z-y plane
            %
            %The first column will be the theta_x angle, and the second
            %column will be the theta_y angle.  Angles are in degrees.
            
            angles = zeros(size(obj.direction,1), 2);
            y = obj.direction(:, 2);
            x = obj.direction(:, 1);
            z = obj.direction(:, 3);
            theta_x = atan(x./z) * 180/pi;
            theta_y = atan(y./z) * 180/pi;
            
            angles(:,1) = theta_x;
            angles(:,2) = theta_y;
        end       
        
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
                    case 'waveindex'
                        % The indexing that specifies the wavelength.  use
                        % obj.get('wavelength') to get a vector of
                        % wavelengths
                        obj.waveIndex = varargin{ii+1};
                    case 'wave'
                        % The wavelengths used for this ray object.
                        obj.wave = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
           end
        end
        
        function val = get(obj,param,varargin)
            % Get parameters about the rays
            p = ieParamFormat(param);
            switch p
                case 'nrays'
                    val = size(obj.origin,1);
                case 'sphericalangles'
                    val = obj.toSphericalAngles();
                case 'projectedangles'
                    val = obj.toProjectedAngles();
                case 'wave'
                    val = obj.wave;
                case 'waveindex'
                    val = obj.waveIndex;
                case 'wavelength'
                    val = zeros(size(obj.waveIndex));
                    val(isnan(obj.waveIndex)) = NaN;
                    liveInd = obj.get('liveIndices');
                    val(~isnan(obj.waveIndex)) = obj.wave(obj.waveIndex(liveInd));
                    val = val';
                case 'liveindices'  %return the indices of rays that are still alive
                    val = ~isnan(obj.waveIndex);
                otherwise
                    error('Unknown parameter %s\n',p);
            end
        end
        
        function set(obj,param,val)
        % set(param,val)
        %
        % sets various data members for the ray class
        %
        
            p = ieParamFormat(param);
            switch p
                case 'wave'
                    obj.wave = val;
                otherwise
                    error('Unknown parameter %s\n',p);
            end
        end
        
        function obj = projectOnPlane(obj, planeLocation)
            % ray.projectOnPlane(planeLocation)
            %
            % planeLocation: the z coordinate of a plane that is parallel to the
            %                x-y plane.
            % nLines: number of lines to draw for the illustration
            % intersects the rays with a plane at a specified location.
            % This is meant to make it easier to analyze the lightfield function.
            
            %remove dead rays
            
            liveIndices = ~isnan(obj.waveIndex);
            
            intersectZ = repmat(planeLocation, [size(obj.waveIndex(liveIndices, 1), 1) 1]);
            intersectT = (intersectZ - obj.origin(liveIndices, 3))./obj.direction(liveIndices, 3);
            intersectPosition = obj.origin(liveIndices, :) + obj.direction(liveIndices, :) .* repmat(intersectT, [1 3]);
            
            obj.origin(liveIndices, :) = intersectPosition;
        end
        
        function obj = recordOnFilm(obj, film, nLines)
        %records the ray on the film surface
        %nLines: number of lines to draw for the illustration
        
            if (ieNotDefined('nLines'))
                nLines = 0;
            end
            
            %parameters for plotting from lens to sensor
            lWidth = 0.1; lColor = [0 0.5 1]; lStyle = '-';
            
            
            %make a clone of the rays
            liveRays = rayC();
            liveRays.makeDeepCopy(obj);
            
            %calculate intersection point at sensor
            intersectZ = repmat(film.position(3), [size(liveRays.origin, 1) 1]);
            intersectT = (intersectZ - liveRays.origin(:, 3))./liveRays.direction(:, 3);
            liveRays.origin = liveRays.origin + liveRays.direction .* repmat(intersectT, [1 3]);
            
            %plot phase space and for the future - ray-raced lines
            if (nLines > 0)
                liveRays.plotPhaseSpace();

                 %draw debug lines
                % draw lines...TODO: figure this out...
                samps = obj.drawSamples;
                
                xCoordVector = [obj.origin(samps,3) liveRays.origin(samps,3) NaN([nLines 1])]';
                yCoordVector = [obj.origin(samps,2) liveRays.origin(samps,2) NaN([nLines 1])]';
                xCoordVector = real(xCoordVector(:));
                yCoordVector = real(yCoordVector(:));

                if isempty(obj.plotHandle), obj.plotHandle = vcNewGraphWin; end
                figure(obj.plotHandle);
                line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth,'LineStyle',lStyle);
                pause(0.1);
            end

            %remove dead rays
            deadIndices = isnan(obj.waveIndex);
            liveRays.origin(deadIndices, : ) = [];
            liveRays.direction(deadIndices, : ) = [];
            %liveRays.wavelength(deadIndices) = [];
            liveRays.waveIndex(deadIndices) = [];

            intersectPosition = liveRays.origin;
            
            % Record the real rays - if there are any
            if(~isempty(liveRays.origin))
                %imagePixel is the pixel that will gain a photon due to the traced ray
                imagePixel.position = [intersectPosition(:,2) intersectPosition(:, 1)];
                imagePixel.position = real(imagePixel.position); %add error handling for this
                imagePixel.position = round(imagePixel.position * film.resolution(2)/film.size(2) + ...
                    repmat(-film.position(2:-1:1)*film.resolution(2)/film.size(2)  + (film.resolution(2:-1:1) + 1)./2, [size(imagePixel.position,1) 1]));   %
                              
                imagePixel.wavelength = liveRays.get('wavelength');
                
                convertChannel = liveRays.waveIndex;
                %wantedPixel is the pixel that we wish to add 1 photon to
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
                
                %special case for length 1.  For some reason, hist has
                %issues with length 1.
                if (length(serialUniqueIndex(:)) == 1)
                    serializeFilm = film.image(:);
                    %when there is only 1 bin, it doesn't matter how many photons, so just add one
                    serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + 1; 
                    film.image = reshape(serializeFilm, size(film.image));
                elseif(length(serialUniqueIndex(:) > 0))
                    [countEntries] = hist(serialWantedPixel, serialUniqueIndex);
                    %serialize the film, then the indices, then add by countEntries
                    serializeFilm = film.image(:);
                    serializeFilm(serialUniqueIndex) = serializeFilm(serialUniqueIndex) + countEntries';  %this line might be problematic...
                    film.image = reshape(serializeFilm, size(film.image));
                else
                     warning('No photons were collected on film!');
                end
            end
        end
        
        function obj = expandWavelengths(obj, wave, waveIndex)
        % Replicates the ray bundle to cover a series of wavelengths
        %
        % The first ray trace step(from point source to the first lens element).
        % is performed without knowledge of wavelength - to save
        % computation and memory.  Once the ray enters the lens, then
        % wavelength information is necessary.  So this function replicates
        % every ray and assigns the wavelength information for each ray.  
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
%             tmp = (wave' * ones(1, subLength))'; 
%             obj.wavelength = tmp(:);
            
            obj.set('wave', wave);
            
            
            %assign the indices of this wavelength expansion 
            %  TODO: maybe make this simpler and cleaner somehow ...
            tmp = (waveIndex' * ones(1, subLength))';
            obj.waveIndex = tmp(:);
        end
        
        function removeDead(obj, deadIndices)
            %removeDead(deadIndices)
            
            %Removes dead rays (these are usually those that do not make it
            %out an aperture) by setting these dead indices to Nan.
            obj.origin(deadIndices, : ) = NaN;
            obj.direction(deadIndices, : ) = NaN;
            obj.waveIndex(deadIndices) = NaN;
        end
    end
end