classdef rayC <  clonableHandleObject
    % Create a ray object
    %
    % ray = rayC('origin', origin,'direction', direction, 'waveIndex', waveIndex, 'wave', wavelength)
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
            % Available parameters:
            %   nRays, sphericalAngles, projectedAngles, wave, waveIndex,
            %   wavelength, liveindices
            
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
                    if (mod(length(varargin), 2) ~= 0)
                        error('Incorrect parameter request. \n');
                    end
                    if (~isempty(varargin))
                        % this part deals with customized gets for specific
                        % wave indices and survived rays
                        val = obj.waveIndex;
                        for ii=1:2:length(varargin)
                            p = ieParamFormat(varargin{ii});
                            switch p
                                case 'survivedraysonly'
                                    survivedFlag = varargin{ii+1};
                                    if(survivedFlag)
                                       survivedRays = ~isnan(val); 
                                       val =  val(survivedRays);
                                    end
                                otherwise
                                    error('Unknown parameter %s\n',varargin{ii});
                            end
                        end
                    end
                    
                case 'wavelength'
                    val = zeros(size(obj.waveIndex));
                    val(isnan(obj.waveIndex)) = NaN;
                    liveInd = obj.get('liveIndices');
                    val(~isnan(obj.waveIndex)) = obj.wave(obj.waveIndex(liveInd));
                    val = val';
                case 'liveindices'  
                    % Rays with a waveIndex made it through the tracing
                    % path. We return the indices of rays that are still
                    % alive. We aren't sure why waveIndex is the right slot
                    % to check ... but it appears to be (BW).
                    val = ~isnan(obj.waveIndex);
                case 'liverays'
                    % Set the rays without a wavelength to empty  These
                    % remaining rays are the live rays.
                    val = rayC();
                    val.makeDeepCopy(obj);
                    
                    liveIndices = val.get('live indices');
                    val.origin(~liveIndices, : ) = [];
                    val.direction(~liveIndices, : ) = [];
                    val.waveIndex(~liveIndices) = [];
                    
                case 'origin'

                    %if no additional parameters are given, return raw
                    %origin matrix
                    val = obj.origin;
                    if (mod(length(varargin), 2) ~= 0)
                        error('Incorrect parameter request. \n');
                    end
                    if (~isempty(varargin) )
                        % this part deals with customized gets for specific
                        % wave indices and survived rays
                        
                        for ii=1:2:length(varargin)
                            p = ieParamFormat(varargin{ii});
                            switch p
                                case 'waveindex'
                                    wantedWaveIndex = varargin{ii+1};
                                    wantedWave = obj.get('waveIndex');
                                    if(~ieNotDefined('survivedFlag') && survivedFlag) %handles case if survivedrays called first
                                        wantedWave = wantedWave(survivedRays);
                                    end
                                    wantedWave = (wantedWave == wantedWaveIndex);
                                    val = val(:, wantedWave);
                                case 'survivedraysonly'
                                    survivedFlag = varargin{ii+1};
                                    if(survivedFlag)
                                       survivedRays = ~isnan(val(1,:)); %removes nans based off first coordinate
                                       val =  val(:, survivedRays);
                                    end
                                otherwise
                                    error('Unknown parameter %s\n',varargin{ii});
                            end
                        end
                    end
                    
                case 'direction'
                    
                    %if no additional parameters are given, return raw
                    %direction matrix
                    %consider putting this in a function so we don't need
                    %to define twice
                    val = obj.direction;
                    if (mod(length(varargin), 2) ~= 0)
                        error('Incorrect parameter request. \n');
                    end
                    if (~isempty(varargin))
                        % this part deals with customized gets for specific
                        % wave indices and survived rays
                        
                        for ii=1:2:length(varargin)
                            p = ieParamFormat(varargin{ii});

                            switch p
                                case 'waveindex'
                                    wantedWaveIndex = varargin{ii+1};
                                    wantedWave = obj.get('waveIndex');
                                    if(~ieNotDefined('survivedFlag') && survivedFlag) %handles case if survivedrays called first
                                        wantedWave = wantedWave(survivedRays);
                                    end
                                    wantedWave = (wantedWave == wantedWaveIndex);
                                    val = val(:, wantedWave);
                                case 'survivedraysonly'
                                    survivedFlag = varargin{ii+1};
                                    if(survivedFlag)
                                       survivedRays = ~isnan(val(1,:)); %removes nans based off first coordinate
                                       val =  val(:, survivedRays);
                                    end
                                otherwise
                                    error('Unknown parameter %s\n',varargin{ii});
                            end
                        end
                    end
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
        
        function obj = normalizeDir(obj)
           %obj = normalizeDir(obj)
           %normalizes all direction vectors so that the 2 norm is 1 
           obj.direction = normvec(obj.direction,'p',2,'dim',2);
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