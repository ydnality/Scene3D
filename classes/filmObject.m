classdef filmObject <  handle
    % Create a film object
    %
    %   film = filmObject(position,size,wave,waveConversion,resolution);
    %
    % Spatial units throughout are mm
    %
    % For ray tracing, the sensor plane is called the 'film'.  At some
    % point we will need to be able to convert data from the ISET sensor
    % object and this film object.  In the fullness of time, they may be
    % closely coordinated.
    %
    % The film properties are
    %
    %   position - relative to lens
    %   size     - size in millimeters (height, width)
    %   wave     - sample wavelengths
    %   waveConversion - we will see
    %   resolution - Number of samples (pixels) in the film plane
    %
    % AL Vistasoft Copyright 2014
    
    properties
        position;
        size;        %in mm
        wave;
        waveConversion;
        resolution;  %decide whether this includes the 3rd dimension or not
        image;
    end
    
    methods
        
        %default constructor
        function obj = filmObject(position, size,  wave, waveConversion, resolution)
            
            if (ieNotDefined('position')), obj.position = [0 0 100];
            else                           obj.position = position;
            end
            
            if (ieNotDefined('size')),     obj.size = [48 48];
            else                           obj.size = size;
            end
            
            if (ieNotDefined('wave')),     obj.wave = [400 550 700];  % in nm;
            else                           obj.wave = wave;
            end
            
            %this field might go away soon
            if (ieNotDefined('waveConversion')),     obj.waveConversion = [400 1; 550 2; 700 3];  % in nm;
            else                           obj.waveConversion = waveConversion;
            end
            
            if (ieNotDefined('resolution')),     obj.resolution = [200 200 length(obj.wave)];
            else                           obj.size = resolution;
            end
            %TODO: error checking.  make sure all dimensions are good
            
            obj.image = zeros(obj.resolution);
        end
        
        
    end
    
end