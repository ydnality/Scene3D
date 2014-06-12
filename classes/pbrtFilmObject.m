classdef pbrtFilmObject <  handle
    % Create a pbrtFilmObject
    %
    %   film = pbrtFilmObject(position,size,wave,waveConversion,resolution);
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
        position = [ 0 0 100];
        size = [48 48];        %in mm
        wave = 400:50:700;
%         waveConversion;
        resolution = [200 200 7];  %decide whether this includes the 3rd dimension or not
        image;
    end
    
    methods
        

        %default constructor
        function obj = pbrtFilmObject(varargin)
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'position'
                        obj.position = varargin{ii+1};
                    case 'size'
                        obj.size = varargin{ii+1};  %must be a 2 element vector
                    case 'wave'
                        obj.wave = varargin{ii+1};
                    case 'resolution'
                        obj.resolution = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
%             (position, size,  wave, waveConversion, resolution)
            
%             if (ieNotDefined('position')), obj.position = [0 0 100];
%             else                           obj.position = position;
%             end
%             
%             if (ieNotDefined('size')),     obj.size = [48 48];
%             else                           obj.size = size;
%             end
%             
%             if (ieNotDefined('wave')),     obj.wave = [400 550 700];  % in nm;
%             else                           obj.wave = wave;
%             end
%             
%             %this field might go away soon
%             if (ieNotDefined('waveConversion')),     obj.waveConversion = [400 1; 550 2; 700 3];  % in nm;
%             else                           obj.waveConversion = waveConversion;
%             end
%             
%             if (ieNotDefined('resolution')),     obj.resolution = [200 200 length(obj.wave)];
%             else                           obj.resolution = resolution;
%             end
            %TODO: error checking.  make sure all dimensions are good
            
            %reassign resolution and image size after wave is set
            obj.resolution(3) = length(obj.wave);
            obj.image = zeros(obj.resolution);
        end
        
        
    end
    
end