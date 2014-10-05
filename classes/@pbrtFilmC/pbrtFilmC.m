classdef pbrtFilmC <  clonableHandleObject
    % Create a pbrtFilmC
    %
    % Initiated by property/val pairs
    %
    %  film = pbrtFilmC('position',val,'size',val,'wave',val,'resolution',val);
    %
    % Spatial units throughout are mm
    %
    % The film properties are
    %
    %   position   - film position in lens coordinate frame
    %   size       - film size in millimeters (height, width)
    %   wave       - wavelengths
    %   resolution - Number of spatial samples (pixels) in the film plane
    %
    % Example:
    % wave = 500; sz = [10,10]; pos = [0 0 20]; res = [150 150 1];
    % smallFilm = pbrtFilmC('position', pos, 'size', sz, 'wave', wave, 'resolution', res);
    %
    % AL Vistasoft Copyright 2014
    
    % This object has no methods yet.  It gets attached to other objects that need this
    % information, like psfCameraC.
    
    properties
        position = [ 0 0 100];
        size = [48 48];            % in mm (Too big).
        wave = 400:50:700;
        resolution = [200 200];  
        image;
    end
    
    methods
        
        %default constructor
        function obj = pbrtFilmC(varargin)
            
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

            %TODO: error checking.  make sure all dimensions are good
            
            % Assign resolution to (nX,nY,nW) 
            % Not sure this is a good idea.  Leave wave off, I think. (BW)
            obj.resolution(3) = length(obj.wave);
            obj.image = zeros(obj.resolution);
        end
        
        
    end
    
end