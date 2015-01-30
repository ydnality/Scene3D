classdef filmC <  clonableHandleObject
    % Create a filmC
    %
    % Initiated by property/val pairs
    %
    %  film = filmC('position',val,'size',val,'wave',val,'resolution',val);
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
    % smallFilm = filmC('position', pos, 'size', sz, 'wave', wave, 'resolution', res);
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
        function obj = filmC(varargin)
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                
                val = varargin{ii+1};
                obj = obj.parameterAssign(p, val);
                
            end

            %TODO: error checking.  make sure all dimensions are good
            
            % Assign resolution to (nX,nY,nW) 
            % Not sure this is a good idea.  Leave wave off, I think. (BW)
            obj.resolution(3) = length(obj.wave);
            obj.image = zeros(obj.resolution);
        end
        
        function obj = parameterAssign(obj, p, val)
            %helper function for default constructor, assigns values (val)
            %to parameters (p) of class
            switch p
                case 'position'
                    obj.position = val;
                case 'size'
                    obj.size = val;  %must be a 2 element vector
                case 'wave'
                    obj.wave = val;
                case 'resolution'
                    obj.resolution = val;
                otherwise
                    error('Unknown parameter %s\n',val);
            end
        end
        
        function clear(obj)
            %clears the content of the film
            %
            %Example: film.clear();
            %
            
           obj.image = zeros(size(obj.image)); 
        end
        
    end
    
end