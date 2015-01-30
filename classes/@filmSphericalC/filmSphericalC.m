classdef filmSphericalC < filmC
    %FILMSPHERICALC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius = -25;
        % how will the film size translate to spherical film size?? let's
        % make that the perimeter length of the sensor
    end
    
    methods
        %default constructor
        function obj = filmSphericalC(varargin)
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                
                val = varargin{ii+1};
                obj = obj.parameterAssign(p, val);
            end
        end
        
        function obj = parameterAssign(obj, p, val)
            %helper function for default constructor, assigns values (val)
            %to parameters (p) of class
            switch p
                case 'radius'
                    obj.radius = val;
                otherwise
                    parameterAssign@filmC(obj, p, val);  %call super class parameter assigning function
            end
        end
        
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'angularsize'
                    
                    % theta *r = perimeterSize
                    % theta = perimeterSize/r
                    angleSize = obj.size./obj.radius;
                    res = angleSize;
                case 'sphericalcenter'
                    res = [ 0 0 obj.position(3) + obj.radius]; %TODO: let's figure out sign convention...
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
    end
    
end

