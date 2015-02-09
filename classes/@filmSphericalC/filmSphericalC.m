classdef filmSphericalC < filmC
    %FILMSPHERICALC Summary of this class goes here
    %
    % The spherical sensor is defined as a set of sample points on the
    % surface of a sphere. The center of spherical sensor is assigned an
    % angle of (0,0). The pixels than sample an angular extent in latitude
    % and longitude. These extents can differ.
    %
    % The pixel positions are sampled uniformly in angles.
    % The positions in three space are calculated using simple geometry.
    %
    % AL. Vistasoft 2015.
    
    properties
        % Radius of curvature (same convention as lens). Units are mm.
        radius = -25;   
        
        % Note about the 'size' parameter for the spherical sensor
        %
        %
        % We don't have a row and column any more for the
        % size. Instead we use the length of the exposed perimeter of the
        % sphere as the size variable.  The surface area of the sensor can
        % be compute computed from the radius of curvature and the
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
                    %call super class parameter assigning function
                    parameterAssign@filmC(obj, p, val);  
            end
        end
        
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'angularsize'
                    %
                    % theta * r = perimeterSize
                    % theta     = perimeterSize/r
                    angleSize = obj.size./obj.radius;
                    res = angleSize;
                case 'steradiansize'
                    % The size (in steradians) of the sensor.  A sphere has
                    % 4pi steradians. The sensor surface covers a fraction
                    % of the sphere.  That fraction times 4pi is the
                    % steradian size of the sensor
                    disp('Steradian size not yet implemented.')
                case 'sphericalcenter'
                    %TODO: let's figure out sign convention...
                    res = [ 0 0 obj.position(3) + obj.radius]; 
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
    end
    
end

