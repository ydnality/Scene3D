classdef lensSurfaceObject <  handle
    % Create a lens object
    %
    %   lens = lensSurfaceObject( parameter, value, ....);
    %
    % Presently we only represent spherical lenses and apertures.
    %
    % This is meant to be the superclass for additional lens objects.  This
    % is not meant to be a stand-alone object.  This object contains some
    % basic properties common to almost all lenses.  Also - it contains the
    % calculateApertureSamples function which outputs samples on the
    % aperture to aid in ray-tracing.
    %
    % TODO: determine which aperture is it best to sample?
    %
    % Is it possible to write this in the form
    %  lensSurfaceObject(lensType, relevant parameter list ... )
    %
    % Examples:
    %
    % AL Vistasoft Copyright 2014
    
    properties
        
        % These are surfaces of spherical lenses for now
        name = 'default';
        type = 'surface';
        
        sRadius = 1;      % Sphere's radius
        sCenter = [0 0 0];      % Sphere's center position
        wave = 400:50:700;         % nm
        apertureD = 1;        % mm diameter
%         centerPosition = [0 0 0];  % almost always on axis
        n = ones(7,1);   %index of refraction   %check dimension order
    end
    
    methods
        
        % %%%%% Lens surface object constructor %%%%%%%
        function obj = lensSurfaceObject(varargin)
            %zpos must be assigned AFTER sCenter is assigned (after sCenter
            %in parameter declaration order).  Zpos assumes that lenses are
            %centered on z-axis.
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'apertured'
                        % Units are mm
                        obj.apertureD = varargin{ii+1};
                        
%                     case 'centerPosition'
%                         obj.centerPosition = varargin{ii+1};
                        
%                     case 'diffractionEnabled'
%                         obj.diffractionEnabled = varargin{ii+1};
%                         
                    case 'sradius'
                        obj.sRadius = varargin{ii+1};
                        
                    case 'scenter'
                        obj.sCenter = varargin{ii+1};
                        
                    case 'wave'
                        obj.wave = varargin{ii+1};
                        
                    case 'zpos' %**MUST be assigned after sCenter is assigned
                        %assumes that lenses are centered on z axis
                        zPos = varargin{ii+1};
                        obj.sCenter = [ 0 0 obj.centerComputeFromZSRadius(zPos)];
                        
                    case 'n' % Index of refraction
                        obj.n = varargin{ii+1};
                        
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
        end
        
        
        function center = centerComputeFromZSRadius(obj, zPos)
            %computes the spherical center given the z position and the
            %spherical radius.
            
            center = zPos + obj.sRadius;
        end
        
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'name'
                    res = obj.name;
                case 'type'
                    res = obj.type;
                case 'zintercept'
                    %assumes a centered lens on the y = 0 axis;
                    res = obj.sCenter(3) - obj.sRadius;
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end

    end
    
end