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
    %   lSurf = lensObject; lens.apertureSample = [51,51];
    %   fGrid = lSurf.fullGrid;
    %   aMask = lSurf.apertureMask; vcNewGraphWin; imagesc(aMask)
    %   aGrid = lSurf.apertureGrid; vcNewGraphWin; plot(aGrid.X(:),aGrid.Y(:),'o')
    %
    %   pointSource = [0 0 -50];
    %   rays = lSurf.rtSourceToEntrance(pointSource);
    %   ro =  rays.origin;
    %   rd =  rays.origin + rays.direction;
    %   vcNewGraphWin; hold on;
    %   for ii=1:10:size(rd,1), line([ro(ii,1),rd(ii,1)],[ro(ii,2),rd(ii,2)],[ro(ii,3),rd(ii,3)]); end
    %
    %  To get a variable, both actual slots and derived, use:
    %   lens.get('name')
    %
    %  To set a variable use.  (We may want to create a set() function, but
    %  we haven't yet).
    %   lens.name = 'test';
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

            for ii=1:2:length(varargin)
%                 p = ieParamFormat(varargin{ii});
                switch varargin{ii}
                    case 'apertureD'
                        % Units are mm
                        obj.apertureD = varargin{ii+1};
                        
%                     case 'centerPosition'
%                         obj.centerPosition = varargin{ii+1};
                        
%                     case 'diffractionEnabled'
%                         obj.diffractionEnabled = varargin{ii+1};
%                         
                    case 'sRadius'
                        obj.sRadius = varargin{ii+1};
                        
                    case 'sCenter'
                        obj.sCenter = varargin{ii+1};
                        
                    case 'wave'
                        obj.wave = varargin{ii+1};
                        
                    case 'n' % Index of refraction
                        obj.n = varargin{ii+1};
                        
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
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