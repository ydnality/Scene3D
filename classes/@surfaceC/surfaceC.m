classdef surfaceC <  handle
    % Create a lens surface object
    %
    %   lens = surfaceC(parameter, value, ....);
    %
    % Presently we only represent spherical surfaces and apertures.  Multi
    % element lenses (lensC) consist of a set of these surfaces and
    % apertures.
    %
    % This object is meant to be a class within the multi-element lens
    % objects, not a meant to be a stand-alone object.  This object
    % contains basic properties common to almost all lenses.
    %
    % get - beginning of a get function
    % calculateApertureSamples - outputs samples on the aperture to aid in
    %                            ray-tracing
    %
    % Parameter/vals:
    %    apertured, sradius, scenter, wave, zpos, n
    %
    % Examples:
    %    s = surfaceC;
    %    s = surfaceC('sRadius',2,'sCenter',[0 0 -2]);
    %
    % AL Vistasoft Copyright 2014
    
    properties
        
        % These are spherical surface properties
        name = 'default';
        type = 'surface';
        subtype = 'refractive';
        
        sRadius = 1;                % Sphere's radius
        sCenter = [0 0 0];          % Sphere's center position
        wave = 400:50:700;          % nm
        apertureD = 1;              % mm diameter
        n =  ones(7,1)              % index of refraction
        
    end
    
    methods (Access = public)
        
        % %%%%% Lens surface object constructor %%%%%%%
        function obj = surfaceC(varargin)
            %zpos must be assigned AFTER sCenter is assigned (after sCenter
            %in parameter declaration order).  Zpos assumes that lenses are
            %centered on z-axis.
            
            if isempty(varargin), return; end
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case {'apertured','aperturediameter'}
                        % Units are mm
                        obj.apertureD = varargin{ii+1};
                        
                    case 'sradius'
                        obj.sRadius = varargin{ii+1};
                        
                    case 'scenter'
                        obj.sCenter = varargin{ii+1};
                        
                    case 'wave'
                        obj.wave = varargin{ii+1};
                        
                    case {'zpos','zposition'}
                        %**MUST be assigned after sCenter is assigned
                        %assumes that lenses are centered on z axis
                        zPos = varargin{ii+1};
                        obj.sCenter = [ 0 0 obj.centerComputeFromZSRadius(zPos)];
                        
                    case 'n' % Index of refraction
                        % There should be one index of refraction for each
                        % wavelength
                        obj.n = varargin{ii+1};
                        if length(obj.n) ~= length(obj.wave)
                            warning('Index of refraction vector does not match wavelength');
                        end
                        
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
                case 'wave'
                    res = obj.wave;
                case 'zintercept'
                    %assumes a centered lens on the y = 0 axis;
                    res = obj.sCenter(3) - obj.sRadius;
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
        
         function set(obj,pName,val,varargin)
            pName = ieParamFormat(pName);
            switch pName
                case 'wave'
                    % The wavelength is annoying.
                    prevWave = obj.wave;
                    obj.wave = val;
                    obj.n = interp1(prevWave, obj.n, obj.wave, 'linear', 'extrap');
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
         end
    end
    
end