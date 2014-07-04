classdef psfCameraObject <  handle
    % Create a film object
    %
    %   film = filmObject(position,size,wave,waveConversion,resolution);
    %
    % Spatial units throughout are mm
    %
    % For ray tracing, the sensor plane is called the 'film'.  At some
    % point we will need to convert data between the ISET sensor object and
    % this film object.  In the fullness of time, they may be closely
    % coordinated.
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
        lens;
        film;
        pointSource;
        rays;
    end
    
    methods
        
        %default constructor
        function obj = psfCameraObject(varargin)
            % psfCameraObject('lens',lens,'film',film,'point source',point);
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        obj.lens = varargin{ii+1};
                    case 'film'
                        obj.film = varargin{ii+1};
                    case 'pointsource'
                        obj.pointSource = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
        end
        
        function val = get(obj,param)
            % psfCamera.get('parameter name')
            % Start to set up the gets for this object
            val = [];
            param = ieParamFormat(param);
            switch param
                case 'spacing'
                    % Millimeters per sample
                    r = obj.film.resolution(1);
                    s = obj.film.size(1);
                    val = s/r;
                otherwise
                    error('unknown parameter %s\n',param)
            end
            
        end
        
        function oi = estimatePSF(obj,nLines, jitterFlag)
            %oi = estimatePSF(obj)
            %
            %estimates the PSF given the point source, lens, and film. 
            %returns the optical image of the film
            
            if ieNotDefined('nLines'),     nLines = false;     end
            if ieNotDefined('jitterFlag'), jitterFlag = false; end
            
            ppsfObjectFlag = true;
            obj.rays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfObjectFlag, jitterFlag);

            %duplicate the existing rays for each wavelength
            % Note that both lens and film have a wave, sigh.
            % obj.rays.expandWavelengths(obj.film.wave);
            obj.rays.expandWavelengths(obj.lens.wave);

            %lens intersection and raytrace
            obj.lens.rtThroughLens(obj.rays, nLines);

            %intersect with "film" and add to film
            obj.rays.recordOnFilm(obj.film);
            
            oi = obj.oiCreate();
            
        end
        
        function oi = oiCreate(obj)
        
            % Create an optical image from the camera (film) image data.
            oi = oiCreate;
            oi = initDefaultSpectrum(oi);
            oi = oiSet(oi,'wave', obj.film.wave);            
            oi = oiSet(oi,'photons',obj.film.image);
            
            % The photon numbers do not yet have meaning.  This is a hack,
            % that should get removed some day, to give the photon numbers
            % some reasonable level. 
            oi = oiAdjustIlluminance(oi,1);
            
            % Set focal length in meters
            oi = oiSet(oi,'optics focal length',obj.lens.focalLength/1000);
            
            % This isn't exactly the fnumber.  Do we have the aperture in
            % there?  Ask MP what to use for the multicomponent system.
            % This is just a hack to get something in there
            % Maybe this should be obj.lens.apertureMiddleD?
            fN = obj.lens.focalLength/obj.lens.surfaceArray(1).apertureD;
            oi = oiSet(oi,'optics fnumber',fN);
            
            % Estimate the horizontal field of view
            hfov = rad2deg(2*atan2(obj.film.size(1)/2,obj.lens.focalLength));
            oi = oiSet(oi,'hfov', hfov);

            % Set the name based on the distance of the sensor from the
            % final surface.  But maybe the obj has a name, and we should
            % use that?  Or the film has a name?
            temp = obj.film.position;
            filmDistance = temp(3);
            oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
        
        end
        
        function oi = showFilm(obj)
            % Create a new oi object, store it in the ISET database, and
            % show the oiWindow
            %
            oi = obj.oiCreate;
            vcAddAndSelectObject(oi); 
            oiWindow; 
        end
        
        function obj = recordOnFilm(obj)
             % Record the psf onto the film object in the camera.
             %
             % Uses the current ppsfRays and the film objects.  These are
             % both part of the camera object.
             % 
             % obj = recordOnFilm(obj)
             %
             
            obj.rays.recordOnFilm(obj.film); 
         end
    end
    
end