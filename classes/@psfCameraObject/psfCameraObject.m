classdef psfCameraObject <  handle
    % Create a point spread camera object
    %
    % Spatial units throughout are mm
    %
    % AL Vistasoft Copyright 2014
    
    % Figure out the relationship between these rays and the ppsfRays in
    % the ppsfCameraObject.
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
                case 'imagecentroid'
                    % obj.get('image centroid')
                    % x,y positions (0,0) is center of the image centroid.
                    % Used for calculating centroid of the psf
                    % Could use obj.film.image for the data, rather than oi
                    % p_renderOiMatlabToolFull
                    % Figure out center pos by calculating the centroid of illuminance image
                    flm = obj.film;
                    img = flm.image;  img = sum(img,3);

                    % Force to unit area and flip up/down for a point spread
                    img = img./sum(img(:));
                    img = flipud(img);
                    % vcNewGraphWin; mesh(img);

                    % Calculate the weighted centroid/center-of-mass
                    xSample = linspace(-flm.size(1)/2, flm.size(1)/2, flm.resolution(1));
                    ySample = linspace(-flm.size(2)/2, flm.size(2)/2, flm.resolution(2));
                    [filmDistanceX, filmDistanceY] = meshgrid(xSample,ySample);
                    
                    % distanceMatrix = sqrt(filmDistanceX.^2 + filmDistanceY.^2);
                    val.X = sum(sum(img .* filmDistanceX));
                    val.Y = sum(sum(img .* filmDistanceY));
                    
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
            
            % Trace from the point source to the entrance aperture of the
            % multielement lens
            ppsfObjectFlag = true;
            obj.rays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfObjectFlag, jitterFlag);

            % Duplicate the existing rays for each wavelength
            % Note that both lens and film have a wave, sigh.
            % obj.rays.expandWavelengths(obj.film.wave);
            obj.rays.expandWavelengths(obj.lens.wave);

            %lens intersection and raytrace
            obj.lens.rtThroughLens(obj.rays, nLines);
            
            % Something like this?  Extend the rays to the film plane?
            % if nLines > 0; obj.rays.draw(obj.film); end

            %intersect with "film" and add to film
            obj.rays.recordOnFilm(obj.film);
            
            oi = obj.oiCreate();
            
        end
        
        function draw(obj,toFilm,nLines,apertureD)
            % psfCamera.draw(toFilm,nLines,apertureD)
            %
            % toFilm:    Logical that specifies whether to draw to film
            %            surface.  Default is false.
            % nLines:    Number of rays to use for the drawing
            %            Default is 200
            % apertureD: Sets the size of the film, when toFilm is true.
            %            Default is 100mm (really big)
            %
            % Show the ray trace lines to the film (sensor) plane
            
            if ieNotDefined('toFilm'), toFilm = false; end
            if ieNotDefined('nLines'), nLines = 200; end
            
            
            % Always true, for nicer picture, I guess.
            jitterFlag = true;
            
            
            % Not sure what to do here
            ppsfObjectFlag = false;

            % If toFilm is true, add the film surface as if it is an
            % aperture.  This will force the ray trace to continue to that
            % plane
            sArray = obj.lens.surfaceArray;  % Store the original
            
            wave      = obj.lens.wave;

            if toFilm
                disp('Drawing to film surface')
                
                % SHOULD BE Planar object.  But it won't draw to that
                sRadius   = 1e5;  % Many millimeters
                zPosition = obj.film.position(3);

                % We need a principled way to set this.
                if ieNotDefined('apertureD'), apertureD = 100; end
                
                obj.lens.surfaceArray(end+1) = ...
                    lensSurfaceObject('wave',wave,...
                    'aperture diameter',apertureD,...
                    'sRadius',sRadius,...
                    'zPosition',zPosition);
            end
            
            % Not sure why we need the ppsfObjectFlag (BW)
            obj.rays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfObjectFlag, jitterFlag);
            
            % Duplicate the existing rays for each wavelength
            % Note that both lens and film have a wave, sigh.
            % obj.rays.expandWavelengths(obj.film.wave);
            obj.rays.expandWavelengths(wave);

            %lens intersection and raytrace
            obj.lens.rtThroughLens(obj.rays, nLines);
            
            % Put it back the way you found it.
            if toFilm
                obj.lens.surfaceArray = sArray;
            end
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