classdef psfCameraObject <  handle
    % Create a film object
    %
    %   film = filmObject(position,size,wave,waveConversion,resolution);
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
        lens;
        film;
        pointSource;
        rays;
    end
    
    methods
        
        %default constructor
        function obj = psfCameraObject(varargin)
            
%             ( lens, film, pointSource)
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        obj.lens = varargin{ii+1};
                    case 'film'
                        obj.film = varargin{ii+1};  %must be a 2 element vector
                    case 'pointsource'
                        obj.pointSource = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
             
            
%             if (ieNotDefined('lens')),     error('Input lens are required');
%             else                           obj.lens = lens;
%             end
%             
%             if (ieNotDefined('film')),     error('Input film are required');
%             else                           obj.film = film;
%             end
% 
%             if (ieNotDefined('pointSource')),     obj.pointSource = [0 0 -20000]; %default point source at -20m
%             else                           obj.pointSource = pointSource;  %units in mm
%             end
        end
        
        function oi = estimatePSF(obj)
            obj.rays = obj.lens.rayTraceSourceToLens(obj.pointSource);

            %duplicate the existing rays, and creates one for each
            %wavelength
            obj.rays.expandWavelengths(obj.film.wave);

            %lens intersection and raytrace
            obj.lens.rayTraceThroughLens(obj.rays);

            %intersect with "film" and add to film
            obj.rays.recordOnFilm(obj.film);
            
            %show the oi from the film
            oi = oiCreate;
            oi = initDefaultSpectrum(oi);
            oi = oiSet(oi, 'wave', obj.film.wave);
            oi = oiSet(oi,'photons',obj.film.image);
            optics = oiGet(oi,'optics');
            optics = opticsSet(optics,'focal length',obj.lens.focalLength/1000);
            optics = opticsSet(optics,'fnumber', obj.lens.focalLength/(2*1));
            oi = oiSet(oi,'optics',optics);
            hfov = rad2deg(2*atan2(obj.film.size(1)/2,obj.lens.focalLength));
            oi = oiSet(oi,'hfov', hfov);
        end
        
        
        function showFilm(obj)
            oi = oiCreate;
            oi = initDefaultSpectrum(oi);
            oi = oiSet(oi, 'wave', obj.film.wave);
            oi = oiSet(oi,'photons',obj.film.image);

            optics = oiGet(oi,'optics');
            optics = opticsSet(optics,'focal length',obj.lens.focalLength/1000);
            optics = opticsSet(optics,'fnumber', obj.lens.focalLength/(2*1));
            oi = oiSet(oi,'optics',optics);
            hfov = rad2deg(2*atan2(obj.film.size(1)/2,obj.lens.focalLength));
            oi = oiSet(oi,'hfov', hfov);

            temp = obj.film.position;
            filmDistance = temp(3);
            oi = oiSet(oi, 'name', ['filmDistance: ' num2str(filmDistance)]);
            vcAddAndSelectObject(oi); oiWindow; 
        end
        
        function obj = recordOnFilm(obj)
             %records the psf onto film of the current ppsfRays and the
             %film
             %obj = recordOnFilm(obj)
             %
             
            obj.rays.recordOnFilm(obj.film); 
         end
    end
    
end