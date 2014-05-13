classdef ppsfCameraObject <  psfCameraObject
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
%         lens;
%         film;
%         pointSource;
        ppsfRays;
    end
    
    methods
         function obj = ppsfCameraObject( lens, film, pointSource)
             %deals with omitted inputs - super class decides the default
             %values
             if (ieNotDefined('lens')), lens = []; end
             if (ieNotDefined('film')), film = []; end
             if (ieNotDefined('pointSource')), pointSource = []; end
             
             obj = obj@psfCameraObject(lens,film, pointSource);
         end
         
         function ppsfReturn = estimatePPSF(obj)
            %calculate the origin and direction of the rays
            %     rays.traceSourceToLens(pointSources(curInd, :), lens);

            disp('-----trace source to lens-----');
            tic
            obj.ppsfRays = obj.lens.rayTraceSourceToLens(obj.pointSource(1, :));
            apertureSamples = obj.lens.apertureSample;
            obj.ppsfRays = ppsfObject(obj.ppsfRays.origin, obj.ppsfRays.direction, obj.ppsfRays.wavelength, 0, 0, apertureSamples);  %think of a best way to put in aperture sample location
            toc

            %duplicate the existing rays, and creates one for each
            %wavelength
            disp('-----expand wavelenghts-----');
            tic
            obj.ppsfRays.expandWavelengths(obj.film.wave);
            toc

            %lens intersection and raytrace
            disp('-----rays trace through lens-----');
            tic
            obj.lens.rayTraceThroughLens(obj.ppsfRays);
            toc
            
            ppsfReturn = obj.ppsfRays;
         end
    end
    
end