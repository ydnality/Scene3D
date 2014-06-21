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
         function obj = ppsfCameraObject(varargin)
             
             
%              ( lens, film, pointSource)
             %deals with omitted inputs - super class decides the default
             %values
             
             for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        lens = varargin{ii+1};
                    case 'film'
                        film = varargin{ii+1};  %must be a 2 element vector
                    case 'pointsource'
                        pointSource = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
             end
             
             %              if (ieNotDefined('lens')), lens = []; end
             %              if (ieNotDefined('film')), film = []; end
             %              if (ieNotDefined('pointSource')), pointSource = []; end
             obj = obj@psfCameraObject('lens', lens,'film', film, 'pointSource', pointSource);  %this lets the psfCameraObject do error handling
         end
         
         function ppsfReturn = estimatePPSF(obj,nLines)
            % calculate the origin and direction of the exiting rays
            %
            % nLines is the number of lines to draw on the diagram.
            % For no diagram set nLines to 0 (false).  This is the default.
            % ppsfCamera.estimatePPSF(nLines)    
            
            if ieNotDefined('nLines'), nLines = false; end
            
            disp('-----trace source to lens-----');
            tic
            ppsfObjectFlag = true;
            obj.ppsfRays = obj.lens.rtSourceToEntrance(obj.pointSource, ppsfObjectFlag);
            toc
            
            %duplicate the existing rays, and creates one for each
            %wavelength
            disp('-----expand wavelenghts-----');
            tic
            obj.ppsfRays.expandWavelengths(obj.lens.wave);
            toc

            %lens intersection and raytrace
            disp('-----rays trace through lens-----');
            tic
            obj.lens.rtThroughLens(obj.ppsfRays,nLines);
            toc
            
            %project rays onto the z = 0 plane for a proper light field
            obj.ppsfRays.projectOnPlane(0);
            obj.ppsfRays.pointSourceLocation = obj.pointSource;
            ppsfReturn = obj.ppsfRays;
            
         end
         
         
         function obj = recordOnFilm(obj)
             %records the psf onto film of the current ppsfRays and the
             %film
             %obj = recordOnFilm(obj)
             %
             
            obj.ppsfRays.recordOnFilm(obj.film); 
         end
    end
    
end