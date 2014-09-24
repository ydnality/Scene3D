classdef ppsfCameraC <  psfCameraC
    % The superclass is psfCameraC.
    %
    % This ppsfRay object is added to this subclass of psfCameraC.
    %
    %   ppsfCamera = ppsfCameraC;
    %
    % Spatial units throughout are mm
    %
    % The camera specifies a lens, film, and input point.  
    %
    % Examples:
    %    camera = ppsfCameraC;
    %    camera = ppsfCameraC('lens',lens,'film',film,'point source',ps);
    %
    % AL Vistasoft Copyright 2014
    %
    % See also: psCreate, lensC, pbrtFilmC (which may become filmC).
    
    properties
        ppsfRays;
    end
    
    methods (Access = public)
         function obj = ppsfCameraC(varargin)
             % Initialize so ppsf = ppsfCameraC; will work
             lens = [];
             film = [];
             pointSource = [];
             
             % Set the parameters
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
             
             %Use the psfCameraC constructor
             obj = obj@psfCameraC('lens', lens,...
                 'film', film, ...
                 'pointSource', pointSource);
             
             % Should we clear out rays on return and only keep ppsfRays?
             % What is the plan here?
             
         end
         
         function obj = recordOnFilm(obj, nLines)
             %records the psf onto film 
             % The ppsfRays 
             % film
             % nLines: >0 for drawing existing debug lines.  0 or false for
             % not drawing them
             %obj = recordOnFilm(obj)
             %
             
            if ieNotDefined('nLines'), nLines = false; end
             
            obj.ppsfRays.recordOnFilm(obj.film, nLines); 
         end
    end
    
end