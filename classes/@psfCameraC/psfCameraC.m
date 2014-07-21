classdef psfCameraC <  handle
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
    
    methods (Access = public)
        
        %default constructor
        function obj = psfCameraC(varargin)
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
        
    end
    
end