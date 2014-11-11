classdef VoLTCameraC <  handle
    % Create a point spread camera object
    %
    % Spatial units throughout are mm
    %
    % AL Vistasoft Copyright 2014
    % See Also: ppsfCameraC.
    
    % Figure out the relationship between these rays and the ppsfRays in
    % the ppsfCameraC.
    properties
        scene;
        VoLTObject;
        lens;
        film;
        oi = [];  %derived data member
    end
    
    methods (Access = public)
        
        %default constructor
        function obj = VoLTCameraC(varargin)
            % psfCameraC('lens',lens,'film',film,'point source',point);
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        obj.lens = varargin{ii+1};
                    case 'voltobject'
                        obj.VoLTObject = varargin{ii+1};
                    case 'film'
                        obj.film = varargin{ii+1};
                   case 'scene'
                        obj.scene = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            
        end
        
    end
end