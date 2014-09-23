classdef LFC
    %LFC (light field class) consists of a light field and operations that
    %can be performed on a light field.
    %   
    % See Also: ppsfRayC
    
    properties
        %A 4xn matrix.  The first two rows signify the X and Y positions of
        %rays.  The last two rows signifies the x and y coordinates of a
        %unit direction vector for the rays.
        LF;  
        %The corresponding waveIndex for each ray.  nx1 in size.  n must
        %match for LF and waveIndex.
        waveIndex;
        %The wavelengths referenced by waveIndex.  wave(waveIndex) will
        %give all the wavelengths of each ray of the light field.
        wave;
    end
    
    methods
        function obj = LFC(varargin)
        % Initialization of the Light Field Class Object
        
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lf'  %TODO: error checking for LF matching waveIndex
                        obj.LF = varargin{ii+1};
                    case 'waveindex'
                        obj.waveIndex = varargin{ii+1};
                    case 'wave'
                        obj.wave = varargin{ii+1};
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
        end     
        
                % Get VoLT properties
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'lf'
                    res = obj.LF;
                case 'waveindex'
                    res = obj.waveIndex;
                case 'wave'
                    res = obj.wave;
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
    end
    
end

