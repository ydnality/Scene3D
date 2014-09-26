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
        %
        %  LFC('LF',lf,'wave',wave,'waveIdx',waveIdx);
        %
        
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
                    %if no additional parameters are given, return raw
                    %lf matrix
                    
                    res = obj.LF;
                    if (mod(length(varargin), 2) ~= 0)
                        error('Incorrect parameter request. \n');
                    end
                    if (~isempty(varargin))
                        % this part deals with customized gets for specific
                        % wave indices and survived rays
                        for ii=1:2:length(varargin)
                            p = ieParamFormat(varargin{ii});
                            switch p
                                case 'waveindex'
                                    wantedWaveIndex = varargin{ii+1};
                                    wantedWave = obj.get('waveIndex');
                                    if(~ieNotDefined('survivedFlag') && survivedFlag) %handles case if survivedrays called first
                                        wantedWave = wantedWave(survivedRays);
                                    end
                                    wantedWave = (wantedWave == wantedWaveIndex);
                                    res = res(:, wantedWave);
                                case 'survivedraysonly'
                                    survivedFlag = varargin{ii+1};
                                    if(survivedFlag)
                                       survivedRays = ~isnan(res(1,:)); %removes nans based off first coordinate
                                       res =  res(:, survivedRays);
                                    end
                                otherwise
                                    error('Unknown parameter %s\n',varargin{ii});
                            end
                        end
                    end
                case 'survivedrays'
                    res = isnan(obj.waveIndex);
                case 'waveindex'
                    %if no additional parameters are given, return raw
                    %waveindex
                    
                    res = obj.waveIndex;
                    %if additional parameters are given, "filter out"
                    %specific wavindex entries.
                    if (mod(length(varargin), 2) ~= 0)
                        error('Incorrect parameter request. \n');
                    end
                    if (~isempty(varargin))
                        % this part deals with customized gets for specific
                        % wave indices and survived rays
                        for ii=1:2:length(varargin)
                            p = ieParamFormat(varargin{ii});
                            switch p
                                case 'survivedraysonly'
                                    survivedFlag = varargin{ii+1};
                                    if(survivedFlag)
                                       survivedRays = ~isnan(res(:)); %removes nans based off first coordinate
                                       res =  res(survivedRays);
                                    end
                                otherwise
                                    error('Unknown parameter %s\n',varargin{ii});
                            end
                        end
                    end

                case 'wave'
                    res = obj.wave;
                case 'ray'
                    % Convert the light field to a ray representation with
                    % origin and direction
                    % This should probably become a rayC object that knows
                    % about its wavelength.
                    rayOrigin = zeros(3, size(obj.LF, 2));
                    rayDir = rayOrigin;
                    
                    rayOrigin(1,:) = obj.LF(1,:);
                    rayOrigin(2,:) = obj.LF(2,:);
                    rayOrigin(3,:) = 0;        % Probably should be a variable name
                    
                    rayDir(1,:) = obj.LF(3,:);
                    rayDir(2,:) = obj.LF(4,:);
                    rayDir(3,:) = 1 - rayDir(1,:).^2 + rayDir(2,:).^2;
                    
                    %create a ray object given this information
                    res = rayC('origin', rayOrigin, 'direction', rayDir, 'waveIndex', obj.get('waveIndex'), 'wave', obj.get('wave'));
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
    end
    
end

