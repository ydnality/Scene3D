classdef VoLTC < clonableHandleObject
    %VOLTC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
       AMatricesUpdated = false; 
    end
    
    properties
        lens; %lens to use for the VoLT model
        film; %see if you can remove this later - a film should not be necessary to calculate this linear transform... think about this...
        
        wave = 400:100:700;
        depths = -100;
        % If fieldPositions is 
        %   1-dimension implies rotationally symmetric, and use of rotation
        % transform 
        %   2-dimension implies not rotationally symemetric 
        % Only 1-dimension currently implemented
        fieldPositions = linspace(.00001,2, 5);  
        % Must have the same dimensions as depths, wavelengths, field position
        ACollection; 
        A1stCollection;
        A2ndCollection;
    end
    
    methods (Access = private)
        function pSLocations = getPSLocations(obj, depthIndex)
            %Returns a collection of point source locations (n x 3, where n
            %is the number of point sources) at a specific depth, specified
            %by depthIndex
            %
            % pSLocations = getPSLocations(obj, depthIndex)
            % depthIndex: the depth index that the user desires to extract
            % PSLocations at.  the resulting collection will assume an x
            % coordinate of 0.  Rotation will be applied later to evaluate
            % PSF's that are off-axis.
            %
            %pSY = 0.01:.3:2;
            %pSZ = -102 * ones(length(pSY), 1);
            %pSLocations = [zeros(length(pSY), 1) pSY' pSZ];
            
            if (ieNotDefined('depthIndex')), depthIndex = 1;
            end
            
            %error checking
            assert(depthIndex > 0 && depthIndex <= length(obj.depths), 'depthIndex out of range');
            
            %old code based off cartesian coordinates
            %pSZ = obj.depths(depthIndex) * ones(length(obj.fieldPositions), 1);
            %pSLocations =  [zeros(length(obj.fieldPositions), 1) obj.fieldPositions'  pSZ];
            
            %here we will transform from spherical to cartesian coordinates
            pSZ = -obj.depths(depthIndex) * cos(obj.fieldPositions/180 * pi)';
            pSLocations =  [zeros(length(obj.fieldPositions), 1) obj.depths(depthIndex) * sin(obj.fieldPositions/180 * pi)'  pSZ];
        end
    end
    
    methods
        
        function obj = VoLTC(varargin)
        % Initialization of the Volume of Linear Transformations object    
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'lens'
                        obj.lens = varargin{ii+1};
                    case 'film'
                        obj.film = varargin{ii+1};  %must be a 2 element vector
                    case 'fieldpositions'
                        obj.fieldPositions = varargin{ii+1};
                    case 'depths'
                        obj.depths = varargin{ii+1};
                    case 'wave'  %somehow syncrhonize this...
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
                case 'lens'
                    res = obj.lens;
                case 'wave'
                    res = obj.wave;
                case 'depths'
                    res = obj.depths;
                case 'numdepths'
                    res = length(obj.depths);
                case 'numfieldpositions'
                    % Number of field heights
                    res = length(obj.fieldPositions);
                case 'acollection'
                    % A transforms from entrance to exit pupil
                    if(obj.AMatricesUpdated)
                        res = obj.ACollection;
                    else
                        error('A Matrices have not been recalculated.  Run VoLTC.calculateMatrices() first');
                    end
                case 'a1stcollection'
                    % A transforms from entrance to middle aperture
                    if(obj.AMatricesUpdated)
                        res = obj.A1stCollection;
                    else
                        error('A Matrices have not been recalculated.  Run VoLTC.calculateMatrices() first');
                    end
                case 'a2ndcollection'
                    % A transforms from middle aperture to exit pupil
                    if(obj.AMatricesUpdated)
                        res = obj.A2ndCollection;
                    else 
                        error('A Matrices have not been recalculated.  Run VoLTC.calculateMatrices() first');
                    end
                case 'pslocations'
                    % VoLT.get('psLocations')
                    if (length(varargin) >= 1)
                        res = obj.getPSLocations(varargin{1});
                    else
                        res = obj.getPSLocations();
                    end
                case 'fieldpositions'
                    res = obj.fieldPositions;
                %TODO (maybe not): have a get A matrix for a particular field pos,
                %depth, wavelength.
                %coefValues dimensions: (4,4, #fieldPositions, #depths, #wavelengths);
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
        
       
    end    
end

