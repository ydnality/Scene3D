classdef lensC <  handle
    % Create a multiple element lens object
    %
    %   lens = lensC(elOffset, elRadius, elAperture, elN, aperture, focalLength, center);
    %
    % Distance units, when not specified, are millimeters.
    %
    % Presently we represent multi-element lenses as a set of spherical
    % lenses and circular apertures. The multiple elements are defined by a
    % series of surfaces with a curvature, position, and index of
    % refraction (as a function of wavelength).
    %
    % We code the offset to each surface radius to the center of the
    % spherical lens.  Positive means sphere center to the right. Aperture
    % parameters (a single number is a diameter in mm). index of refraction
    % (n) for the material to the left of the surface.
    %
    % Lens files and surface arrays are specified from the front most
    % element (closest to the scene), to the back-most element (closest to
    % the sensor)
    %
    % pinhole cameras have no aperture and the pinhole lens will inherit
    % this superclas. This will be a superclass that will be inherited by
    % other classes in the future
    %
    % We aim to be consistent with the PBRT lens files, and maybe the Zemax
    % as far possible ?
    %
    % This could become a camera, or we could make a camera object that has
    % a lens and film.
    %
    % Todo: consider a set function for better error handling.
    %
    % Example:
    %   lensC()
    %
    % AL Vistasoft Copyright 2014
    
    properties
        name = 'default';
        type = 'multi element lens';   % We might be removing this soon (why? - BW)
        surfaceArray = surfaceC();     % Set of spherical surfaces and apertures
        diffractionEnabled = false;    % Not implemented yet
        wave = 400:50:700;         % nm
        focalLength = 50;          % mm, focal length of multi-element lens
        apertureMiddleD = 1;       % mm, diameter of the middle aperture
        apertureSample = [11 11];  % Number of spatial samples in the aperture.  Use odd number
        centerZ = 0;               % Theoretical center of lens (length-wise) in the z coordinate
        
        % Black Box Model
        BBoxModel=[]; % Empty
    end
    
    properties (SetAccess = private)
        centerRay = [];   %for use for ideal Lens
        inFocusPosition = [0 0 0];
    end
    methods (Access = public)
        %Multiple element lens constructor
        %TODO: error handling
        function obj = lensC(varargin)
            % Constructor
            
            if isempty(varargin)
                lensFileName = fullfile(s3dRootPath,'data', 'lens', '2ElLens.dat');                
                obj = lensC('apertureSample', [151 151], ...
                    'fileName', lensFileName, ...
                    'apertureMiddleD', 8);
            else
                
                for ii=1:2:length(varargin)
                    p = ieParamFormat(varargin{ii});
                    switch p
                        case 'name'
                            obj.name = varargin{ii+1};
                        case 'type'
                            obj.type = varargin{ii+1};
                        case 'surfacearray'
                            obj.surfaceArray = varargin{ii+1};
                        case 'aperturesample'
                            obj.apertureSample = varargin{ii+1};  %must be a 2 element vector
                        case 'aperturemiddled'
                            obj.apertureMiddleD = varargin{ii+1};
                        case 'apertureindex'
                            obj.apertureIndex(varargin{ii+1});
                        case 'focallength'
                            obj.focalLength = varargin{ii+1};
                        case 'diffractionenabled'
                            obj.diffractionEnabled = varargin{ii+1};
                        case 'wave'
                            obj.set('wave', varargin{ii+1});
                        case 'filename'
                            obj.fileRead(varargin{ii+1});
                            
                        case {'blackboxmodel';'blackbox';'bbm'} % equivalent BLACK BOX MODEL
                            obj.BBoxModel = varargin{ii+1};
                            
                        otherwise
                            error('Unknown parameter %s\n',varargin{ii});
                    end
                end
            end
            
        end
        
        % Get properties
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'name'
                    res = obj.name;
                case 'wave'
                    res = obj.wave;
                case 'nwave'
                    res = length(obj.wave);
                case 'type'
                    res = obj.type;
                case {'nsurfaces','numels'}
                    % Should be nsurfaces
                    res = length(obj.surfaceArray);
                case 'totaloffset'
                    % This is the size (in mm) from the front surface to
                    % the back surface.  The last surface is at 0, so the
                    % position of the first surface is the total size.
                    sArray = obj.surfaceArray;
                    res    = -1*sArray(1).get('zpos');
                    % res = -(obj.surfaceArray(1).sCenter(3) - obj.surfaceArray(1).sRadius);
                case 'offsets'
                    % Offsets format (like PBRT files) from center/zPos
                    % data
                    res = obj.offsetCompute();
                case 'surfacearray'
                    % lens.get('surface array',[which surface])
                    if isempty(varargin), res = obj.surfaceArray;
                    else                  res = obj.surfaceArray(varargin{1});
                    end
                    
                case {'indexofrefraction','narray'}
                    nSurf = obj.get('nsurfaces');
                    sWave  = obj.surfaceArray(1).wave;
                    res = zeros(length(sWave),nSurf);
                    for ii=1:nSurf
                        res(:,ii) = obj.surfaceArray(ii).n(:)';
                    end
                case {'refractivesurfaces'}
                    % logicalList = lens.get('refractive surfaces'); 
                    % Returns
                    %  1 at the positions of refractive surfaces, and
                    %  0 at diaphgrams
                    nSurf = obj.get('nsurfaces');
                    res = ones(nSurf,1);
                    for ii=1:nSurf
                        if strcmp(obj.surfaceArray(ii).subtype,'diaphragm')
                            res(ii) = 0;
                        end
                    end
                    res = logical(res);
                case {'nrefractivesurfaces'}
                    % nMatrix = lens.get('n refractive surfaces')
                    %
                    % The refractive indices for each wavelength of each
                    % refractive surface.  The returned matrix has size
                    % nWave x nSurface
                    lst = find(obj.get('refractive surfaces'));
                    nSurfaces = length(lst);
                    nWave = obj.get('nwave');
                    res = zeros(nWave,nSurfaces);
                    for ii = 1:length(lst)
                        res(:,ii) = obj.surfaceArray(lst(ii)).n(:);
                    end
                    
                case 'sradius'
                    % spherical radius of curvature of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).sRadius;
                case 'sdiameter'
                    % lens.get('s diameter',nS);
                    % Aperture diameter of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).apertureD;
                case {'aperture','diaphragm'}
                    % lens.get('aperture')
                    % Returns the surface number of the aperture
                    % (diaphragm)
                    s = obj.surfaceArray;
                    for ii=1:length(s)
                        if strcmp(s(ii).subtype,'diaphragm')
                            res = ii;
                            return;
                        end
                    end
                case 'aperturesample'
                    res =  obj.apertureSample ;    
                case {'middleapertured','aperturemiddled'}
                    % The middle aperture is the diameter of the diaphragm,
                    % which is normally the middle aperture.  We should
                    % change this somehow for clarity.  Or we should find
                    % the diaphragm and return its diameter.
                    %
                    % The diameter of the middle aperture
                    % units are mm
                    res = obj.apertureMiddleD;
                    
                case {'blackboxmodel';'blackbox';'bbm'} % equivalent BLACK BOX MODEL
                    param=varargin{1};  %witch field of the black box to get
                    res = obj.bbmGetValue(param);
                    
                case {'lightfieldtransformation';'lightfieldtransf';'lightfield'} % equivalent BLACK BOX MODEL
                    if nargin >2
                        param = varargin{1};  %witch field of the black box to get
                        param = ieParamFormat(param);
                        switch param
                            case {'2d'}
                                res = obj.bbmGetValue('abcd');
                            case {'4d'}
                                abcd = obj.bbmGetValue('abcd');
                                nW=size(abcd,3);
                                dummy=eye(4);
                                for li=1:nW
                                    abcd_out(:,:,li)=dummy;
                                    abcd_out(1:2,1:2,li)=abcd(:,:,li);
                                end
                                res=abcd_out;                                
                            otherwise
                                error(['Not valid :',param ,' as type for  Light Field tranformation']); 
                        end
                    else
                        res = obj.bbmGetValue('abcd');
                    end
                        
                
                case {'opticalsystem'; 'optsyst';'opticalsyst';'optical system structure'} 
                    % Get the equivalent optical system structure generated
                    % by Michael's script      
                    % Can be specify refractive indices for object and
                    % image space as varargin {1} and {2}
                    if nargin >2
                        n_ob=varargin{1};    n_im=varargin{2};
                        OptSyst=obj.bbmComputeOptSyst(n_ob,n_im);
                    else
                        OptSyst=obj.bbmComputeOptSyst();
                    end
                    res = OptSyst;
                    
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
        
        %% Set

        function set(obj,pName,val,varargin)
            pName = ieParamFormat(pName);
            switch pName
                case 'wave'
                    % lens.set('wave',val);
                    % The wavelength is annoying.
                    % This could go away some day.  But for now, the wavelength set is
                    % ridiculous because there are so many copies of wave.  So, we
                    % should set it here rather than addressing every surface element.
                    % (BW)
                    obj.wave = val;
                    nSurfaces = obj.get('n surfaces');
                    for ii=1:nSurfaces
                        obj.surfaceArray(ii).set('wave', val);
                    end
                case 'surfacearray'
                    % lens.set('surface array',val);
                    obj.surfaceArray = val;
                case 'aperturesample'
                    obj.apertureSample = val; 
                case 'surfacearrayindex'
                    index=varargin{1};
                    % lens.set('surface array',val);
                    obj.surfaceArray(index) = val;                    
                case 'apertureindex'
                    % lens.set('aperture index',val);
                    % Indicates which of the surfaces is the aperture.
                    index=val;
                    obj.apertureIndex(index); % Set the surface (specify by varargin) as aperture                    
                case 'nall'
                    % Set the index of refraction to all the surfaces
                    nSurfaces = obj.get('n surfaces');
                    for ii=1:nSurfaces
                        obj.surfaceArray(ii).n = val;
                    end
                case 'nrefractivesurfaces'
                    % lens.set('n refractive surfaces',nMatrix)
                    %
                    % The nMatrix should be nWave x nSurfaces where
                    % nSurfaces are only the refractive surfaces.  This
                    % sets the index of refraction to all those surfaces
                    
                    % Indices of the refractive surfaces
                    lst = find(obj.get('refractive surfaces'));  
                    
                    % For every column in val, put it in the next
                    % refractive index of the lst.
                    kk = 0;   % Initiate the counter
                    for ii = 1:length(lst)
                        kk = kk + 1;
                        obj.surfaceArray(lst(ii)).n(:) = val(:,kk);
                    end
                    
                case {'effectivefocallength';'efl';'focalradius';'imagefocalpoint';...
                        'objectfocalpoint';'imageprincipalpoint';'objectprincipalpoint';...
                        'imagenodalpoint';'objectnodalpoint';'abcd';'abcdmatrix'}
                    % Build the field to append
                    obj.bbmSetField(pName,val);
                    
                case {'blackboxmodel';'blackbox';'bbm'}
                    % Get the parameters from the optical system structure
                    % to build an  equivalent Black Box Model of the lens.
                    % The OptSyst structure has to be built with the
                    % function 'paraxCreateOptSyst' Get 'new' origin for
                    % optical axis
                    OptSyst=val;
                    %                     z0 = OptSyst.cardPoints.lastVertex;
                    z0=paraxGet(OptSyst,'lastVertex');
                    % Variable to append
                    %                     efl=OptSyst.cardPoints.fi; %focal lenght of the system
                    efl=paraxGet(OptSyst,'efl');
                    obj.bbmSetField('effectivefocallength',efl);
                    %                     pRad = OptSyst.Petzval.radius; % radius of curvature of focal plane
                    pRad = paraxGet(OptSyst,'focalradius'); % radius of curvature of focal plane
                    obj.bbmSetField('focalradius',pRad);
                    %                     Fi=OptSyst.cardPoints.dFi;     %Focal point in the image space
                    Fi= paraxGet(OptSyst,'imagefocalpoint')-z0;     %Focal point in the image space
                    obj.bbmSetField('imagefocalpoint',Fi);
                    %                     Hi=OptSyst.cardPoints.dHi; % Principal point in the image space
                    Hi= paraxGet(OptSyst,'imageprincipalpoint')-z0; % Principal point in the image space
                    obj.bbmSetField('imageprincipalpoint',Hi);
                    %                     Ni=OptSyst.cardPoints.dNi;     % Nodal point in the image space
                    Ni=paraxGet(OptSyst,'imagenodalpoint')-z0;    % Nodal point in the image space
                    obj.bbmSetField('imagenodalpoint',Ni);
                    %                     Fo=OptSyst.cardPoints.dFo-z0; %Focal point in the object space
                    Fo=paraxGet(OptSyst,'objectfocalpoint')-z0; %Focal point in the object space
                    obj.bbmSetField('objectfocalpoint',Fo);
                    %                     Ho=OptSyst.cardPoints.dHo-z0; % Principal point in the object space
                    Ho=paraxGet(OptSyst,'objectprincipalpoint')-z0; % Principal point in the object space
                    obj.bbmSetField('objectprincipalpoint',Ho);
                    %                     No=OptSyst.cardPoints.dNo-z0; % Nodal point in the object space
                    No=paraxGet(OptSyst,'objectnodalpoint')-z0; % Nodal point in the object space
                    obj.bbmSetField('objectnodalpoint',No);
                    M = paraxGet(OptSyst,'abcd'); % The 4 coefficients of the ABCD matrix of the overall system
                    obj.bbmSetField('abcd',M);
                    
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
            
        end
        
    end
    
    
    % Not sure why these are private. (BW).
    methods (Access = private)
        
        function apertureMask = apertureMask(obj)
            % Identify the grid points on the circular aperture.
            %
            
            % Here is the full sampling grid for the resolution and
            % aperture radius
            aGrid = obj.fullGrid;
            
            % We assume a circular aperture. This mask is 1 for the sample
            % points within a circle of the aperture radius
            firstApertureRadius = obj.surfaceArray(1).apertureD/2;
            apertureMask = (aGrid.X.^2 + aGrid.Y.^2) <= firstApertureRadius^2;
            % vcNewGraphWin;  mesh(double(apertureMask))
            
        end
        
        function aGrid = fullGrid(obj,randJitter, rtType)
            % Build the full sampling grid, possibly adding a little jitter
            % to avoid aliasing artifacts
            
            if (ieNotDefined('randJitter')), randJitter = false; end
            if (ieNotDefined('rtType')), rtType = 'realistic'; end
            
            % If an ideal rt type, use the middle aperture as the front
            % aperture because there will really be only 1 aperture, the
            % middle one.
            if (strcmp(rtType, 'ideal'))
                firstApertureRadius = obj.apertureMiddleD/2;
            else
                firstApertureRadius = obj.surfaceArray(1).apertureD/2;
            end
            
            % First make the rectangular samples.
            xSamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(1));
            ySamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(2));
            [X, Y] = meshgrid(xSamples,ySamples);
            
            %Add a random jitter.  Uniform distribution.  Scales to plus or
            %minus half the sample width
            if(randJitter)
                X = X + (rand(size(X)) - .5) * (xSamples(2) - xSamples(1));
                Y = Y + (rand(size(Y)) - .5) * (ySamples(2) - ySamples(1));
            end
            aGrid.X = X; aGrid.Y = Y;
            
        end
        
        function aGrid = apertureGrid(obj,randJitter, rtType)
            % Find the full grid, mask it and return only the (X,Y)
            % positions inside the masked region.  This is the usual set of
            % positions that we use for calculating light fields.
            if ieNotDefined('randJitter'), randJitter = false; end
            
            aGrid = fullGrid(obj,randJitter, rtType);
            aMask = apertureMask(obj);
            aGrid.X = aGrid.X(aMask);
            aGrid.Y = aGrid.Y(aMask);
            
            %debug check
            % figure; plot(aGrid.X, aGrid.Y, 'o');
        end
        
        function centers = centersCompute(obj, sOffset, sRadius)
            % Computes the center position of each spherical lens element
            % from the array of element offsets and radii.
            %
            % N.B.  obj is not used here.  Maybe this shouldn't be in
            % the class definition?
            nSurfaces = length(sOffset);
            
            % The lens file includes the offsets between each surface and its
            % previous one (millimeters).  The location of the first surface is
            % at negative of the total offset.
            zIntercept = zeros(nSurfaces,1);
            zIntercept(1) = -sum(sOffset);
            
            % Cumulte the offsets to find the intercept of the following
            % surface
            for ii=2:length(sOffset)
                zIntercept(ii) = zIntercept(ii-1) + sOffset(ii);
            end
            
            % Centers of the spherical surfaces are stored in this matrix
            % (nSurfaces x 3)
            z = zeros(nSurfaces,1);
            centers = [z z zIntercept + sRadius];
        end

                
        function offsets = offsetCompute(obj)
            % Computes offsets from the object's z positions

            %get the radii
            nEls = obj.get('nsurfaces');
            zPos = zeros(1, nEls);
            offsets = zeros(1,nEls);
            sArray = obj.get('surfaceArray');

            for i = 1:nEls
                zPos(i) = sArray(i).get('zpos');
                if (i > 1)
                    offsets(i-1) = zPos(i) - zPos(i-1);
                else
                    offsets(i) = 0;
                end
            end
        end
        
        function elementsSet(obj, sOffset, sRadius, sAperture, sN)
            %Copies the lens elements into the surfaceArray slot of the
            %multielement lens object.
            %
            %The input vectors and should have equal length, and be of type
            %double.
            %
            % Should this just be called from lens.set('surface array')?
            
            
            % Check argument size
            if (length(sOffset)     ~= length(sRadius) || ...
                    length(sOffset) ~= length(sAperture) || ...
                    length(sOffset) ~= size(sN,1))
                error('Input vectors must be of equal length');
            end
            
            % If no wavelength dependence of index of refraction specified,
            % we assume constant across all measurement wavelengths
            % sN is nSurfaces x nWave
            if (size(sN,2) == 1)
                sN = repmat(sN, [1 length(obj.wave)]);
            end
            
            
            %create array of surfaces
            obj.surfaceArray = surfaceC();
            
            %compute surface array centers
            centers = obj.centersCompute(sOffset, sRadius);
            
            for i = 1:length(sOffset)
                obj.surfaceArray(i) = ...
                    surfaceC('sCenter', centers(i, :), ...
                    'sRadius', sRadius(i), ...
                    'apertureD', sAperture(i), 'n', sN(i, :));
            end
            
        end
        
        
        
        function obj = rtIdealThroughLens(obj, rays, nLines)
            % Traces lens.rays through the ideal lens - what is that?
            % This is different from the Realistic one above, but we should
            % figure out how to make them pretty much the same.
            
            lWidth = 0.5; lColor = [0 0.5 1];
            
            %---------------- lens refraction code  -----
            %when intersecting ideal lens, change the direction to intersect the
            %inFocusPosition, and update the origin
            lensIntersectT = (obj.centerZ - rays.origin(:,3))./ rays.direction(:,3);
            lensIntersectPosition = rays.origin +  repmat(lensIntersectT, [1 3]) .* rays.direction;
            
            % Debug visualization
            if (nLines)
                vcNewGraphWin;
                samps = randi(size(rays.origin,1),[nLines,1]);
                
                xCoordVector = [rays.origin(samps,3) lensIntersectPosition(samps,3) NaN([nLines 1])]';
                yCoordVector = [rays.origin(samps,2) lensIntersectPosition(samps,2) NaN([nLines 1])]';
                xCoordVector = real(xCoordVector(:));
                yCoordVector = real(yCoordVector(:));
                line(xCoordVector,  yCoordVector ,'Color',lColor,'LineWidth',lWidth);
                pause(0.2);
                
            end
            
            
            %calculate new direction
            %             newRays = rayObject(); % added
            rays.origin = lensIntersectPosition;
            rays.direction = repmat(obj.inFocusPosition , [size(rays.origin,1) 1 ]) - rays.origin;
            %rays.wavelength = rays.wavelength;
            
            % diffraction HURB calculation
            if (obj.diffractionEnabled)
                obj.rtHURB(rays, lensIntersectPosition, obj.apertureMiddleD/2);
            end
            %---------------- lens refraction code  -----
            
        end
    end
    
end
