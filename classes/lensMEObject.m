classdef lensMEObject <  handle
    % Create a multiple element lens object
    %
    %   lens = lensMEObject(elOffset, elRadius, elAperture, elN, aperture, focalLength, center);  
    %
    % Distance units, when not specified, are millimeters.
    %
    % Presently we represent multielement lenses as a set of spherical
    % lenses and circular apertures. The multiple elements are defined by a
    % series of surfaces with a curvature, position, and index of
    % refraction (as a function of wavelength).
    %
    % We code the offset to each surface radius to the center of the
    % spherical lens.  Positive means sphere center to the right. Aperture
    % parameters (a single number is a diameter in mm). index of refraction
    % (n) for the material to the left of the surface.
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
    % Example:
    %   lensRealisticObject()
    %   lensObject(2, 67, 3, 1, 3, 50, 0)
    %
    % AL Vistasoft Copyright 2014
    
    properties
        name = 'default';
        type = 'multi element lens';
        surfaceArray;              % Set of spherical surfaces and apertures
        diffractionEnabled = false;% Not implemented yet
        wave = 400:50:700;         % nm
        focalLength = 50;          % mm, focal length of multi-element lens
        apertureMiddleD = 1;       % mm, diameter of the middle aperture
        apertureSample = [11 11];  % Number of spatial samples in the aperture.  Use odd number
    end
    
    % Private functions only called within this class
    methods (Access = private)
        
        function centers = centersCompute(obj, sOffset, sRadius)
            % Used when we read a lens file. The function computes the center
            % position of each spherical lens element from the array of element
            % offsets and radii.
            
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
        
        function elementsSet(obj, sOffset, sRadius, sAperture, sN)
            %Copies the lens elements into the surfaceArray slot of the
            %multielement lens object.
            %
            %The input vectors and should have equal length, and be of type
            %double.
            
            
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
            
            %compute surface array centers
            centers = obj.centersCompute(sOffset, sRadius);
            
            %create array of surfaces
            obj.surfaceArray = lensSurfaceObject();
            
            for i = 1:length(sOffset)
                obj.surfaceArray(i) = ...
                    lensSurfaceObject('sCenter', centers(i, :), ...
                    'sRadius', sRadius(i), ...
                    'apertureD', sAperture(i), 'n', sN(i, :));
            end
            
        end
    end
    
    methods
        %Multiple element lens constructor
        %TODO: error handling
        function obj = lensMEObject(varargin)
            %  surfaceList = lensReadFile(fName);
            %  lensME = lensMEObject('surfaceArray',surfList);
            
            for ii=1:2:length(varargin)
                p = ieParamFormat(varargin{ii});
                switch p
                    case 'surfacearray'
                        obj.surfaceArray = varargin{ii+1};
                    case 'aperturesample'
                        obj.apertureSample = varargin{ii+1};  %must be a 2 element vector
                    case 'aperturemiddled'
                        obj.apertureMiddleD = varargin{ii+1};
                    case 'focallength'
                        obj.focalLength = varargin{ii+1};
                    case 'diffractionenabled'
                        obj.diffractionEnabled = varargin{ii+1};
                    case 'wave'
                        obj.wave = varargin{ii+1};
                    case 'filename'
                        obj.fileRead(varargin{ii+1});
                    otherwise
                        error('Unknown parameter %s\n',varargin{ii});
                end
            end
            %             obj.setElements(elOffset, elRadius, elAperture, elN);
        end
        
        function res = get(obj,pName,varargin)
            % Get various derived lens properties though this call
            pName = ieParamFormat(pName);
            switch pName
                case 'name'
                    res = obj.name;
                case 'type'
                    res = obj.type;
                case 'numels'
                    % Should be nsurfaces
                    res = length(obj.surfaceArray);
                case 'totaloffset'
                    res = -(obj.surfaceArray(1).sCenter(3) - obj.surfaceArray(1).sRadius);
                case 'surfacearray'
                    res = obj.surfaceArray;
                case 'sradius'
                    % spherical radius of curvature of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).sRadius; 
               case 'sdiameter'
                    % Aperture diameter of this surface.
                    % lens.get('sradius',whichSurface)
                    if isempty(varargin), this = 1;
                    else this = varargin{1};
                    end
                    res = obj.surfaceArray(this).apertureD;
                otherwise
                    error('Unknown parameter %s\n',pName);
            end
        end
        
        %%
        function fileRead(obj, fullFileName)
            % Reads PBRT lens file matrix of data
            % The file has focal length information added to the header
            % This function converts the PBRT matrix of data into the format
            % that we use for setting up the multielement lens.
            
            % Open the lens file
            %fid = fopen(fullfile(dataPath, 'rayTrace', 'dgauss.50mm.dat'));
            fid = fopen(fullFileName);
            
            % Read each of the lens and close the file
            import = textscan(fid, '%s%s%s%s', 'delimiter' , '\t');
            fclose(fid);
            
            %first find the start of the lens line demarked by #   radius
            firstColumn = import{1};
            
            %Scan past the comment lines, looking for the data
            % Check the first entry in each line for #
            % When it is not, read the rest of the file as a
            % d = fread(mumble,double)
            % d= reshape(d,X,4);
            %
            % radius axpos N aperture
            % The current test is just to see if the word radius is on the
            % line.  This might be improved.
            continu = true;
            i = 1;
            while(continu && i <= length(firstColumn))
                compare = regexp(firstColumn(i), 'radius');
                if(~(isempty(compare{1})))
                    continu = false;
                end
                i = i+1;
            end
            % i is the row where the data begin
            
            % The next lens are the matrix data
            % put data into lens object
            radius = str2double(import{1});
            radius = radius(i:length(firstColumn));
            
            
            % Change from pbrt Scene3D format to raytrace Scene3D format
            % In PBRT, the row has the offset from the previous surface.  In
            % PBRT the data are read from the bottom up.  The last row has no
            % offset.
            % In PBRT, we trace from the sensor to the scene.
            % In Scene3D we trace from the scene to the sensor.
            % So, the offsets are shifted down.  This means:
            %
            offset = str2double(import{2});
            offset = offset(i:length(firstColumn));
            offset = [0; offset(1:(end-1))];
            
            % Index of refraction in the 3rd column
            N = str2double(import{3});
            N = N(i:length(firstColumn));
            
            % Diameter of the aperture (or maybe radius.  Must determine).
            aperture = str2double(import{4});
            aperture = aperture(i:length(firstColumn));
            
            %modify the object and reinitialize
            obj.elementsSet(offset, radius, aperture, N);
        end
        
        %%
        function numEls = numEls(obj)
            numEls = length(obj.surfaceArray);
        end
        
        %%
        function obj =  draw(obj)
            %draws the illustration of the lens on a figure - you must declare
            %a new graphwin first!
            
            vcNewGraphWin;
            %             prevSurfaceZ = -obj.get('totalOffset');
            %             prevAperture = 1;
            
            for lensEl = 1:obj.numEls
                curEl = obj.surfaceArray(lensEl);
                
                
                %illustrations for debug
                if (curEl.sRadius ~=0)
                    %draw arcs if radius is nonzero
                    nextEl = obj.surfaceArray(min(lensEl+1, end));
                    prevEl = obj.surfaceArray(max(lensEl-1, 1));
                    
                    %lens elements do NOT always end when the neighboring
                    %element begins.  this allows for a fudge factor.  This
                    %won't matter too much because the aperture radius will
                    %be the limiting factor.
                    delta = 10;
                    
                    if (curEl.sRadius > 0 )
                        leftBoundary = curEl.get('zIntercept');
                        rightBoundary = nextEl.get('zIntercept') + delta;
                    else
                        leftBoundary = prevEl.get('zIntercept') - delta;
                        rightBoundary = curEl.get('zIntercept');
                    end
                    zPlot = linspace(leftBoundary, rightBoundary, 1000);
                    
                    yPlot  = sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
                    yPlotN = -sqrt(curEl.sRadius^2 - (zPlot - curEl.sCenter(3)) .^2);
                    
                    
                    %TODO:find a better way to plot the arcs later - this one is prone to potential problem
                    %  withinRange = and(and((yPlot < curEl.apertureD),(zPlot < prevSurfaceZ + curEl.offset + arcZone)), (zPlot > prevSurfaceZ + curEl.offset - arcZone));
                    % withinRange = and(and((yPlot < curEl.apertureD/2),(zPlot <curEl.get('zIntercept') + arcZone)), (zPlot > curEl.get('zIntercept') - arcZone));
                    
                    %may have problems with a concave lens - but these are
                    %rare
                    withinRange = (yPlot < curEl.apertureD/2);
                    line(zPlot(withinRange), yPlot(withinRange));
                    line(zPlot(withinRange), yPlotN(withinRange));
                else
                    %draw the aperture opening if radius = 0
                    
                    %TODO: draw the difference between specified aperture
                    %from file and specified aperture from object instance
                    
                    %right now: take the minimum value
                    curAperture = min(curEl.apertureD/2, obj.apertureMiddleD/2);
                    
                    line(curEl.sCenter(3) * ones(2,1), [-curEl.apertureD/2 -curAperture]);
                    line(curEl.sCenter(3) * ones(2,1), [curAperture curEl.apertureD/2]);
                end
                
                
            end
        end
        
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
        
        function aGrid = fullGrid(obj,randJitter)
            % Build the full sampling grid, possibly adding a little jitter
            % to avoid aliasing artifacts
            
            if (ieNotDefined('randJitter')), randJitter = false; end
            
            % First make the rectangular samples.
            firstApertureRadius = obj.surfaceArray(1).apertureD/2;
            xSamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(1));
            ySamples = linspace(-firstApertureRadius, firstApertureRadius, obj.apertureSample(2));
            [X, Y] = meshgrid(xSamples,ySamples);
            
            %Add a random jitter.  Uniform distribution.  Scales to plus or
            %minus half the sample width
            if(randJitter)
                X = X + (rand(size(X)) - .5) * (xSamples(2) - xSamples(1));
                Y = Y + (rand(size(Y)) - .5) * (ySamples(2) - ySamples(2));
            end
            aGrid.X = X; aGrid.Y = Y;
            
        end
        
        function aGrid = apertureGrid(obj,randJitter)
            % Find the full grid, mask it and return only the (X,Y)
            % positions inside the masked region.  This is the usual set of
            % positions that we use for calculating light fields.
            if ieNotDefined('randJitter'), randJitter = false; end
            
            aGrid = fullGrid(obj,randJitter);
            aMask = apertureMask(obj);
            aGrid.X = aGrid.X(aMask);
            aGrid.Y = aGrid.Y(aMask);
        end
        
        function obj = rtThroughLens(obj, rays, nLines)
            % lens.rtThroughLens(rays,true/false)
            % 
            % Inputs:  lens object,rays at the entrance aperture,
            % calculates the rays at the exit aperture.  The exiting rays
            % are represented by the 
            %
            % On return, the input variable rays, which starts out
            % representing the rays at the entrance aperture, is changed to
            % be the position and direction of the rays at the exit
            % aperture.
            
            % The order is from furthest from film to film, which is also
            % how the rays pass through the optics.
            passedCenterAperture = false;  %true if rays are traced through lens aperture
            if (ieNotDefined('nLines')), nLines = false; end
            
            % prevSurfaceZ = -obj.get('totalOffset');
            prevN = ones(length(rays.origin), 1);
            
            % For each surface element
            for lensEl = 1:obj.numEls
                
                % Get the surface data
                curEl = obj.surfaceArray(lensEl);
                curAperture = curEl.apertureD/2;
                
                %  ----vectorized
                
                %remove the dead rays and use tmp data structure, liveRays
                %                 liveRays = rayObject();
                %                 liveRays.makeDeepCopy(rays);
                %                 deadIndices = isnan(rays.origin);
                %
                %                 liveRays.origin(deadIndices, : ) = [];
                %                 liveRays.direction(deadIndices, : ) = [];
                %                 liveRays.wavelength(deadIndices) = [];
                %                 liveRays.waveIndex(deadIndices) = [];
                %                 liveRays.apertureSamples.X(deadIndices) = [];
                %                 liveRays.apertureSamples.Y(deadIndices) = [];
                %                 liveRays.apertureLocation(deadIndices, :) = [];
                
                %ray trace through a single element
                
                % If spherical, calculate ray intersection with lens element 
                if (curEl.sRadius ~= 0)
                    
                    repCenter = repmat(curEl.sCenter, [length(rays.origin) 1]);
                    repRadius = repmat(curEl.sRadius, [length(rays.origin) 1]);
                    radicand = dot(rays.direction, rays.origin - repCenter, 2).^2 - ...
                        ( dot(rays.origin - repCenter, rays.origin -repCenter, 2)) + repRadius.^2;
                    
                    % Calculate something about the ray angle with respect
                    % to the current surface.  AL to figure this one out
                    % and put in a book reference.
                    if (curEl.sRadius < 0)
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) + sqrt(radicand));
                    else
                        intersectT = (-dot(rays.direction, rays.origin - repCenter, 2) - sqrt(radicand));
                    end
                    
                    %make sure that intersectT is > 0
                    if (min(intersectT(:)) < 0)
                        disp('intersectT less than 0 for lens %d',lensEl);
                    end
                    
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + repIntersectT .* rays.direction;
                    
                    %illustrations for debugging
                    % Let's try to make the lines 3D.  That would be cool!
                    if (nLines)
                        if lensEl == 1
                            vcNewGraphWin; obj.draw();
                            % Find some samples.  Let's get smarter at how
                            % we find them.  For now, pick 20 random
                            % values.
%                             samps = 1:size(rays.origin,1);
%                             nLines = length(samps);
                            samps = randi(size(rays.origin,1),[nLines,1]);
                        end
                        xCoordVector = [rays.origin(samps,3) intersectPosition(samps,3) NaN([nLines 1])]';
                        yCoordVector = [rays.origin(samps,2) intersectPosition(samps,2) NaN([nLines 1])]';
                        xCoordVector = real(xCoordVector(:));
                        yCoordVector = real(yCoordVector(:));
                        line(xCoordVector,  yCoordVector ,'Color','b','LineWidth',1);
                        pause(0.2);
                    end
                    
                    %                         if (isnan(intersectPosition))
                    %                             disp('nan value');
                    %                         end
                else
                    %plane intersection with lens aperture - TODO: maybe put
                    %in function?
                    %assumes that aperture is perfectly perpendicular with optical axis
                    intersectZ = repmat(curEl.sCenter(3), [length(rays.origin) 1]); 
                    intersectT = (intersectZ - rays.origin(:, 3))./rays.direction(:, 3);
                    repIntersectT = repmat(intersectT, [1 3]);
                    intersectPosition = rays.origin + rays.direction .* repIntersectT;
                    curAperture = min(curEl.apertureD, obj.apertureMiddleD)/2;
                    
                    %added for ppsfObject apertureTracking
                    if(isa(rays, 'ppsfObject'))
                        rays.aMiddleInt.XY = 0;
                        rays.aMiddleInt.XY = intersectPosition(:,1:2);  %only X-Y coords
                        rays.aMiddleInt.Z = intersectZ;    %aperture Z
                        passedCenterAperture = true;
                    end
                    %                         if (isnan(intersectPosition))
                    %                             disp('nan value');
                    %                         end
                end
                
                % remove rays that land outside of the aperture
                %TODO: consider making these set functions later
                outsideAperture = intersectPosition(:, 1).^2 + intersectPosition(:, 2).^2 >= curAperture^2;
                rays.origin(outsideAperture, : ) = NaN;
                rays.direction(outsideAperture, : ) = NaN;
                rays.wavelength(outsideAperture) = NaN;
                rays.waveIndex(outsideAperture) = NaN;
                intersectPosition(outsideAperture, :) = NaN;
                prevN(outsideAperture) = NaN;
                
                
                %special case with ppsfObjects
                if(isa(rays,'ppsfObject'))
                    rays.aEntranceInt.XY(outsideAperture, :) = NaN;
                    rays.aMiddleInt.XY(outsideAperture, :) = NaN;                
                    rays.aExitInt.XY(outsideAperture, :) = NaN;    
                end
                
                % Spherical surface so apply Snell's law
                if(curEl.sRadius ~= 0)
                    
                    %in bounds case - perform vector Snell's law
                    repCenter = repmat(curEl.sCenter, [length(rays.origin) 1]);
                    normalVec = intersectPosition - repCenter;  %does the polarity of this vector matter? YES
                    normalVec = normalVec./repmat(sqrt(sum(normalVec.*normalVec, 2)),[1 3]); %normalizes each row
                    
                    if (curEl.sRadius < 0)  %which is the correct sign convention? This is correct
                        normalVec = -normalVec;
                    end
                    
                    %modify the index of refraction depending on wavelength
                    %TODO: have this be one of the input parameters (N vs. wavelength)
                    %                     if (curEl.n ~= 1)
                    %                         curN = (rays.wavelength - 550) * -.04/(300) + curEl.n;
                    %                     else
                    %                         curN = ones(length(rays.wavelength), 1);
                    %                     end
                    
                    %                     waveIndex = find([obj.wave' obj.nWave']== rays.wavelength);  %this doesn't work yet
                    %this was supposed to conver wavelength to waveIndex -
                    %but I coudln't find a way to vectorize it, so instead
                    %we precompute it
                    
                    liveIndices = ~isnan(rays.wavelength);
                    curN = ones(size(prevN));
                    curN(liveIndices) = curEl.n(rays.waveIndex(liveIndices));  %deal with nans
                    curN(~liveIndices) = nan;
                    
                    %                     curN = ones(length(rays.wavelength), 1) * curEl.n;
                    ratio = prevN./curN;    %snell's law index of refraction
                    
                    %Vector form of Snell's Law
                    c = -dot(normalVec, rays.direction, 2);
                    repRatio = repmat(ratio, [1 3]);
                    newVec = repRatio .* rays.direction + repmat((ratio.*c -sqrt(1 - ratio.^2 .* (1 - c.^2))), [1 3])  .* normalVec;
                    newVec = newVec./repmat(sqrt(sum(newVec.*newVec, 2)), [1 3]); %normalizes each row
                    
                    %update the direction of the ray
                    rays.origin = intersectPosition;
                    rays.direction = newVec;
                    prevN = curN;  %note: curN won't change if the aperture is the overall lens aperture
                end
                
                % HURB diffraction calculation
                if (obj.diffractionEnabled)
                    obj.rtHURB(rays, intersectPosition, curEl.aperture);
                end
                
                
                % iterate previous z
                % prevSurfaceZ = prevSurfaceZ + curEl.offset;
            end
        end
        
        function rays = rtSourceToEntrance(obj, pointSource, ppsfObjectFlag)
            % Ray trace from a point to the aperture grid
            %
            % A couple of more comments.  THis is from the source to the
            % first surface, right?
            %
            % See also: rtEntranceToExit in meLens
            
            if (ieNotDefined('ppsfObjectFlag'))
                ppsfObjectFlag = false;
            end
            
            % Define rays object
            if (~ppsfObjectFlag)
                rays = rayObject();
            else
                rays = ppsfObject();
            end
            
            % Create rays from the point source to each aperture grid point
            %the new origin is the position of the current point source
            aGrid = obj.apertureGrid(false);   %TODO: set randJitter somehow
            nPts  = numel(aGrid.X(:));
            aGrid.Z = repmat(-obj.get('totalOffset'),[nPts,1]);
            ePoints = [aGrid.X(:),aGrid.Y(:),aGrid.Z(:)];  %e stands for end... endPoionts
            
            % Computes the directions between the origin and the end points
            rays.origin = repmat(pointSource, [nPts, 1, 1] );
            rays.direction = rayDirection(rays.origin,ePoints);
            
            if(isa(rays,'ppsfObject'))
                rays.aEntranceInt.XY = zeros(length(aGrid.X), 2);
                rays.aEntranceInt.XY(:,1) = aGrid.X(:);
                rays.aEntranceInt.XY(:,2) = aGrid.Y(:);
                rays.aEntranceInt.Z = aGrid.Z(1);
                
                %initialization of these various intersection parameters
                %consider putting in function
                rays.aMiddleInt.XY = zeros(length(aGrid.X), 2);
                rays.aExitInt.XY = zeros(length(aGrid.X), 2);
                
                %                 rays.aEntranceInt.X = aGrid.X(:);
                %                 rays.aEntranceInt.Y = aGrid.Y(:);
            end
            
            
            % rays.direction = [(aGrid.X(:) -  rays.origin(:,1)) (aGrid.Y(:) -  rays.origin(:,2)) (obj.centerPosition(3) - rays.origin (:,3)) .* ones(size(obj.apertureSample.Y(:)))];
            % rays.direction = rays.direction./repmat( sqrt(rays.direction(:, 1).^2 + rays.direction(:, 2).^2 + rays.direction(:,3).^2), [1 3]); %normalize direction
        end
        
        function obj = rtHURB(obj, rays, lensIntersectPosition, curApertureRadius)
            %Performs the Heisenburg Uncertainty Ray Bending method on the
            %rays, given a circular aperture radius, and lens intersection
            %position This function accepts both vector forms of inputs, or
            %individual inputs
            
            % Look for cases when you can use: bsxfun ...
            
            %potentially vectorize later for speed
            %             for i = 1:size(rays.direction, 1)
            %calculate the distance of the intersect point to the center of the lens
            
            %                 curLensIntersectPosition = lensIntersectPosition(i, :);
            
            %we don't care about the z coordinate so we remove it
            %                 curLensIntersectPosition = curLensIntersectPosition(1:2);
            
            %                 curRay.origin = rays.origin(i, :);
            %                 curRay.direction = rays.direction(i, :);
            %                 curRay.wavelength = rays.wavelength(i);
            
            ipLength = sqrt(sum(dot(lensIntersectPosition(:, (1:2)), lensIntersectPosition(:, (1:2)), 2), 2));
            
            %calculate directionS and orthogonal directionL
            
            directionS = [lensIntersectPosition(:, 1) lensIntersectPosition(:,2) zeros(length(lensIntersectPosition), 1)];
            directionL = [-lensIntersectPosition(:,2) lensIntersectPosition(:,1) zeros(length(lensIntersectPosition), 1)];
            
            
            
            normS = repmat(sqrt(sum(dot(directionS, directionS, 2), 2)), [1 3]);
            normL = repmat(sqrt(sum(dot(directionL, directionL, 2), 2)), [1 3]);
            divideByZero = sum(normS,2) ==0;
            
            
            directionS(~divideByZero, :) = directionS(~divideByZero, :)./normS(~divideByZero, :);
            directionL(~divideByZero, :) = directionL(~divideByZero, :)./normL(~divideByZero, :);
            directionS(divideByZero, :) = [ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 2)];
            directionL(divideByZero, :) = [zeros(sum(divideByZero == 1), 1) ones(sum(divideByZero == 1), 1) zeros(sum(divideByZero == 1), 1)];
            
            %                 if (norm(directionS)~= 0)
            %                     directionS = directionS./norm(directionS);
            %                     directionL = directionL./norm(directionL);
            %                 else
            %                     %this is the case that the ray hits the absolute center
            %                     %Then, there is a special case to avoid division by 0
            %                     directionS = [1 0 0];
            %                     directionL = [0 1 0];
            %                 end
            
            pointToEdgeS = curApertureRadius - ipLength;   %this is 'a' from paper  //pointToEdgeS stands for point to edge short
            pointToEdgeL = sqrt((curApertureRadius* curApertureRadius) - ipLength .* ipLength);  %pointToEdgeS stands for point to edge long
            
            lambda = rays.wavelength * 1e-9;  %this converts lambda to meters
            sigmaS = atan(1./(2 * pointToEdgeS *.001 * 2 * pi./lambda));  %the .001 converts mm to m
            sigmaL = atan(1./(2 * pointToEdgeL * .001 * 2 * pi./lambda));
            
            %this function regenerates a 2D gaussian sample and
            %returns it randOut
            %gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseSPointer, noiseLPointer);    %experiment for now
            [randOut] = randn(length(sigmaS),2) .* [sigmaS sigmaL];
            
            %calculate component of these vectors based on 2 random degrees
            %assign noise in the s and l directions according to data at these pointers
            noiseS = randOut(:,1);
            noiseL = randOut(:,2);
            
            %project the original ray (in world coordinates) onto a new set of basis vectors in the s and l directions
            projS = (rays.direction(: , 1) .* directionS(: ,1) + rays.direction(: , 2) .* directionS(:,2))./sqrt(directionS(:,1) .* directionS(:,1) + directionS(:,2) .* directionS(:,2));
            projL = (rays.direction(: , 1) .* directionL(:, 1) + rays.direction(: , 2 ) .* directionL(:,2))./sqrt(directionL(:,1) .* directionL(:,1) + directionL(:,2) .* directionL(:,2));
            thetaA = atan(projS./rays.direction(: , 3));   %azimuth - this corresponds to sigmaS
            thetaE = atan(projL./sqrt(projS.*projS + rays.direction(:, 3).* rays.direction(: , 3)));   %elevation - this corresponds to sigmaL
            
            %add uncertainty
            thetaA = thetaA + noiseS;
            thetaE = thetaE + noiseL;
            
            %convert angles back into cartesian coordinates, but in s,l space
            newprojL = sin(thetaE);
            smallH = cos(thetaE);   %smallH corresponds to the projection of the ray onto the s-z plane
            newprojS = smallH .* sin(thetaA);
            rays.direction(:, 3) = smallH .* cos(thetaA);
            
            %convert from s-l space back to x-y space
            rays.direction(:, 1) = (directionS(:, 1) .* newprojS + directionL(:, 1) .* newprojL)./sqrt(directionS(:,1) .* directionS(:,1) + directionL(:,1) .* directionL(:,1));
            rays.direction(:, 2) = (directionS(:, 2) .* newprojS + directionL(:, 2) .* newprojL)./sqrt(directionS(:,2) .* directionS(:,2) + directionL(:,2) .* directionL(:,2));
            normDirection = repmat(sqrt(sum(dot(rays.direction, rays.direction, 2),2)), [1 3]);
            rays.direction = rays.direction./normDirection;
            
            %reassign ray
            %                 rays.origin(i,:) = curRay.origin;
            %                 rays.direction(i, :) = curRay.direction;
            %                 rays.wavelength(i,:) = curRay.wavelength;
            %             end
        end
        
        
    end
    
end