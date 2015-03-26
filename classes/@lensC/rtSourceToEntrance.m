function rays = rtSourceToEntrance(obj, pointSource, ppsfCFlag, jitterFlag, rtType, subSection, depthTriangles)
% Ray trace from a point to the aperture grid on the first lens surface
% (the one furthest from the sensor).
%
%  rays = rtSourceToEntrance(obj, pointSource, ppsfCFlag, jitterFlag, rtType)
%
% Rays have 1 wavelength assigned to them.  There is no particular
% wavelength dependence in air, so there is no need to have multiple
% indices of refraction or wavelength at this moment.
%
% The rays will be "expanded" later to save on computations to
% handle wavelength differences.
%
% The object is a lensC
%     pointSource    -  a 3 vector
%     ppsfCFlag -  specifies whether we are computing the
%     plenoptic PSF (pPSF) or not.  Default: false (PSF only)
%
% See also: rtEntranceToExit

%% Parameter checking

if ieNotDefined('ppsfCFlag'), ppsfCFlag = false; end
if ieNotDefined('jitterFlag'),     jitterFlag = false;     end
if ieNotDefined('rtType'),         rtType = 'realistic';   end
if ieNotDefined('subSection'),      subSection = [];    end
if ieNotDefined('depthTriangles'),    depthTriangles = []; end

% Could pass as an argument
debugOn = false;

% Define rays object
if (~ppsfCFlag), rays = rayC();
else             rays = ppsfC();
end

%% Find distance to first surface for different ray trace types
rtType = ieParamFormat(rtType);
switch rtType
    case 'realistic'
        frontRTSurface = -obj.get('totalOffset');
    case 'ideal'
        frontRTSurface = obj.centerZ;
        
        % --- center ray calculation-------------
        
        % Isn't this in a function somewhere?  If not, let's put it in.
        
        %trace ray from point source to lens center, to image.
        %This helps determine the point of focus
        obj.centerRay.origin = pointSource;
        
        % This could be a call to rayDirection
        obj.centerRay.direction = [ 0 0 obj.centerZ] - obj.centerRay.origin;
        obj.centerRay.direction = obj.centerRay.direction./norm(obj.centerRay.direction);
        
        %calculate the z-position of the in-focus plane using
        %thin lens equation.
        inFocusDistance = 1/(1/obj.focalLength - -1/pointSource(3));
        
        % Calculates the 3-vector for the in-focus position.
        % The in-focus position is the intersection of the
        % in-focus plane and the center-ray
        inFocusT = (inFocusDistance - obj.centerRay.origin(3))/obj.centerRay.direction(3);
        obj.inFocusPosition = obj.centerRay.origin + inFocusT .* obj.centerRay.direction;
        
        % --------center ray calculation -------
    case 'linear'
        % VoLT Method.  Not yet implemented or maybe not needed.  AL to
        % check.
        error('Linear method not implemented yet')
    otherwise
        error('Unknown ray trace method:  %s\n',rtType)
end


%% Do slightly different things for full or partial aperture.
if (isempty(subSection))
    % Case 1: full aperture is sampled
    
    % Create rays from the point source to each aperture grid point
    % The new origin is the position of the current point source
    % The jitterFlag scatters the points on the front surface a bit
    % to improve the rendering.
    aGrid   = obj.apertureGrid(jitterFlag, rtType);
else
    % Case 2: only a subSection of the aperture is sampled.
    
    % for now, it's a sub section of the aperture, defined by 'subSection',
    % which is a 4 element vector, providing the lowerLeft, and upperRight
    % boundaries
    %
    % [leftX lowerY rightX upperY]
    
    leftX = subSection(1);
    lowerY = subSection(2);
    rightX = subSection(3);
    upperY = subSection(4);
    
    % make the rectangular samples.
    firstApertureRadius = obj.surfaceArray(1).apertureD/2;
    xSamples = linspace(firstApertureRadius *leftX, firstApertureRadius*rightX, obj.apertureSample(1));
    ySamples = linspace(firstApertureRadius *lowerY, firstApertureRadius*upperY, obj.apertureSample(2));
    [X, Y] = meshgrid(xSamples,ySamples);
    
    %Add a random jitter.  Uniform distribution.  Scales to plus or
    %minus half the sample width
    if(jitterFlag)
        X = X + (rand(size(X)) - .5) * (xSamples(2) - xSamples(1));
        Y = Y + (rand(size(Y)) - .5) * (ySamples(2) - ySamples(1));
    end
    aGrid.X = X(:); aGrid.Y = Y(:);
end

%% Set Z to the position of the front surface
nPts    = numel(aGrid.X(:));
aGrid.Z = repmat(frontRTSurface,[nPts,1]);

% These are the end points of the ray in the aperture plane
ePoints = [aGrid.X(:),aGrid.Y(:),aGrid.Z(:)];

% These are the directions from the point source to the end
% points in the aperture
rays.origin    = repmat(pointSource, [nPts, 1, 1] );
rays.direction = rayDirection(rays.origin,ePoints);
rays.distance  = ptDistances(pointSource,ePoints);


%% Project occlude - Experimental code. Does not run normally
%
% Need to deal with removing dead rays differently below.
%
% Some rays will be occluded from the aperture by scene objects.  If we
% have a triangle mesh, then we can check for which rays will be occluded
% by the objects.  We do that here.

% If this variable is included, we process for scene depth occlusions
if(~ieNotDefined('depthTriangles'))
%if(false)
    
    %setup origin and direction of rays
    epsilon = repmat([0 0 2], [size(rays.origin, 1) 1]);
    orig = single(rays.origin + epsilon);
    dir  = single(ePoints - rays.origin);
    
    vert1 = single(depthTriangles.vert1);
    vert2 = single(depthTriangles.vert2);
    vert3 = single(depthTriangles.vert3);
    
    %debug visualization
    newE = orig + dir;    
   % vcNewGraphWin;  
    hold on;
    samps = 1:10:size(orig(:,1));
    nSamps = length(samps);
    line([orig(samps,1) newE(samps,1) NaN(nSamps, 1)]',  [orig(samps,2) newE(samps,2) NaN(nSamps, 1)]', [orig(samps,3) newE(samps,3) NaN(nSamps, 1)]');
    
    % Setup for the function TriangleRayIntersection intersections with all
    % combinations of rays and triangles We need to provide every
    % combination of line and triangles.
    %
    % We need to make this faster and simpler and perhaps parallel.  It
    % seems to run out of memory.
    % 
    % Repeat origin so can potentially intersect with all triangles
    origExp = repmat(orig, [size(vert1,1) 1]);  
    dirExp = repmat(dir, [size(vert1,1) 1]);    %cycles once before repeating ... ie 1 2 3 4 5 1 2 3 4 5
    
    %repeats each triangle size(orig) times
    vert1Exp = repmat(vert1(:), [1 size(orig, 1)])';    
    vert1Exp = reshape(vert1Exp, size(dirExp,1), [] );     %ex.  1 1 2 2 3 3 4 4 ... the 1st one represents potential intersection with 1st ray... 2nd repetition intersection with 2nd ray etc.
    
    vert2Exp = repmat(vert2(:), [1 size(orig, 1)])';
    vert2Exp = reshape(vert2Exp, size(dirExp,1), [] );
    
    vert3Exp = repmat(vert3(:), [1 size(orig, 1)])';
    vert3Exp = reshape(vert3Exp, size(dirExp,1), [] );
    
    %perform intersection
    [intersect,~,~,~,xcoor] = TriangleRayIntersection(origExp, dirExp, ...
        vert1Exp, vert2Exp, vert3Exp, 'lineType', 'segment'); %, 'lineType' , 'line');
    
    %separate out and figure out which rays intersected
    % index = 1:length(intersect);
    % blockedIndex = index' .* intersect;
    % blockedIndex(blockedIndex==0) = [];  %remove 0's
    % blockedRays = mod(blockedIndex, size(orig,1));  %do a mod
    % blockedRays(blockedRays ==0) = size(orig, 1);  % so 0's are length rays instead of 0
    % %blockedRays = blockedRays + 1;  %add 1 so index starts at 1 instead of 0
    
    %figure out which rays intersected
    index = repmat(1:size(orig, 1), [1 size(vert1,1)]);  %is this line wrong??
    blockedRays = index(intersect);
    
    %if(false)
    if (~isempty(blockedRays))
        
        %debug visualization
        
        debugOn = true;
        if (pointSource(2) < 10 && debugOn)
            disp('blocked Rays');
            %plot origin
            hold on;
            scatter3(origExp(1,1), origExp(1,2), origExp(1,3), 100, 'r', 'o', 'filled');
            %plot intersections
            hold on;
            scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o',  'filled');
            %plot end points
            hold on;
            scatter3(ePoints(:,1), ePoints(:,2), ePoints(:,3) , 'r', 'o',   'filled');
            
            %plot blocked rays
            hold on;
            endP = rays.origin(blockedRays,:) + rays.direction(blockedRays, :) * 10;
            scatter3(endP(:,1), endP(:,2), endP(:,3), 'r', 'o', 'filled');
        end
        
        %    Remove blocked rays from the ray bundle
        rays.origin(blockedRays, :) = [];
        rays.direction(blockedRays, :) = [];   %perhaps make this more elegant in the future...
        
    end
    

end


%% The ray is of type ppsf.  These are properties it will need 

if(isa(rays,'ppsfC'))
      
    % If the rays are a plenoptic point spread, AL thinks we
    % should store the XY positions at the front aperture, middle aperture,
    % and exit aperture.  The slots for this information are created and
    % initialized here.
    %
    % BW isn't sure why this is useful, but he believes it.
    rays.aEntranceInt.XY = zeros(length(aGrid.X), 2);
    rays.aEntranceInt.XY(:,1) = aGrid.X(:);
    rays.aEntranceInt.XY(:,2) = aGrid.Y(:);
    rays.aEntranceInt.Z = aGrid.Z(1);
    
    % The intersection positions (XY) of the rays in the middle
    % aperture and the exit aperture will be stored here.
    rays.aMiddleInt.XY = zeros(length(aGrid.X), 2);
    rays.aExitInt.XY = zeros(length(aGrid.X), 2);
    
    %if (false)
    if(~ieNotDefined('depthTriangles'))        
        if (~isempty(blockedRays > 0))
            rays.aEntranceInt.XY(blockedRays,:) = [];
            rays.aMiddleInt.XY(blockedRays,:) = [];
            rays.aExitInt.XY(blockedRays,:) = [];
        end
    end
    
end

end