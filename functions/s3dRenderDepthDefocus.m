function [oi, oiD, D] = s3dRenderDepthDefocus(scene, oi, imgPlaneDist, depthEdges, cAberration)
% Computes the blurred optical image given a scene, optics,a nd  image plane distance
%
% OBSOLETE , but we still use this in oiDepthCompute (AL)
%
%  [oi, oiD, D]=s3dRenderDepthDefocus(scene,oi,imgPlaneDist,depthEdges,cAberration)
%
% Arguments:
%   scene:         ISET scene structure
%   oi:            ISET oi structure
%
%
% IMPLEMENTATION
%  We find the pixels in different depth planes and form an OI of them
%  (oiD). We blur with the appropriate OTF for that depth plane.  We add up
%  the results.
%
%  This algorithm fails to account for depth edges properly. At such edges,
%  light from the background is inappropriately added to the forward
%  occluding plane. Fix this.  Perhaps we blur, find the locations that are
%  spread into the foreground plane, and reduce or zero the photons there.
%
% ALTERNATIVE IMPLEMENTATION
% Another approach is to compute a table of PSFs that depend on wavelength
% and depth.  This could have the form of a 4D matrix in which
% (x,y,lambda,depth). 
%   * The (x,y) dimensions must match the OI spatial dimensions. 
%   * For each scene point we form a (x,y,lambda) matrix that is
% interpolated from the depth dimension.  
%   * We then compute a weighted sum of these, where the weights are the
%   wavelength SPD of the point.  We add these to the output photons.
%
% Example:
% 
% See also:  s3d_Debug.m
%
%
% TODO:
%   This routine is too long. It should rely as much as possible on
%   external functions. and we should probably build the alternative
%   implementation.
%
% Copyright, Stanford, 2011

if ieNotDefined('scene'), error('Scene required'); end
if ieNotDefined('oi'),    error('oi required'); else optics =  oiGet(oi,'optics');end  
if ieNotDefined('imgPlaneDist')
    optics = oiGet(oi,'optics');
    imgPlaneDist = opticsGet(optics,'focal length'); 
end

% Set up the depthCenters and Edges
dMap = sceneGet(scene,'depth map');
if ieNotDefined('depthEdges')
    depthEdges = [min(dMap(:)),max(dMap(:))];
elseif length(depthEdges) == 1
    depthEdges = [min(dMap(:)),depthEdges, max(dMap(:))];
end

depthCenters = zeros(length(depthEdges)-1,1);
for ii=1:(length(depthEdges)-1)
    depthCenters(ii) = depthEdges(ii) + (depthEdges(ii+1) - depthEdges(ii))/2;
end

% The model of chromatic aberration is
%  Defocus(varies SF,constant lambda) + chromatic(constant SF,varies lambda)
%
wave     = sceneGet(scene, 'wave');
if ieNotDefined('cAberration'), cAberration = zeros(length(wave),1); end

% Compute the defocus (diopters, 1/m) at specified depth centers. This
% defocus uses the imgPlaneDist, which may not be equal to the focal
% length.
D = opticsDepthDefocus(depthCenters,optics,imgPlaneDist);
% vcNewGraphWin; 
% plot(depthCenters,D); xlabel('Object distance (m)'); ylabel('Defocus (diopters)')

% Figure out how many deg per samples.  We use the scene, rather than the
% oi, because we haven't computed the oi yet.  The value is the same in
% cpd.
maxSF = sceneGet(scene,'maxfreqres','cpd');
nSteps = min(ceil(maxSF),70);              % Round up, but don't go too high.
sampleSF = linspace(0, maxSF, nSteps);     % cyc/deg
maxVal = 2^16- 1;

%% Attach the optics to the oi
oi = oiSet(oi,'optics',optics);

% Initialize the optic images for each of the depths
oiD    = cell(1,length(depthCenters));
alphaMask = cell(1, length(depthCenters));   %alpha mask for each depth plane, for more realistic blur combination
originalAlphaMask = cell(1, length(depthCenters));


%% From the furthest distance to the nearest distance.
for dd = length(depthCenters):-1:1

    fprintf('Depth %.2f Defocus %.2f (%d of %d)\n',depthCenters(dd),D(dd),dd,length(depthCenters));
    
    % Use defocus of first item in array for all wavelengths
    defocus = cAberration + ones(size(cAberration))*D(dd);
    
    % Here, we build up the OI
    % The returned otf is otf(wave,sf) 
    [otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,defocus);
    % mesh(sampleSFmm,wave,abs(otf))
    
    % Convert the otf into its 2D format.
    optics = opticsBuild2Dotf(optics,otf,sampleSFmm);
    oi = oiSet(oi,'optics',optics);
    % plotOTF(oi,'ls wavelength')
    % plotOTF(oi,'otf wavelength')
    % vcAddAndSelectObject(scene); sceneWindow
    % oi = oiCompute(scene,oi); vcAddAndSelectObject(oi); oiWindow
    
    % Make a scene with only the data from within the depth conditions.
    % This might be
    depthRange = [depthEdges(dd) depthEdges(dd+1)];
    if depthRange(1) == depthRange(2),    sceneD = scene;
    else             sceneD     = sceneDepthRange(scene,depthRange);
    end
    % vcAddAndSelectObject(sceneD); sceneWindow
    % vcAddAndSelectObject(scene); sceneWindow
    
    %-----identify hole regions
    dPlane = (depthEdges(dd) <= dMap) .* (dMap < depthEdges(dd+1)) * 1.0;  
    hRegions = 1.-dPlane;
    alphaMask{dd} = dPlane * 1.0 * maxVal;  %multiply by a constant so that we have less quantization errors
    originalAlphaMask{dd} = alphaMask{dd};
    % perform hole-filling operations
     depthBorder = (imfilter(dPlane * 1.0, [-1 -1 -1; -1 8 -1; -1 -1 -1], 'same', 'symmetric') > 0) * 1.0;
   % depthBorder = (imfilter(dPlane * 1.0, [-1 -1 -1 -1 -1 ; -1 -1 -1 -1 -1; -1 -1 -24 -1 -1; -1 -1 -1 -1 -1; -1 -1 -1 -1 -1], 'same', 'symmetric') > 0) * 1.0;
    borderPixels = sceneGet(sceneD, 'photons') .* repmat(depthBorder, [1 1 31]);
%     borderScene = sceneD;
%     borderScene = sceneSet(borderScene', 'photons', borderPixels);
%     vcAddAndSelectObject(borderScene); sceneWindow;

    filteredBorderPixels = imfilter(borderPixels, fspecial('gaussian', 100,3.5));    %TODO: add the filter in as a parameter!
    filteredWeight = imfilter(depthBorder, fspecial('gaussian', 100,3.5));
    
    inventedPixels = filteredBorderPixels./repmat(filteredWeight, [1 1 31]) .* repmat(hRegions, [1 1 31]);
    inventedPixels(isnan(inventedPixels)) = 0;
    combinedPixels = sceneGet(sceneD, 'photons').* repmat(dPlane, [1 1 31]) + inventedPixels;
    
    filledScene = sceneD;
    filledScene = sceneSet(filledScene, 'photons', combinedPixels);
    %sceneD = filledScene;
%     vcAddAndSelectObject(filledScene); sceneWindow;
    
    %----- blur the depth plane
    oiD{dd} = oiCompute(filledScene,oi);
    
    %----- Compute a transparency for the alphaMask layers.  Transparency will
    % be close to the blurring of the depth plane by the appropriate blur.
    alphaScene = sceneD;
    alphaScene = sceneSet(alphaScene, 'wave', [550]); %set alphaScene to a single wave for faster computation
    alphaScene = sceneSet(alphaScene, 'photons', alphaMask{dd}); 
    %old code (inefficient), but left in for reference
    %alphaScene = sceneSet(alphaScene, 'photons', repmat(alphaMask{dd}, [1 1 sceneGet(alphaScene, 'nwave')]));   %this may be inefficient for now, but we can simplify to 1 channel later for speed
    oiAlphaScene = oiCompute(alphaScene, oi);
    tempAlpha = oiGet(oiAlphaScene, 'photons');
    %tempAlpha =  tempAlpha/max(max(tempAlpha));
    tempAlpha =  tempAlpha/1.8e03;    %this 1.8e03 constant was chosen experimentally to give something that looks good. TODO: more physically accurate methods.  
    tempAlpha(tempAlpha>1) = 1;
    tempAlpha(isnan(tempAlpha)) = 0; 
    alphaMask{dd} = tempAlpha;
    
    % vcAddAndSelectObject(oiD{dd}); oiWindow
    % Then we should combine the oiD structures into the final oi output.
    % That would be outside this loop.  
end

% oiWindow
% for ii=1:length(oiD), vcAddAndSelectObject(oiD{ii}); end; oiWindow

% Combine the multiple defocused optical images into a single OI.
if length(oiD) > 1,  oi = oiCombineDepthsNew(oiD, dMap, depthEdges, depthCenters, alphaMask);
else                 oi = oiD{1};
end

% Should this be in the oiCombineDepths?
oi = oiSet(oi,'illuminance',oiCalculateIlluminance(oi));

return

