function oi = s3dRenderOIAndDepthMap(pbrt, oiName, dockerFlag)
%This function renders an oi AND the depth map, given a pbrt object.
%
%   oi = s3dRenderOIAndDepthMap(pbrt, oiName, dockerFlag)
%
%  pbrt:        The pbrt structure set up elsewhere
%  focalLength: In millimeters of the lens assumed in the oi
%  oiName:      Output name of the oi
%  dockerFlag:  Run as a docker container from docker hub.
%
% See also:  s_s3dRenderHDRBenchLF.m
%
% AL, VISTASOFT, 2014

%% Input argument checking
if (ieNotDefined('pbrt')),        error('pbrt object required.'); end
if (ieNotDefined('dockerFlag')),  dockerFlag = true; end
if (ieNotDefined('oiName')),      oiName = 'unamedScene'; end

%%
if (isa(pbrt, 'pbrtObject'))
    %% render oi irradiance
    
    %Must make a copy because the docker container will overwrite the
    %object.
    radianceRenderPbrt = pbrtObject;
    radianceRenderPbrt.makeDeepCopy(pbrt);
    
    oi = s3dRenderOI(radianceRenderPbrt, oiName, dockerFlag);
    % vcAddObject(oi); oiWindow;
    
    %% Render Depth map
    %change the sampler to stratified for non-noisy depth map
    depthRenderPbrt = pbrtObject; depthRenderPbrt.makeDeepCopy(pbrt);
    numRenders = 1;
    groundTruthDepthMap = s3dRenderDepthMap(depthRenderPbrt, numRenders, oiName, dockerFlag);
    oi = oiSet(oi, 'depth map', groundTruthDepthMap);
    % vcAddObject(oi); oiWindow;

elseif (ischar(pbrt))
    % Renders from a pbrt text file
    oi = s3dRenderOI( pbrt, oiName, dockerFlag);
    [directory, fileName, extension] = fileparts(pbrt);
    
    %depth map pbrt file must have a _depth appended to name
    depthPbrtFile = fullfile(directory, [fileName '_depth', extension]);
    numRenders = 1;
    groundTruthDepthMap = s3dRenderDepthMap(depthPbrtFile, numRenders, oiName, dockerFlag);
    oi = oiSet(oi, 'depthmap', groundTruthDepthMap);
else
    error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject');
end

end