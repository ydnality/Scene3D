%% s_fourFlashRenderForwardModel
%This script models the 4 flash problem as a forward model.  In other
%words, we are given our data, and we will try to reproduce the 4 flash
%images given our input data (depth, normals, reflectance, light)

%% first render the images, depth map, and normals
s_fourFlashRenderTestSpheresWideFOV;

% important output: groundTruthDepthMap, normalVector

%% calculate light vectors
% assume flash is in line with camera for now

% assume 1 for reflectance, and 1 channel for now
reflectance = ones(size(groundTruthDepthMap));

% calculate light vectors (assumes square image for now)
xVector = linspace(-fieldOfView/2, fieldOfView/2, size(groundTruthDepthMap,2));
xVector = tan(xVector * pi/180);
yVector = -linspace(-fieldOfView/2, fieldOfView/2, size(groundTruthDepthMap,1));
yVector = tan(yVector * pi/180);

% put into a grid
[xGrid yGrid] = meshgrid(xVector, yVector);
zGrid = ones(size(xGrid));

% concatenate and normalize
lightVector = cat(3, xGrid, yGrid, zGrid);
lightVector = normvec(lightVector, 'p', 2, 'dim', 3);
scaledLightVector = lightVector/2 + .5;

vcNewGraphWin; imshow(scaledLightVector);  %**normalize the vectors first

%%  calculate dot product between normals and light vector


dotProduct = dot(lightVector, normalVector, 3);
vcNewGraphWin; imshow(dotProduct);

%% put in depth information

%intensity drops off as 1/depth^2

depthFactor = 1./groundTruthDepthMap.^2;

finalRender = (dotProduct .* depthFactor .* reflectance);

% normalize
finalRender = finalRender/max(finalRender(:));
figure; imshow(finalRender * .5 + .5);