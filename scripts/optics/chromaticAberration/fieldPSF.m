%% Renders a point spread function at different field locations
% Andy Lin
% 
% Also renders this same scene at different depths


reRenderScenes = true;
inputFile = 'testTargetCA/realisticStarField9Offset.pbrt';

%%% Scene rendering
%%% Ambient only scene generation
sceneName = 'Starfield Image';
oi = s3dRenderOI(inputFile, .050);
oi = oiSet(oi, 'name', sceneName);
% strip the file name from the path and assign that as the name of the
% object  ... vcSaveObject(oi,);
vcAddAndSelectObject(oi);
oiWindow;


%% generate a starfield image for use as a texture

imageSize = [800 800];
horSpacing = 50;
verSpacing = 50;
midPoint = round(imageSize/2);
xStart = 1;
yStart = 1;

sFImage = zeros(imageSize);
sFImage(yStart:verSpacing:imageSize(1), xStart:horSpacing:imageSize(2)) = 1;
figure; imshow(sFImage);
imwrite(sFImage, 'starField.tif');


%% loop through several depths

%depthList = [4 5 6 7 8 9 10 12 14 16 24];
% depthList = [4  7 8 9 11 14 24];
depthList = [9 24];
%depthList = 4;

for i = depthList  
    reRenderScenes = true;
    inputFile = ['testTargetCA/realisticStarField' int2str(i) '.pbrt']

    %%% Scene rendering
    %%% Ambient only scene generation
    sceneName = ['Starfield Image; Depth: ' int2str(i)];
    oi = s3dRenderOI(inputFile, .050);
    oi = oiSet(oi, 'name', sceneName);
    % strip the file name from the path and assign that as the name of the
    % object  ... vcSaveObject(oi,);
    vcAddAndSelectObject(oi);
    oiWindow;
end