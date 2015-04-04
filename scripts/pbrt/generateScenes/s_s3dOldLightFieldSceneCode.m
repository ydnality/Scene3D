%% render light field


numPinholesW = 80;  %these 2 parameters must be even (for now)
numPinholesH = 80;
superPixelPitch = filmDiag/sqrt(2)/numPinholesW;  %assumes square aspect ratio
pinholeArrayDist = .3981;  %calculated using equation

%todo: don't hard code 9,9, 31
lightField = zeros(numPinholesW, numPinholesH, 9, 9, 31);

for i = 1:numPinholesW
    for j = 1:numPinholesH
        i
        j
        % identify center position of superPixel
        centerPos = [(i - numPinholesW/2 - .5) * superPixelPitch (j- numPinholesH/2 - .5) * superPixelPitch filmDist];
        
        % trace chief ray from center position to center of exit aperture
        chiefRayDir = [ 0 0 0 ] - centerPos;
        chiefRayDir = chiefRayDir/norm(chiefRayDir, 2); 
        
        t = (filmDist - pinholeArrayDist)/chiefRayDir(3);
        pinholePos = [t * chiefRayDir(1) t * chiefRayDir(2)];
        
        pinholeExitApLoc = [pinholePos -(filmDist - pinholeArrayDist)];  %note that sensor is negative in pbrt world
        filmDiagSmall = filmDiag/numPinholesW;  
        %assign pinhole position to PBRT, and figure out correct cropWindow
        lens = pbrtLensRealisticObject(filmDist, filmDiagSmall, specFile, apertureDiameter, ...
            diffraction, chromaticAberration, [], pinholeExitApLoc, [centerPos(1) centerPos(2)]);
        curPbrt.camera.setLens(lens);

        
        
        % we need to reset the sampler because of some cloning limitations   
        % Sampler
        sampler = pbrtSamplerObject();
        samples = curPbrt.sampler.removeProperty();
        samples.value = 128;
        curPbrt.setSampler(sampler);

        %curPbrt.camera.setResolution(720, 720);
        curPbrt.camera.setResolution(9, 9);
        
        %widthUnit = 1/numPinholesW;
        %heightUnit = 1/numPinholesH;
        
             
        %order for cropWindow is x y
        %curPbrt.camera.setCropWindow(0, 1, 0, 1);
        %curPbrt.camera.setCropWindow((numPinholesW - i)/numPinholesW, (numPinholesW - i + 1)/numPinholesW, (j-1)/numPinholesH,  (j)/numPinholesH);
        
        %uncomment to use a 2 element lens instead of a pinhole
        % curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));

        oi = s3dRenderSceneAndDepthMap(curPbrt, 'simpleScene', true);
        %vcAddObject(scene); sceneWindow;
        
        fullName = vcSaveObject(oi, fullfile(dataPath, 'pbrtScenes', 'benchScene', 'LF', ['superpixel' int2str(i) '_' int2str(j) '.mat']));
        photons = sceneGet(oi, 'photons');
        lightField(i,j, :,:, :) = photons(1:9, 1:9, :);  %there is an annoying rounding bug in pbrt.  This could be a problem.  hopefully not... 
    end
end
%instead of visualizing scenes, we will save it
%fullName = vcSaveObject(scene, fullfile(dataPath, 'pbrtScenes', 'benchScene', 'HDRVideo', ['frame' int2str(i) '.mat']));



%% try rendering some light field sub images at specific aperture locations.  
% This renders a pinhole image at a particular aperture location. 

testScene = sceneCreate;
centerImagePhotons = lightField(:,:, 7, 7, :);  %7,7 is a specific aperture location.  The middle would be 5,5.
centerImagePhotons = reshape(centerImagePhotons, [80 80 31]);
centerImagePhotons = permute(centerImagePhotons, [2 1 3]);
centerImagePhotons = centerImagePhotons(:,end:-1:1, :);
testScene = sceneSet(testScene, 'photons', centerImagePhotons);
vcAddObject(testScene); sceneWindow;

%% try rendering a basic light field blurred image

%sum all the sub aperture views (3rd and 4th dimensions)  
summedimage = sum(sum(lightField, 3), 4);
summedimage = reshape(summedimage, [80 80 31]);
summedimage = permute(summedimage, [2 1 3]);
summedimage = summedimage(:,end:-1:1, :);

testScene = sceneSet(testScene, 'photons', summedimage);
vcAddObject(testScene); sceneWindow;



 
%% Load scenes from file and store it as a light field
numPinholesW = 80;  %these 2 parameters must be even (for now)
numPinholesH = 80;

%todo: don't hard code 9,9, 31
lightField = zeros(numPinholesW, numPinholesH, 9, 9, 31);

for i = 1:numPinholesW
    i
    for j = 1:numPinholesH
        loadedScene = load(fullfile(dataPath, 'pbrtScenes', 'benchScene', 'LF', ['superpixel' int2str(i) '_' int2str(j) '.mat']));
        photons = sceneGet(loadedScene.scene, 'photons');
        %vcAddObject(loadedScene.scene); sceneWindow;
        lightField(i,j, :,:, :) = photons(1:9, 1:9, :); 
    end
end
%%
