%% An algorithm that calculates depth using 2 flash-only images.  
% Andy Lin
% 
% The algorithm, inspired by Hany Farid, takes the ratio of the 2 images,
% and solves for the depth, using ray geometry, providing a crude depth.  
% Next, normal vectors are calculated using this crude depth map, and used
% to correct for the effects of lambertian surface reflectance, 
% providing a more accurate depth map.  
%
% **This script requires renderTwoFlashDepthScenePbrtObject.m to be run
% first.**  Either that, or the proper images loaded in


%% load 1st image - uncomment this section if you wish to load saved vcimages
%fullName = '2FlashDepth/indObject/idealDownFrontFlashImage.mat';  
%fullName = 'twoFlashDepth/depthTargetDepths/50mmFront.pbrt.image.mat';  
% fullName = 'twoFlashDepth/depthTargetDepths/50mmFront10s.pbrt.image.mat';  
name = 'frontFlashIp.mat';  

%load([s3dRootPath '/data/' fullName],'vci');
load([s3dRootPath '/papers/ReflectanceAndDepth/Data/03192015_depthEstimation/' name]);
vciFlash = vci;
vcAddObject(vciFlash);
ipWindow;

%% load 2nd flash image (flash now placed in back)
%fullName = '2FlashDepth/indObject/idealDownBackFlashImage.mat';
%fullName = 'twoFlashDepth/depthTargetDepths/50mmBack.pbrt.image.mat';  
% fullName = 'twoFlashDepth/depthTargetDepths/50mmBack10s.pbrt.image.mat';  
name = 'backFlashIp.mat';

% 
% multiplicationFactor = 1;  %to account for differences in exposure
% load([s3dRootPath '/data/' fullName],'vci');
load([s3dRootPath '/papers/ReflectanceAndDepth/Data/03192015_depthEstimation/' name]);
vciFlashBack = vci;
vcAddObject(vciFlashBack);
ipWindow;

%% load depth map - uncomment if you wish to load a depth map from file
%GTDepthMap = '2FlashDepth/indObject/depthMaps/groundTruthDepthMapDown.mat';
% GTDepthMap = 'twoFlashDepth/depthTargetDepths/GTDepthMap.mat';
% import = load(GTDepthMap);
% groundTruthDepthMap = import.depthMap; %ProcessedMedian;

%% Assign camera parameters

%if we are running this script right after
%renderTwoFlashDepthScenePbrtObject.m, it will look at the camera
%parameters from there - if not, assign these values.

if(~exist('curPbrt', 'var'))
    sensorWidth = 30.59;  %36; (used for indObject)
    sensorHeight  = 30.59;  %36; %24;
    
    
    sensorDistance = 140; %get this from the render file
    %f = 50;  % f signifies distance between 2 flashes    
else    
    sensorDistance = curPbrt.camera.lens.filmDistance;
    %height has a ratio of 1
    aspectRatio = curPbrt.camera.film.yresolution/curPbrt.camera.film.xresolution;
    hypotneuseRatio = sqrt(1^2 + aspectRatio^2);
    sensorHeight = (curPbrt.camera.lens.filmDiag/hypotneuseRatio);
    sensorWidth = (curPbrt.camera.lens.filmDiag/hypotneuseRatio) * aspectRatio;
    f = norm( curPbrt.camera.position(1,:) - curPbrt.lightSourceArray{1}.from, 2); %distance between 2 flashes
end

f = 6;  %make sure to set this
fieldOfView = atan(sensorWidth/2/sensorDistance) * 2 * 180/pi;  %this does not need to be changed

%% Process images and obtain ratio image
% The ratio image consists of the image with the front flash divided by the
% image with the back flash.  
flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
flashImageBack = imadjust(imageGet(vciFlashBack, 'results'), [], [], imageGet(vciFlash, 'gamma'));
multiplicationFactor = 1;  %to account for differences in exposure

% calculate ratio image
HSVFlashImage = rgb2hsv(flashImage);
HSVFlashImageBack = rgb2hsv(flashImageBack);

temp = imageGet(vciFlash, 'results');
linearIntensityFlash = sum(temp, 3);
linearIntensityFlash = linearIntensityFlash * multiplicationFactor; %for exposure adjustment!
temp = imageGet(vciFlashBack, 'results');
linearIntensityFlashBack = sum(temp, 3); 

%back flash 2 stuff
% temp = imageGet(vciFlashBack2, 'results');
% linearIntensityFlashBack2 = sum(temp, 3); 
%
% linearIntensityFlashBack = max(linearIntensityFlashBack, linearIntensityFlashBack2);

ratioImage = linearIntensityFlash./linearIntensityFlashBack;



%% Use ratio image to calculate crude depth
% After obtaining the ratio image, we can then go through the algebraic
% steps needed to calculate depth

% figure; imagesc(ratioImage);

% ------- stuff used for an ideal lens ---------
% we are using the magnification ratio to calculate field of view since
% this is a macro photography case.  See
% http://en.wikipedia.org/wiki/Angle_of_view for details.  
% alpha = 2 * arctan (d/(2 * F * (1 + m/P))
% here P = 1 since we are using an ideal lens, and m = S2/S1 = 133.33/80
%TODO: allow this to be calculated automatically
% fieldOfView = 15.4150; % this FOV takes into account magnification in a macro sense (does this work better?)

% fieldOfView = 25; % used for front back 100 experiment
% fieldOfView = 26.9915; % used for front back 100 experiment  %vertical
% fieldOfView = 39.5978;   %horizontal field of view
% ------- stuff used for an ideal lens ---------

% calculate alpha and phi for each pixel
xMatrix = linspace(-sensorWidth/2, sensorWidth/2, size(ratioImage,2));
xMatrix = repmat(xMatrix, [size(ratioImage,1), 1]);
yMatrix = linspace(-sensorHeight/2, sensorHeight/2, size(ratioImage,1))';
yMatrix = repmat(yMatrix, [1 size(ratioImage,2)]);
z = sensorDistance;
fakeD1 = sqrt(xMatrix.^2 + yMatrix.^2 + z.^2);

alpha = asin(xMatrix./fakeD1); 
phi   = asin(z./(fakeD1.*cos(alpha)));
% alpha = asin(xMatrix./fakeD1);
% phi   = atan(yMatrix./sqrt(yMatrix.^2 + z.^2));

%front back flash case
radical = abs(sqrt(4*cos(alpha).^2.*sin(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
d1Test = (2.*f.^2)./(-2.*cos(alpha).*sin(phi).*f + radical);
d1Test1st = d1Test;

figure; imagesc(d1Test);
colorbar; title('Calculated Depth (1st pass)'); %caxis([80 150]);

%% get quadratic coefficients out
c = (1 - ratioImage);
b = 2 .* cos(alpha) .* sin(phi) .*f;
a = ones(size(b)) * (f.^2);

%solves for 1/d1Test

%% general case for when flash in the back can be anywhere

fx = f/2;
fy = 0; %f/2;
fz = f/2 * sqrt(3);

b = -2 .* (fx .* cos(alpha) .* cos(phi) + fy .* sin(alpha) - fz .* cos(alpha) .* sin(phi));
a = fx^2 + fy^2 + fz^2;
c = 1 - ratioImage;

d1Test =1./((-b + abs(sqrt(b.^2 - 4 .*a .*c))) ./ (2.*a));


%% First filter the depth map using a separable median, and bilateral filter
%This will provide better data for calculating the normal map later
pixelUnit = sensorWidth / size(ratioImage,2);

%filter the depth map for better noise characteristics
d1TestMedFiltered = medianFilter(d1Test,5);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',5)';
d1TestFiltered = bilateralFilter(d1TestMedFiltered, 10, 4, 25);  
figure; imagesc(d1TestFiltered);
d1TestFiltered1st = d1TestFiltered;

colorbar; title('Filtered Depth'); caxis([80 150]);

%% Estimate the surface normals using the crude depth map
% We perform surface normal estimation by taking the cross product of 2
% perpendicular vectors on the surface.  This vector is averaged amongst
% the 4 pairs of vectors around 1 point of the surface.

% test case - use GT depth map to calculate depth map
% this allows us to decouple the effect of a bad depth map calculation
% method with the 2flashdepth algorithm

% d1TestFiltered = imresize(groundTruthDepthMap, [400 600]);

%calculate the x,y,z, point cloud from the depth map, then use this point
%cloud to determine the normal vectors
width = size(d1Test, 2);
height = size(d1Test, 1);
theta = linspace(-fieldOfView/2,fieldOfView/2, width);
phi = linspace(-fieldOfView/2,fieldOfView/2, height);

[thetaImage phiImage] = meshgrid(theta, phi);

yImage = sin(phiImage * pi/180) .* d1TestFiltered;
hyp = cos(phiImage * pi/180) .* d1TestFiltered;
xImage = hyp .* sin(thetaImage * pi/180);

zImage = hyp .* cos(thetaImage * pi/180);
% 
% figure; imagesc(zImage);
% figure; imagesc(xImage);
% figure; imagesc(yImage);

% calculate normal vectors using averaged cross products
recordedLateralDistance = zeros(size(d1Test));
normalMap = zeros([size(d1Test,1) size(d1Test,2) 3]);
for i = 2:(size(d1TestFiltered, 2) - 1)
    for j = 2:(size(d1TestFiltered,1) -1)
        
        
        aVectorX = xImage(j-1,i) - xImage(j,i) ;
        aVectorY = yImage(j-1,i) - yImage(j,i) ;
        aVectorZ = zImage(j-1,i) - zImage(j, i); 
        aVector = [aVectorX -aVectorY -aVectorZ];
        
        bVectorX = xImage(j,i+1) - xImage(j,i) ;
        bVectorY = yImage(j,i+1) - yImage(j,i) ;
        bVectorZ = zImage(j,i+1) - zImage(j, i); 
        bVector = [bVectorX -bVectorY -bVectorZ];       
        
        cVectorX = xImage(j+1,i) - xImage(j,i) ;
        cVectorY = yImage(j+1,i) - yImage(j,i) ;
        cVectorZ = zImage(j+1,i) - zImage(j, i); 
        cVector = [cVectorX -cVectorY -cVectorZ];
        
        
        dVectorX = xImage(j,i - 1) - xImage(j,i) ;
        dVectorY = yImage(j,i - 1) - yImage(j,i) ;
        dVectorZ = zImage(j,i - 1) - zImage(j, i); 
        dVector = [dVectorX -dVectorY -dVectorZ];
        
        
        adNormal = cross(aVector, dVector);
        adNormal = adNormal./norm(adNormal);
        
        dcNormal = cross(dVector, cVector);
        dcNormal = dcNormal./norm(dcNormal);
        
        cbNormal = cross(cVector, bVector);
        cbNormal = cbNormal./norm(cbNormal);
        
        baNormal = cross(bVector, aVector);
        baNormal = baNormal./norm(baNormal);
        
        %averageNormal = (adNormal + dcNormal + cbNormal + baNormal);
        averageNormal = median([adNormal',dcNormal',cbNormal',baNormal'],  2);
        averageNormal = averageNormal./norm(averageNormal);
        
        normalMap(j, i,:) = reshape(averageNormal, [1 1 3]);
    end
end
scaledNormalMap = normalMap./2 + .5;



figure; imshow(scaledNormalMap);
title('Calculated Normal Map');


%% use a "prior" as the normal map instead 

%in this case we will use a completely flat surface, with normals facing us
%(0,0,1)
zeroMap = zeros(size(d1Test));
oneMap = ones(size(d1Test));
normalMap = cat(3, zeroMap, zeroMap, oneMap);

%% filter the normal map...
% 
% %not sure what to use here since the 3 dimensions are coupled, but for now,
% %use a bilateral filter?

%filter first 2 dimensions
normalMap1 = medianFilter(normalMap(:,:,1),5);
normalMap1 = medianFilter(normalMap1',5)';
%normalMap1 = normalMap(:,:,1);
%normalMap1 = crossBilateralFilter(normalMap1, d1TestFiltered, 10, 4, 25);  
figure; imagesc(normalMap1);

normalMap2 = medianFilter(normalMap(:,:,2),5);
normalMap2 = medianFilter(normalMap2',5)';
%normalMap2 = normalMap(:,:,2);
%normalMap2 = crossBilateralFilter(normalMap2, d1TestFiltered, 10, 4, 25);  
figure; imagesc(normalMap2);

normalMap3 = ones(size(normalMap1)) - (normalMap1.^2 + normalMap2.^2);
normalMap3 = sqrt(normalMap3);

filtNormalMap = cat(3, normalMap1, normalMap2, normalMap3);

scaledFiltNormalMap = filtNormalMap./2 + .5;
figure; imshow(scaledFiltNormalMap);

%% Correcting for Lambert's Law error
% Depending on the relation of the surface normal, and the light direction,
% the intensity is attenuated according to the dot product (proportional to
% the cosine of the angle between the 2 vectors).  This effect was not
% taken into account during the initial calculation. 
linearIntensityFlashCorrected = linearIntensityFlash;
linearIntensityFlashBCorrected = linearIntensityFlashBack;

rayVectors1 = zeros([size(d1Test,1) size(d1Test,2) 3]);
rayVectors2 = zeros([size(d1Test,1) size(d1Test,2) 3]);

frontDot =  zeros([size(d1Test,1) size(d1Test,2)]);
backDot =  zeros([size(d1Test,1) size(d1Test,2)]);

%assuming field of view given above
fakeDistance = sensorWidth/2 / tan(fieldOfView/2 * pi/180);

numWidth = size(d1Test,2);
numHeight = size(d1Test,1);

for i = 1:(size(d1TestFiltered, 2))
    for j = 1:(size(d1TestFiltered,1))      
        % calculate normal vector from front light
        fakeX = pixelUnit * (i - numWidth/2);
        fakeY = pixelUnit * (j - numHeight/2);
        tempVector = [-fakeX fakeY fakeDistance];
        tempVector = tempVector ./ norm(tempVector);
        rayVectors1(j, i, :) = reshape(tempVector, [1 1 3]);

        %use similar triangles to recalculate rayVectors2 (normal vectors
        %for back light
        fakeH = sqrt(fakeX^2 + fakeY^2 + fakeDistance^2);
        fake2RealRatio = (d1TestFiltered(j,i))/fakeH;
        realX = fake2RealRatio* fakeX;
        realY = fake2RealRatio * fakeY;
        realDistance = fake2RealRatio * fakeDistance;
        
        %for front back
        tempVector = [-realX realY realDistance + f];
        
        tempVector = tempVector ./ norm(tempVector);
        rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);
        frontDot(j,i) = sum(rayVectors1(j,i,:) .* filtNormalMap(j,i,:));
        backDot(j,i) = sum(rayVectors2(j,i,:) .* filtNormalMap(j,i,:));
        
        %regularize the dot products - we will never allow the dot product
        %to be 0, so that the value doesn't blow up.
        %frontDot(j,i) = .2 * exp(-2.*frontDot(j,i)) + frontDot(j,i);
        %backDot(j,i) = .2 * exp(-2.*backDot(j,i)) + backDot(j,i);
        
        %calculate correction factors for both images
        linearIntensityFlashCorrected(j,i, :) = linearIntensityFlashCorrected(j,i, :) ./ frontDot(j,i) ; 
        linearIntensityFlashBCorrected(j,i, :) = linearIntensityFlashBCorrected(j,i, :) ./ backDot(j,i) ; 
    end
end 

ratioImage = linearIntensityFlashCorrected./linearIntensityFlashBCorrected;


%% Use ratio image to calculate distance (2nd pass)
% We've corrected for Lambert's Law, so now we can compute the calculated
% distance once again using the identical math.

xMatrix = linspace(-sensorWidth/2, sensorWidth/2, size(ratioImage,2));
xMatrix = repmat(xMatrix, [size(ratioImage,1), 1]);
yMatrix = linspace(-sensorHeight/2, sensorHeight/2, size(ratioImage,1))';
yMatrix = repmat(yMatrix, [1 size(ratioImage,2)]);
z = sensorDistance;

fakeD1 = sqrt(xMatrix.^2 + yMatrix.^2 + z.^2);
alpha = asin(xMatrix./fakeD1); 
phi   = asin(z./(fakeD1.*cos(alpha)));

% alpha = atan(xMatrix./z);
% phi   = atan(yMatrix./sqrt(yMatrix.^2 + z.^2));

%front back flash case
radical = abs(sqrt(4*cos(alpha).^2.*sin(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
d1Test = (2.*f.^2)./(-2.*cos(alpha).*sin(phi).*f + radical);

figure; imagesc(d1Test);
colorbar; title('Calculated Depth (2nd pass)'); %caxis([80 150]);

%% get quadratic coefficients out
ratioImageCor = ratioImage;
ratioImageCor(ratioImageCor<1) = 1;
c = (1 - ratioImageCor);
b = 2 .* cos(alpha) .* sin(phi) .*f;
a = ones(size(b)) * (f.^2);


%% Filter the depth map using a separable median, and bilateral filter (2nd pass)
%This will provide better data 

%filter the depth map for better noise characteristics
d1TestMedFiltered = medianFilter(d1Test,5);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',5)';
d1TestFiltered = bilateralFilter(d1TestMedFiltered, 10, 4, 25);  
figure; imagesc(d1TestFiltered);
colorbar; title('Filtered Depth (2nd pass)'); caxis([80 150]);

%% Compare the depth map to the ground truth
figure; imagesc(groundTruthDepthMap);
colorbar; title('Ground Truth Depth Map'); caxis([80 150]);

%% Metric for error

%resize images to same size as original
nativeSize = size(groundTruthDepthMap);
d1TestResized = imresize(d1Test1st, nativeSize);
% d1TestFilteredResized = imresize(d1TestFiltered, nativeSize);

error = abs(groundTruthDepthMap - d1TestResized)./groundTruthDepthMap;
realValues = ~isnan(error);
meanPercentError = mean(error(realValues)) * 100


d1TestFilteredResized = imresize(d1TestFiltered1st, nativeSize);
errorFiltered = abs(groundTruthDepthMap - d1TestFilteredResized)./groundTruthDepthMap;
realValues = ~isnan(errorFiltered);
meanPercentErrorFiltered= mean(errorFiltered(realValues)) * 100

%% Summary
% We were able to calculate a reasonable depth map from the 2-flash
% algorithm.  Benefits of this algorithm compared to traditional 2-camera
% stereo algorithms include better accuracy for uniform (non-textured)
% regions and superior resolution.  We have yet to make a formal comparison
% of these 2 techniques. The depths may still be off by a scaling factor
% and an offset, but this initial study shows that the algorithm works
% reasonably well.  

%% Future work
% We plan to show the true benefits of this method of calculating a depth
% map.  The algorithm works identically for cases where the flashes, and
% cameras are in different positions.  We plan to investigate this
% scenario.  Also, we have hypothesized that this algorithm will work, with
% the flashes at different positions (for example, side-by-side), but the
% math involved will be slightly different.  We plan to investigate these
% scenarios as well.  Finally, we also plan to form experiments where we
% use local windows to improve accuracy, while sacrificing resolution.  



% OLD NorMAL VECTOR CALCULATION

% % calculate normal vectors using averaged cross products
% recordedLateralDistance = zeros(size(d1Test));
% normalMap = zeros([size(d1Test,1) size(d1Test,2) 3]);
% for i = 2:(size(d1TestFiltered, 2) - 1)
%     for j = 2:(size(d1TestFiltered,1) -1)
%         aRelief = d1TestFiltered(j - 1,i) - d1TestFiltered(j,i);
%         bRelief = d1TestFiltered(j, i + 1) - d1TestFiltered(j,i);
%         cRelief = d1TestFiltered(j + 1,i) - d1TestFiltered(j,i);
%         dRelief = d1TestFiltered(j,i - 1) - d1TestFiltered(j,i);
%         
%         lateralDistance = d1TestFiltered(j,i) * pixelUnit/sensorDistance;
%         %%this does not make a lot of sense
%         %lateralDistance = d1TestFiltered(j,i) * tan(fieldOfView/size(ratioImage,2) * pi/180);
%         %lateralDistance = 100 * tan(fieldOfView/size(ratioImage,2) * pi/180);
%         recordedLateralDistance(j,i) = lateralDistance;
%         aVector = [0 lateralDistance -aRelief];
%         bVector = [lateralDistance 0 -bRelief];
%         cVector = [0 -lateralDistance -cRelief];
%         dVector = [-lateralDistance 0 -dRelief];
%         
%         adNormal = cross(aVector, dVector);
%         adNormal = adNormal./norm(adNormal);
%         
%         dcNormal = cross(dVector, cVector);
%         dcNormal = dcNormal./norm(dcNormal);
%         
%         cbNormal = cross(cVector, bVector);
%         cbNormal = cbNormal./norm(cbNormal);
%         
%         baNormal = cross(bVector, aVector);
%         baNormal = baNormal./norm(baNormal);
%         
%         %averageNormal = (adNormal + dcNormal + cbNormal + baNormal);
%         averageNormal = median([adNormal',dcNormal',cbNormal',baNormal'],  2);
%         averageNormal = averageNormal./norm(averageNormal);
%         
%         normalMap(j, i,:) = reshape(averageNormal, [1 1 3]);
%     end
% end
% scaledNormalMap = normalMap./2 + .5;
