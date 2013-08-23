%% An algorithm that calculates depth using 2 flash-only images.  
% Andy Lin
% 
% The algorithm, inspired by Hany Farid, takes the ratio of the 2 images,
% and solves for the depth, using ray geometry, providing a crude depth.  
% Next, normal vectors are calculated using this crude depth map, and used
% to correct for the effects of lambertian surface reflectance, 
% providing a more accurate depth map.  

%% A simple test to estimate depth using flash/no-flash
% % Load images
% 
% % load flash image
% % fullName = 'benchFlash.mat';
% fullName = 'indObjFlash.mat';
% % fullName = 'indObjFlashLessStrong.mat';
% % fullName = 'deskFlash.mat';
% 
% load(fullName,'vci');
% vciFlash = vci;
% vcAddAndSelectObject('vcimage',vciFlash);
% vcimageWindow;
% 
% % load no-flash image
% % fullName = 'benchNoFlash.mat';
% % fullName = 'indObjNoFlash.mat';
% fullName = 'indObjNoFlashLessNoise.mat';
% % fullName = 'indObjFlashLessStrong.mat';
% % fullName = 'deskNoFlash.mat';
% %fullName = 'deskTest.mat';   %testing exposure things
% 
% load(fullName,'vci');
% vciNoFlash = vci;
% vcAddAndSelectObject('vcimage',vciNoFlash);
% vcimageWindow;
% 
% % Obtain the ratio image
% 
% flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
% noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
% HSVFlashImage = rgb2hsv(flashImage);
% intensityImage = HSVFlashImage(:,:,3);
% 
% HSVNoFlashImage = rgb2hsv(noFlashImage);
% intensityImageNF = HSVNoFlashImage(:,:,3);
% 
% ratioImage = intensityImage./intensityImageNF;
% figure; imagesc(ratioImage);
% 
% % Obtain flash only image - use this to help with finding depth
% 
% flashExposure = 1.5020e-35;
% noFlashExposure = 0.0138;
% 
% %flash only image can help build reflectance image - problem is that depth
% %should be used because flash intensity diminishes with distance - so it is
% %better to know depth in order to better estimate the reflectance.
% 
% %this experiment was meant to show that flash/no-flash images are
% %inadequate in providing reliable depth from images.  

%% Use 2 flashes (one in back of another) to estimate distance
%% load images

% fullName = 'indObjFlashLessStrong.mat';       %autoExpTime =  0.2230
% fullName = 'indObjFlashLambertianNoAmbient.mat'; 
% fullName = '2FlashDepth/indObject/frontFlashImageLambertian.mat'; 
% fullName = '2FlashDepth/indObject/frontFlashImageLambertianGT.mat'; 
% fullName = '2FlashDepth/indObject/frontFlashImageLambertian0.mat'; 
% fullName = '2FlashDepth/indObject/grayFrontImage.mat'; 

% fullName = '/2FlashDepth/indObject/newFrontFlashImage.mat'; 
fullName = '2FlashDepth/indObject/downFrontFlashImage.mat'; 
% fullName = '/2FlashDepth/indObject/downFrontFlashImage12bit.mat'; 

% fullName = 'floorWallBottomBack/frontFlashDownImage.mat'; 
load([s3dRootPath '/compPhotography/' fullName],'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

% fullName = 'indObjFlashBackLambertianNoAmbient.mat';     %autoExpTime =0.2361
% fullName = 'indObjFlashLessStrongBackLambertian.mat'; 
% fullName = '2FlashDepth/indObject/backFlashImageLambertianCloser.mat'; 
% fullName = '2FlashDepth/indObject/backFlashImageLambertianGT.mat'; 
% fullName = '2FlashDepth/indObject/backFlashImageLambertian0.mat'; 
% fullName = '2FlashDepth/indObject/grayBackImage.mat'; 
%% load 2nd flash image (flash now placed in back)
% fullName = '2FlashDepth/indObject/newBackFlashImage.mat'; 
fullName = '2FlashDepth/indObject/downBackFlashImage.mat'; 
% fullName = '2FlashDepth/indObject/downBackFlashImage12bit.mat'; 

% fullName = 'floorWallBottomBack/backFlashDownImage.mat'; 
% fullName = 'floorWallBottomBack/backFlashDown100Image.mat'; 
% fullName = 'floorWallBottomBack/sideFlashDownImage.mat'; 
% fullName = 'floorWallBottomBack/sideFlashDown25Image.mat'; 

% multiplicationFactor = 16/8; %16 for the change in exposure, 8 for change in samples
% multiplicationFactor = 8/8;
multiplicationFactor = 1;
load([s3dRootPath '/compPhotography/' fullName],'vci');
vciFlashBack = vci;
vcAddAndSelectObject('vcimage',vciFlashBack);
vcimageWindow;

%% Process images and obtain ratio image
% The ratio image consists of the image with the front flash divided by the
% image with the back flash.  
flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
flashImageBack = imadjust(imageGet(vciFlashBack, 'results'), [], [], imageGet(vciFlash, 'gamma'));
% figure; imshow(flashImage); figure; imshow(flashImageBack);

%ratio image
HSVFlashImage = rgb2hsv(flashImage);
HSVFlashImageBack = rgb2hsv(flashImageBack);

% linearIntensityFlash = imadjust(HSVFlashImage(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));
% linearIntensityFlashBack = imadjust(HSVFlashImageBack(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));

temp = imageGet(vciFlash, 'results');
linearIntensityFlash = sum(temp, 3);
linearIntensityFlash = linearIntensityFlash * multiplicationFactor; %for exposure adjustment!
temp = imageGet(vciFlashBack, 'results');
linearIntensityFlashBack = sum(temp, 3); 

% figure; imshow(linearIntensityFlash)
% figure; imshow(linearIntensityFlashBack)
ratioImage = linearIntensityFlash./linearIntensityFlashBack;

%% Use ratio image to calculate crude depth
% After obtaining the ratio image, we can then go through the algebraic
% steps needed to calculate depth

% figure; imagesc(ratioImage);
% fieldOfView = 39.60;
% fieldOfView = 20;
fieldOfView = 25; % used for front back 100 experiment
% fieldOfView = 30;  % used for side experiment

% new code for 2D
% alpha = asin(p(2)/d1); 
% phi   = asin(p(3)/(d1*cos(alpha)));

sensorWidth = 36;
sensorHeight  = 24;
% sensorDistance = 11.8395;

% sensorDistance = 49.9969;
sensorDistance = sensorWidth/2 / tan(fieldOfView/2 * pi/180);
% sensorDistance = 75;  %hmm this changes with changing focus...

% calculate alpha and phi for each pixel
xMatrix = linspace(-sensorWidth/2, sensorWidth/2, size(ratioImage,2));
xMatrix = repmat(xMatrix, [size(ratioImage,1), 1]);
yMatrix = linspace(-sensorHeight/2, sensorHeight/2, size(ratioImage,1))';
yMatrix = repmat(yMatrix, [1 size(ratioImage,2)]);
z = sensorDistance;
fakeD1 = sqrt(xMatrix.^2 + yMatrix.^2 + z.^2);
alpha = asin(xMatrix./fakeD1); 
phi   = asin(z./(fakeD1.*cos(alpha)));

% f = 100; %50; %5;
f = 50;
% f = 25;
%front back flash case
radical = abs(sqrt(4*cos(alpha).^2.*sin(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
d1Test = (2.*f.^2)./(-2.*cos(alpha).*sin(phi).*f + radical);

%side by side flash case
% radical = abs(sqrt(4*cos(alpha).^2.*cos(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
% % figure; imagesc(ratioImage > 1);
% % figure; imagesc(ratioImage)
% 
% d1Test = (2.*f.^2)./(2.*cos(alpha).*cos(phi).*f + radical);
figure; imagesc(d1Test);
colorbar; title('Calculated Depth (1st pass)'); caxis([80 150]);

%% First filter the depth map using a separable median, and bilateral filter
%This will provide better data for calculating the normal map later

pixelUnit = sensorWidth / size(ratioImage,2);
% d1TestFiltered = imfilter(d1Test,fspecial( 'gaussian', 10, 10));
% d1TestFiltered = d1Test;
% d1TestFiltered = bilateralFilter(d1Test, 40, 20, 30);  %filtering operation for better depth map

%filter the depth map for better noise characteristics
d1TestMedFiltered = medianFilter(d1Test,5);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',5)';
d1TestFiltered = bilateralFilter(d1TestMedFiltered, 10, 4, 25);  
figure; imagesc(d1TestFiltered);
colorbar; title('Filtered Depth'); caxis([80 150]);

%% Estimate the surface normals using the crude depth map
% We perform surface normal estimation by taking the cross product of 2
% perpendicular vectors on the surface.  This vector is averaged amongst
% the 4 pairs of vectors around 1 point of the surface.

%ignore boundary cases for now - add back in later if wanted
%**depthMapProcessedMedian is the ground truth depth map
% d1TestFiltered = bilateralFilter(depthMapProcessedMedian, 3, 4, 9);  
% % d1TestFiltered = depthMapProcessedMedian;   % using ground truth depth map!
% d1TestFiltered = imresize(d1TestFiltered, size(d1Test));   % using ground truth depth map!

% d1TestMedFiltered = medianFilter(d1Test,15);
% d1TestMedFiltered = medianFilter(d1TestMedFiltered',15)';
% d1TestFiltered = d1TestMedFiltered;

% figure; imagesc(d1TestMedFiltered);
% d1TestFiltered = d1TestMedFiltered;

% calculate normal vectors using averaged cross products
recordedLateralDistance = zeros(size(d1Test));
normalMap = zeros([size(d1Test,1) size(d1Test,2) 3]);
for i = 2:(size(d1TestFiltered, 2) - 1)
    for j = 2:(size(d1TestFiltered,1) -1)
        aRelief = d1TestFiltered(j - 1,i) - d1TestFiltered(j,i);
        bRelief = d1TestFiltered(j, i + 1) - d1TestFiltered(j,i);
        cRelief = d1TestFiltered(j + 1,i) - d1TestFiltered(j,i);
        dRelief = d1TestFiltered(j,i - 1) - d1TestFiltered(j,i);
        
        lateralDistance = d1TestFiltered(j,i) * pixelUnit/sensorDistance; 
        recordedLateralDistance(j,i) = lateralDistance;
        aVector = [0 lateralDistance -aRelief];
        bVector = [lateralDistance 0 -bRelief];
        cVector = [0 -lateralDistance -cRelief];
        dVector = [-lateralDistance 0 -dRelief];
        
        adNormal = cross(aVector, dVector);
        adNormal = adNormal./norm(adNormal);
        
        dcNormal = cross(dVector, cVector);
        dcNormal = dcNormal./norm(dcNormal);
        
        cbNormal = cross(cVector, bVector);
        cbNormal = cbNormal./norm(cbNormal);
        
        baNormal = cross(bVector, aVector);
        baNormal = baNormal./norm(baNormal);
        
        averageNormal = (adNormal + dcNormal + cbNormal + baNormal);
        averageNormal = averageNormal./norm(averageNormal);
        
        normalMap(j, i,:) = reshape(averageNormal, [1 1 3]);
    end
end

scaledNormalMap = normalMap./2 + .5;
% filteredNormals = imfilter(scaledNormalMap,fspecial( 'gaussian', 5, 5));

figure; imshow(scaledNormalMap);
title('Calculated Normal Map');
% figure; imshow(filteredNormals);

% test = sum(normalMap .* normalMap, 3); %normalization testing

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

        %old code
%         tempVector = [- .5 * pixelUnit * (i - numWidth/2) .5 * pixelUnit * (j - numHeight/2) (.5 * fakeDistance + f)];
%         tempVector = tempVector ./ norm(tempVector);        
%         rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);
        
        %new vector2 estimation technique involving depth information
%         fakeX = pixelUnit * (i - numWidth/2);
%         fakeY = -pixelUnit * (j - numHeight/2);
%         
%         alpha =  asin(fakeY/d1Test(j,i));
%         phi = asin(fakeDistance/(d1Test(j,i)*cos(alpha)));
%         px = d1Test(j,i) * cos(alpha) * cos(phi);
%         py = d1Test(j,i) * sin(alpha);
%         pz = -d1Test(j,i) * cos(alpha) * sin(phi);
%         
%         tempVector = [-px -py -pz];      
%         tempVector = tempVector ./ norm(tempVector);
%         rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);

        %use similar triangles to recalculate rayVectors2 (normal vectors
        %for back light
        fakeH = sqrt(fakeX^2 + fakeY^2 + fakeDistance^2);
%         fake2RealRatio = (d1Test(j,i) - 40)/fakeH;
        fake2RealRatio = (d1TestFiltered(j,i))/fakeH;
        realX = fake2RealRatio* fakeX;
        realY = fake2RealRatio * fakeY;
        realDistance = fake2RealRatio * fakeDistance;
        
        %for front back
        tempVector = [-realX realY realDistance + f];
        
        %for side by side
        % tempVector = [-(realX - f) realY realDistance];
        
        tempVector = tempVector ./ norm(tempVector);
        rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);
        
        frontDot(j,i) = sum(rayVectors1(j,i,:) .* normalMap(j,i,:));
        backDot(j,i) = sum(rayVectors2(j,i,:) .* normalMap(j,i,:));
        
        %calculate correction factors for both images
        linearIntensityFlashCorrected(j,i, :) = linearIntensityFlashCorrected(j,i, :) ./ frontDot(j,i) ; 
        linearIntensityFlashBCorrected(j,i, :) = linearIntensityFlashBCorrected(j,i, :) ./ backDot(j,i) ; 
    end
end 

% debug print
% figure; imshow(rayVectors ./2 + .5);

ratioImage = linearIntensityFlashCorrected./linearIntensityFlashBCorrected;
% figure; imagesc(ratioImage);


%% Use ratio image to calculate distance (2nd pass)
% We've corrected for Lambert's Law, so now we can compute the calculated
% distance once again using the identical math.

% figure; imagesc(ratioImage);
% fieldOfView = 39.60;
% fieldOfView = 20;
fieldOfView = 25; % used for front back 100 experiment
% fieldOfView = 30;  % used for side experiment

% new code for 2D

% alpha = asin(p(2)/d1); 
% phi   = asin(p(3)/(d1*cos(alpha)));

sensorWidth = 36;
sensorHeight  = 24;
% sensorDistance = 11.8395;

% sensorDistance = 49.9969;
sensorDistance = sensorWidth/2 / tan(fieldOfView/2 * pi/180);
% sensorDistance = 75;  %hmm this changes with changing focus...

xMatrix = linspace(-sensorWidth/2, sensorWidth/2, size(ratioImage,2));
xMatrix = repmat(xMatrix, [size(ratioImage,1), 1]);
yMatrix = linspace(-sensorHeight/2, sensorHeight/2, size(ratioImage,1))';
yMatrix = repmat(yMatrix, [1 size(ratioImage,2)]);
z = sensorDistance;

fakeD1 = sqrt(xMatrix.^2 + yMatrix.^2 + z.^2);
alpha = asin(xMatrix./fakeD1); 
phi   = asin(z./(fakeD1.*cos(alpha)));

% f = 100; %50; %5;
f = 50;
% f = 25;
%front back flash case
radical = abs(sqrt(4*cos(alpha).^2.*sin(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
d1Test = (2.*f.^2)./(-2.*cos(alpha).*sin(phi).*f + radical);

%side by side flash case
% radical = abs(sqrt(4*cos(alpha).^2.*cos(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
% % figure; imagesc(ratioImage > 1);
% % figure; imagesc(ratioImage)
% 
% d1Test = (2.*f.^2)./(2.*cos(alpha).*cos(phi).*f + radical);
figure; imagesc(d1Test);
colorbar; title('Calculated Depth (2nd pass)'); caxis([80 150]);

% %% calculate depth map error - perhaps this will help us figure out what is wrong with the algorithm
% 
% % errorMap = depthMapProcessedMedian - imresize(d1Test, [300 450]);
% % figure; imagesc(errorMap);

%% Filter the depth map using a separable median, and bilateral filter (2nd pass)
%This will provide better data 

%filter the depth map for better noise characteristics
d1TestMedFiltered = medianFilter(d1Test,5);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',5)';
d1TestFiltered = bilateralFilter(d1TestMedFiltered, 10, 4, 25);  
figure; imagesc(d1TestFiltered);
colorbar; title('Filtered Depth (2nd pass)'); caxis([80 150]);

%% Compare the depth map to the ground truth
import = load('2FlashDepth/indObject/depthMaps/groundTruthDepthMapDown.mat');
groundTruthDepthMap = import.depthMapProcessedMedian;
figure; imagesc(groundTruthDepthMap);
colorbar; title('Ground Truth Depth Map'); caxis([60 120]);
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


