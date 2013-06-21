%load flash image
% fullName = 'benchFlash.mat';
fullName = 'indObjFlash.mat';
% fullName = 'indObjFlashLessStrong.mat';
% fullName = 'deskFlash.mat';
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

%load no-flash image
% fullName = 'benchNoFlash.mat';
% fullName = 'indObjNoFlash.mat';
fullName = 'indObjNoFlashLessNoise.mat';
% fullName = 'indObjFlashLessStrong.mat';
% fullName = 'deskNoFlash.mat';
%fullName = 'deskTest.mat';   %testing exposure things
load(fullName,'vci');
vciNoFlash = vci;
vcAddAndSelectObject('vcimage',vciNoFlash);
vcimageWindow;


%% obtain the ratio image

flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
HSVFlashImage = rgb2hsv(flashImage);
intensityImage = HSVFlashImage(:,:,3);

HSVNoFlashImage = rgb2hsv(noFlashImage);
intensityImageNF = HSVNoFlashImage(:,:,3);

ratioImage = intensityImage./intensityImageNF;
figure; imagesc(ratioImage);

%% obtain flash only image - use this to help with finding depth

flashExposure = 1.5020e-35;
noFlashExposure = 0.0138;

%flash only image can help build reflectance image - problem is that depth
%should be used because flash intensity diminishes with distance - so it is
%better to know depth in order to better estimate the reflectance.



%% use 2 flash distances to estimate distance
% fullName = 'indObjFlashLessStrong.mat';       %autoExpTime =  0.2230
% fullName = 'indObjFlashLambertianNoAmbient.mat'; 
% fullName = 'indObject2FlashDepth/frontFlashImageLambertian.mat'; 
% fullName = 'indObject2FlashDepth/frontFlashImageLambertianGT.mat'; 
% fullName = 'indObject2FlashDepth/frontFlashImageLambertian0.mat'; 
% fullName = 'indObject2FlashDepth/grayFrontImage.mat'; 
fullName = 'indObject2FlashDepth/newFrontFlashImage.mat'; 
fullName = 'indObject2FlashDepth/downFrontFlashImage.mat'; 
% fullName = 'floorWallBottomBack/frontFlashDownImage.mat'; 
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

% fullName = 'indObjFlashBackLambertianNoAmbient.mat';     %autoExpTime =0.2361
% fullName = 'indObjFlashLessStrongBackLambertian.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertianCloser.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertianGT.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertian0.mat'; 
% fullName = 'indObject2FlashDepth/grayBackImage.mat'; 

fullName = 'indObject2FlashDepth/newBackFlashImage.mat'; 
fullName = 'indObject2FlashDepth/downBackFlashImage.mat'; 
% fullName = 'floorWallBottomBack/backFlashDownImage.mat'; 
% fullName = 'floorWallBottomBack/backFlashDown100Image.mat'; 
% fullName = 'floorWallBottomBack/sideFlashDownImage.mat'; 
% fullName = 'floorWallBottomBack/sideFlashDown25Image.mat'; 
multiplicationFactor = 16/8; %16 for the change in exposure, 8 for change in samples
multiplicationFactor = 8/8;
multiplicationFactor = 1;
load(fullName,'vci');
vciFlashBack = vci;
vcAddAndSelectObject('vcimage',vciFlashBack);
vcimageWindow;

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


f = 100; %50; %5;
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



%% estimate the surface normals using a depth map

%ignore boundary cases for now - add back in later if wanted

pixelUnit = sensorWidth / size(ratioImage,2);
% d1TestFiltered = imfilter(d1Test,fspecial( 'gaussian', 10, 10));
% d1TestFiltered = d1Test;
% d1TestFiltered = bilateralFilter(d1Test, 40, 20, 30);  %filtering operation for better depth map

d1TestMedFiltered = medianFilter(d1Test,5);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',5)';
d1TestFiltered = bilateralFilter(d1TestMedFiltered, 10, 4, 25);  

%**depthMapProcessedMedian is the ground truth depth map
% d1TestFiltered = bilateralFilter(depthMapProcessedMedian, 3, 4, 9);  
% % d1TestFiltered = depthMapProcessedMedian;   % using ground truth depth map!
% d1TestFiltered = imresize(d1TestFiltered, size(d1Test));   % using ground truth depth map!


% d1TestMedFiltered = medianFilter(d1Test,15);
% d1TestMedFiltered = medianFilter(d1TestMedFiltered',15)';
% d1TestFiltered = d1TestMedFiltered;

% figure; imagesc(d1TestMedFiltered);
figure; imagesc(d1TestFiltered);
% d1TestFiltered = d1TestMedFiltered;

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
% figure; imshow(filteredNormals);

% test = sum(normalMap .* normalMap, 3); %normalization testing

%% correcting for lamberts law error

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
        
        fakeX = pixelUnit * (i - numWidth/2);
        fakeY = pixelUnit * (j - numHeight/2);
        tempVector = [-fakeX fakeY fakeDistance];
        tempVector = tempVector ./ norm(tempVector);
        rayVectors1(j, i, :) = reshape(tempVector, [1 1 3]);
        
%         tempVector = [- .5 * pixelUnit * (i - numWidth/2) .5 * pixelUnit * (j - numHeight/2) (.5 * fakeDistance + f)];
%         tempVector = tempVector ./ norm(tempVector);        
%         rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);
        
        
        %new vector2 estimation technique involving depth information
%         fakeX = pixelUnit * (i - numWidth/2);
%         fakeY = -pixelUnit * (j - numHeight/2);
%         
%         
%         alpha =  asin(fakeY/d1Test(j,i));
%         phi = asin(fakeDistance/(d1Test(j,i)*cos(alpha)));
%         
%         px = d1Test(j,i) * cos(alpha) * cos(phi);
%         py = d1Test(j,i) * sin(alpha);
%         pz = -d1Test(j,i) * cos(alpha) * sin(phi);
%         
%         tempVector = [-px -py -pz];      
%         tempVector = tempVector ./ norm(tempVector);
%         rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);


        %use similar triangles to recalculate rayVectors2
        fakeH = sqrt(fakeX^2 + fakeY^2 + fakeDistance^2);
%         fake2RealRatio = (d1Test(j,i) - 40)/fakeH;
        fake2RealRatio = (d1TestFiltered(j,i))/fakeH;
        realX = fake2RealRatio* fakeX;
        realY = fake2RealRatio * fakeY;
        realDistance = fake2RealRatio * fakeDistance;
        
        %for front back
        tempVector = [-realX realY realDistance + f];
        
        %for side by side
%         tempVector = [-(realX - f) realY realDistance];
        
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
figure; imagesc(ratioImage);




%% calculate depth map error - perhaps this will help us figure out what is wrong with the algorithm

errorMap = depthMapProcessedMedian - imresize(d1Test, [300 450]);
figure; imagesc(errorMap);



