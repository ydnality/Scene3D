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
fullName = 'indObject2FlashDepth/grayFrontImage.mat'; 
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

% fullName = 'indObjFlashBackLambertianNoAmbient.mat';     %autoExpTime =0.2361
% fullName = 'indObjFlashLessStrongBackLambertian.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertianCloser.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertianGT.mat'; 
% fullName = 'indObject2FlashDepth/backFlashImageLambertian0.mat'; 
fullName = 'indObject2FlashDepth/grayBackImage.mat'; 
load(fullName,'vci');
vciFlashBack = vci;
vcAddAndSelectObject('vcimage',vciFlashBack);
vcimageWindow;

flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
flashImageBack = imadjust(imageGet(vciFlashBack, 'results'), [], [], imageGet(vciFlash, 'gamma'));
figure; imshow(flashImage); figure; imshow(flashImageBack);

%ratio image
HSVFlashImage = rgb2hsv(flashImage);
HSVFlashImageBack = rgb2hsv(flashImageBack);

% linearIntensityFlash = imadjust(HSVFlashImage(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));
% linearIntensityFlashBack = imadjust(HSVFlashImageBack(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));

temp = imageGet(vciFlash, 'results');
linearIntensityFlash = sum(temp, 3);
temp = imageGet(vciFlashBack, 'results');
linearIntensityFlashBack = sum(temp, 3); 

figure; imshow(linearIntensityFlash)
figure; imshow(linearIntensityFlashBack)


ratioImage = linearIntensityFlash./linearIntensityFlashBack;
figure; imagesc(ratioImage);
fieldOfView = 39.60;



% new code for 2D

% alpha = asin(p(2)/d1); 
% phi   = asin(p(3)/(d1*cos(alpha)));

sensorWidth = 36;
sensorHeight  = 24;
sensorDistance = 11.8395;
% sensorDistance = 75;  %hmm this changes with changing focus...

xMatrix = linspace(-sensorWidth/2, sensorWidth/2, size(ratioImage,2));
xMatrix = repmat(xMatrix, [size(ratioImage,1), 1]);
yMatrix = linspace(-sensorHeight/2, sensorHeight/2, size(ratioImage,1))';
yMatrix = repmat(yMatrix, [1 size(ratioImage,2)]);
z = sensorDistance;

fakeD1 = sqrt(xMatrix.^2 + yMatrix.^2 + z.^2);
alpha = asin(xMatrix./fakeD1); 
phi   = asin(z./(fakeD1.*cos(alpha)));


f = 50; %5;
radical = abs(sqrt(4*cos(alpha).^2.*sin(phi).^2.*f.^2 - 4*f^2.*(1 - ratioImage)));
d1Test = (2.*f.^2)./(-2.*cos(alpha).*sin(phi).*f + radical);

figure; imagesc(d1Test);



%% estimate the surface normals using a depth map

%ignore boundary cases for now - add back in later if wanted

pixelUnit = sensorWidth / size(ratioImage,2);
% d1TestFiltered = imfilter(d1Test,fspecial( 'gaussian', 10, 10));
% d1TestFiltered = d1Test;
% d1TestFiltered = bilateralFilter(d1Test, 40, 20, 30);  %filtering operation for better depth map
d1TestMedFiltered = medianFilter(d1Test,15);
d1TestMedFiltered = medianFilter(d1TestMedFiltered',15)';

figure; imagesc(d1TestMedFiltered);
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
fakeDistance = sensorWidth/2 / atan(fieldOfView/2);

numWidth = size(d1Test,2);
numHeight = size(d1Test,1);



for i = 1:(size(d1TestFiltered, 2))
    for j = 1:(size(d1TestFiltered,1))
        tempVector = [pixelUnit * (i - numWidth/2)  pixelUnit * (j - numHeight/2) fakeDistance];
        tempVector = tempVector ./ norm(tempVector);
        rayVectors1(j, i, :) = reshape(tempVector, [1 1 3]);
        
        tempVector = [pixelUnit * (i - numWidth/2)  pixelUnit * (j - numHeight/2) fakeDistance + f];
        tempVector = tempVector ./ norm(tempVector);        
        rayVectors2(j, i, :) = reshape(tempVector, [1 1 3]);
        
        frontDot(j,i) = sum(rayVectors1(j,i,:) .* normalMap(j,i,:));
        backDot(j,i) = sum(rayVectors2(j,i,:) .* normalMap(j,i,:));
        
        %calculate correction factors for both images
        linearIntensityFlashCorrected(j,i, :) = linearIntensityFlashCorrected(j,i, :) ./ frontDot(j,i,: ) ; 
        linearIntensityFlashBCorrected(j,i, :) = linearIntensityFlashBCorrected(j,i, :) ./ backDot(j,i,:) ; 
    end
end


% debug print
% figure; imshow(rayVectors ./2 + .5);


ratioImage = linearIntensityFlashCorrected./linearIntensityFlashBCorrected;
figure; imagesc(ratioImage);


