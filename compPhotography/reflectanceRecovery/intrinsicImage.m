%% Intrinsic Image Estimation From Depth 
% Andy Lin
%
% This script uses ground truth depth map and normal map to create "intrinsic image"

%% Load image ( we will use the front flash image)

%Potential images to use.  Choose One!
% fullName = 'compPhotography/reflectanceRecovery/indObjOi.mat'; 
fullName = 'reflectanceRecovery/indObjSimpleRadianceFrontFlashImage.mat';

load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

%% Load ground truth depth map

import = load('2FlashDepth/indObject/depthMaps/groundTruthDepthMapDown.mat');
groundTruthDepthMap = import.depthMapProcessedMedian;
figure; imagesc(groundTruthDepthMap);
colorbar; title('Ground Truth Depth Map'); caxis([60 120]);

%% Divide image by the amount of dropoff due to depth

linearFlashImage = imageGet(vciFlash, 'results');
figure; imshow(linearFlashImage);
flashImageDepthCorrected = linearFlashImage.*(repmat(imresize(groundTruthDepthMap, [400 600]).^2, [1 1 3]));
figure; imshow(flashImageDepthCorrected./max(flashImageDepthCorrected(:)));

%% Estimate the surface normals using the ground truth depth map
% We perform surface normal estimation by taking the cross product of 2
% perpendicular vectors on the surface.  This vector is averaged amongst
% the 4 pairs of vectors around 1 point of the surface.

% calculate normal vectors using averaged cross products
d1TestFiltered = imresize(groundTruthDepthMap, [400 600]);
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

normalMapFiltered = normalMap;
for i = 1:3
    normalMapFiltered(:,:,i) = medianFilter(normalMap(:,:,i),5);
    normalMapFiltered(:,:,i) = medianFilter(normalMapFiltered(:,:,i)',5)';    
end

% normalMapFiltered = bilateralFilter(normalMapFiltered , 10, 4, 25);  

scaledNormalMap = normalMapFiltered./2 + .5;
% filteredNormals = imfilter(scaledNormalMap,fspecial( 'gaussian', 5, 5));

figure; imshow(scaledNormalMap);
title('Calculated Normal Map (from groundTruth)');



%% divide image by the amount of dropoff due to lambertian reflection attenuation (Correcting for Lambert's Law error)
% Depending on the relation of the surface normal, and the light direction,
% the intensity is attenuated according to the dot product (proportional to
% the cosine of the angle between the 2 vectors).  We will divide by the
% depth corrected image by this factor to obtain an intrinsic image.

flashImageDepthLambertianCorrected = flashImageDepthCorrected;
rayVectors1 = zeros([size(d1Test,1) size(d1Test,2) 3]);
frontDot =  zeros([size(d1Test,1) size(d1Test,2)]);

%assuming field of view given above
fakeDistance = sensorWidth/2 / tan(fieldOfView/2 * pi/180);
numWidth = size(flashImageDepthLambertianCorrected,2);
numHeight = size(flashImageDepthLambertianCorrected,1);

for i = 1:(size(flashImageDepthLambertianCorrected, 2))
    for j = 1:(size(flashImageDepthLambertianCorrected,1)) 
        % calculate normal vector from front light
        fakeX = pixelUnit * (i - numWidth/2);
        fakeY = pixelUnit * (j - numHeight/2);
        tempVector = [-fakeX fakeY fakeDistance];
        tempVector = tempVector ./ norm(tempVector);
        rayVectors1(j, i, :) = reshape(tempVector, [1 1 3]);
        
        frontDot(j,i) = sum(rayVectors1(j,i,:) .* normalMapFiltered(j,i,:));

        %calculate correction factors for both images
        flashImageDepthLambertianCorrected(j,i, :) = flashImageDepthLambertianCorrected(j,i, :) ./ frontDot(j,i) ; 
    end
end 
flashImageDepthLambertianCorrected(isnan(flashImageDepthLambertianCorrected)) = 0;
flashImageDepthLambertianCorrected(isinf(flashImageDepthLambertianCorrected)) = 0;
sortedIntensities = sort(flashImageDepthLambertianCorrected(:));
figure; imshow(flashImageDepthLambertianCorrected./sortedIntensities(round(length(sortedIntensities) * .9)));



%% load the ground truth intrinsic image by dividing by graycard image

%%% load the radiance
load('compPhotography/reflectanceRecovery/indObjSimpleRadianceFrontFlashOi.mat'); 
% load('compPhotography/reflectanceRecovery/indObjOi.mat'); 
irradianceOi = opticalimage;
irradianceOi = oiSet(irradianceOi,'name','Irradiance Image');
vcAddAndSelectObject(irradianceOi); oiWindow;

%%% load the "graycard" image
load('compPhotography/reflectanceRecovery/indObjSimpleGrayFrontOi.mat');  
illuminantOi = opticalimage;
illuminantOi = oiSet(illuminantOi, 'name', 'Graycard Image');
vcAddAndSelectObject(illuminantOi); oiWindow;

radianceValues = oiGet(irradianceOi , 'photons');
illuminantValues = oiGet(illuminantOi, 'photons');
overallIlluminantMean = mean(illuminantValues(:));

%%% calculate the ground truth reflectance by dividing the radiance by the graycard image 
 %assumptions: all lambertian surfaces, only 1 bounce allowed
reflectance = radianceValues./illuminantValues;
reflectance(isnan(reflectance)) = 0;
reflectance(isinf(reflectance)) = 0;

%show as an optical image of the reflectance
reflectanceOi = irradianceOi;
reflectanceOi = oiSet(reflectanceOi, 'cphotons', double(reflectance));
reflectanceOi = oiSet(reflectanceOi,'name','Reflectance');
vcAddAndSelectObject(reflectanceOi); oiWindow;

