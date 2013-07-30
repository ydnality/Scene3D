%% Example of flash/no-flash surface reflectance estimation from PBRT rendered image
% Andy Lin

%% generate the ground truth reflectance
%load the radiance
%load('compPhotography/reflectanceRecovery/indObjRadianceOi.mat'); 
load('compPhotography/reflectanceRecovery/indObjSimpleRadiance2Oi.mat'); 

irradianceOi = opticalimage;
irradianceOi = oiSet(irradianceOi,'name','Irradiance Image');
vcAddAndSelectObject(irradianceOi); oiWindow;

 %load the "graycard" image
load('compPhotography/reflectanceRecovery/indObjIlluminantOi.mat'); 
illuminantOi = opticalimage;
illuminantOi = oiSet(illuminantOi, 'name', 'Graycard Image');
vcAddAndSelectObject(illuminantOi); oiWindow;

radianceValues = oiGet(irradianceOi , 'photons');
illuminantValues = oiGet(illuminantOi, 'photons');
overallIlluminantMean = mean(illuminantValues(:));

 %calculate the ground truth reflectance by dividing the radiance by the graycard image 
 %assumptions: all lambertian surfaces, only 1 bounce allowed
reflectance = radianceValues./illuminantValues;
reflectance(isnan(reflectance)) = 0;
reflectance(isinf(reflectance)) = 0;

 %show as an optical image of the reflectance
reflectanceOi = irradianceOi;
reflectanceOi = oiSet(reflectanceOi, 'cphotons', double(reflectance));
reflectanceOi = oiSet(reflectanceOi,'name','Reflectance');
vcAddAndSelectObject(reflectanceOi); oiWindow;

%% find the basis functions of this reflectance
samples = reshape(reflectance, [size(reflectance,3) size(reflectance,1) * size(reflectance,2) ]);    %dimensions:  31 x numSamples
[bases,sValues,mn] = lmComputeBases(samples,0,3);
% bases'*bases %should be identity matrix
figure; plot(bases); grid on;


%% Make a sensor image from the irradiance image
%process image 
oi = oiAdjustIlluminance(irradianceOi,.1);
% oi = oiSet(oi, 'photons', oiGet(irradianceOi,'photons') * 10^13);  %some normalization issues
optics = oiGet(oi, 'optics');  %to create proper crop at sensor
optics = opticsSet(optics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', optics);
oi = oiSet(oi,'fov', 39.60); %fov of the scene

%sensor processing
sensor = s3dProcessSensor(oi, 0, [], 0);   
sensor = sensorSet(sensor, 'name', oiGet(oi, 'name'));
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

%% make an image from the sensor image (image processing)
[image, transformMatrices] = s3dProcessImage(sensor, []);
image = imageSet(image, 'name', sensorGet(sensor, 'name'));
vcAddAndSelectObject(image); vcimageWindow;
% imageData = imadjust(imageGet(image, 'result'), [],[],imageGet(image,
% 'gamma'));   %gamma is not needed in this case!
pixel = sensorGet(sensor, 'pixel');
imageData = imageGet(image, 'result'); % / pixelGet(pixel, 'conversionGain');
[imageData,r,c] = RGB2XWFormat(imageData);

% old code that RGB2XWFormat fixed
% reshapedImageData = reshape(imageData, [size(imageData,1) * size(imageData,2) size(imageData, 3) ]);
% reshapedImageData = reshapedImageData';

%% reflectance calculation
% photodetectorResponse = pixelGet(pixel, 'pdspectralsr');
% sensorResponse = sensorGet(sensor, 'colorfilters') .* repmat(photodetectorResponse, [1 3]);
% colorTransform = imageGet(image, 'combinedTransform');
sensorResponse = sensorGet(sensor, 'spectral QE') * pixelGet(pixel, 'fill factor') * .7;  
%overall sensor response is spectral QE * fill factor 
%.7 is a fudge factor for now since the unit conversion from PBRT to Iset
%is a mystery.  This will change for different exposures... must adjust for
%it

%relative spectral response of illuminant is difficult to determine
% illuminant = .255 * ones(31, 1);
illuminant = diag(ones(31,1)) * overallIlluminantMean;% * sensorGet(sensor, 'gain') * sensorGet(sensor, 'analog gain');
surfaceReflectanceCalc = bases * (((sensorResponse)' * illuminant * bases)^-1 * imageData');
surfaceReflectanceCalc = XW2RGBFormat(surfaceReflectanceCalc',r,c);

% surfaceReflectanceCalc = reshape(surfaceReflectanceCalc', [size(imageData,1) size(imageData,2) size(surfaceReflectanceCalc, 1)]);
% surfaceReflectance(surfaceReflectance < 0) = 0;
% surfaceReflectance(surfaceReflectance > 1) = 1;
% surfaceReflectanceCalc = imresize(surfaceReflectanceCalc, [size(reflectance,1) size(reflectance,2)]);
oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',surfaceReflectanceCalc);
oi = oiSet(oi, 'name', 'Calculated Reflectance');
vcAddAndSelectObject(oi); oiWindow;

% reflectanceCalcOi = irradianceOi;
% reflectanceCalcOi = oiSet(reflectanceCalcOi, 'cphotons', double(surfaceReflectanceCalc));
% vcAddAndSelectObject(reflectanceCalcOi); oiWindow;

%% flash/no-flash recovery of reflectance 

%first load the flash + ambient image
load([s3dRootPath '/compPhotography/reflectanceRecovery/indObjSimpleRadiance2FlashAmbientOi.mat']); 
ambientFlashOi = opticalimage;
ambientFlashOi = oiSet(ambientFlashOi,'name','Flash + Ambient Irradiance Image');
vcAddAndSelectObject(ambientFlashOi); oiWindow;

%load the ambient only image
load([s3dRootPath '/compPhotography/reflectanceRecovery/indObjSimpleRadiance2AmbientOi.mat']); 
ambientOi = opticalimage;
ambientOi = oiSet(ambientOi,'name','Ambient Irradiance Image');
vcAddAndSelectObject(ambientOi); oiWindow;

%calculate the flash only image from the existing 2
%note this is just a concept - we will do this again with the actual
%processed images
calculatedFlashOnlyOi = opticalimage;
subtractionValue = double(oiGet(ambientFlashOi, 'photons') - oiGet(ambientOi, 'photons'));
calculatedFlashOnlyOi = oiSet(calculatedFlashOnlyOi, 'photons', subtractionValue);
vcAddAndSelectObject(calculatedFlashOnlyOi); oiWindow;

%compare the computed flash only image with the ground truth
load([s3dRootPath '/compPhotography/reflectanceRecovery/indObjSimpleRadiance2FlashOi.mat']); 
gtFlashOnlyOi = opticalimage;
gtFlashOnlyOi = oiSet(gtFlashOnlyOi,'name','Ground Truth Flash Irradiance Image');
vcAddAndSelectObject(gtFlashOnlyOi); oiWindow;

%% process ambient + flash image
%form sensor images
%adjust for pbrt constant 
% originalMeanIlluminance  = oiGet(ambientFlashOi, 'meanIlluminance');
% ambientFlashOi = oiAdjustIlluminance(ambientFlashOi,.1);
% conversionRatio = .1/originalMeanIlluminance;  %this gives a conversion from PBRT to ISET illuminance

conversionRatio = 10^13;
ambientFlashOi = oiSet(ambientFlashOi, 'photons', oiGet(ambientFlashOi,'photons') * 10^13);  %some normalization issues

optics = oiGet(ambientFlashOi, 'optics');  %to create proper crop at sensor
optics = opticsSet(optics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
ambientFlashOi = oiSet(ambientFlashOi, 'optics', optics);
ambientFlashOi = oiSet(ambientFlashOi,'fov', 39.60); %fov of the scene
ambientFlashSensor = s3dProcessSensor(ambientFlashOi, 0, [], 0);   
ambientFlashSensor = sensorSet(ambientFlashSensor, 'name', oiGet(ambientFlashOi, 'name'));
ambientFlashExposureTime = sensorGet(ambientFlashSensor, 'exposuretime');
vcAddAndSelectObject(ambientFlashSensor); sensorImageWindow;

%process sensor images
[ambientFlashImage, transformMatrices] = s3dProcessImage(ambientFlashSensor, []);
ambientFlashImage = imageSet(ambientFlashImage, 'name', sensorGet(ambientFlashSensor, 'name'));
vcAddAndSelectObject(ambientFlashImage); vcimageWindow;

%% process ambient image
%form sensor images
%adjust for pbrt constant 
% ambientOi = oiAdjustIlluminance(ambientOi,.1);
ambientOi = oiSet(ambientOi, 'photons', oiGet(ambientOi, 'photons') * conversionRatio);
% oi = oiSet(oi, 'photons', oiGet(irradianceOi,'photons') * 10^13);  %some normalization issues
optics = oiGet(ambientOi, 'optics');  %to create proper crop at sensor
optics = opticsSet(optics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
ambientOi = oiSet(ambientOi, 'optics', optics);
ambientOi = oiSet(ambientOi,'fov', 39.60); %fov of the scene
ambientSensor = s3dProcessSensor(ambientOi, 0, [], 0);   
ambientSensor = sensorSet(ambientSensor, 'name', oiGet(ambientOi, 'name'));
ambientExposureTime = sensorGet(ambientSensor, 'exposuretime');
vcAddAndSelectObject(ambientSensor); sensorImageWindow;

%process sensor images
[ambientImage, transformMatrices] = s3dProcessImage(ambientSensor, []);
ambientImage = imageSet(ambientImage, 'name', sensorGet(ambientSensor, 'name'));
vcAddAndSelectObject(ambientImage); vcimageWindow;


%% compute flash only image
equalizationRatio = ambientFlashExposureTime/ambientExposureTime; 
%this is the ratio of the ambient + flash exposure to the ambient exposure.
%To convert the ambient only image to the ambient + flash exposure, we must multiply
%the ambient image by this ratio to equalize it.  
calculatedFlashOnlyImage = ambientImage;
calculatedFlashOnlyResult = imageGet(ambientFlashImage, 'result') - imageGet(ambientImage, 'result') * equalizationRatio;
calculatedFlashOnlyImage = imageSet(calculatedFlashOnlyImage, 'result', calculatedFlashOnlyResult);
calculatedFlashOnlyImage = imageSet(calculatedFlashOnlyImage, 'name', 'Calculated Flash Only Image');
vcAddAndSelectObject(calculatedFlashOnlyImage); vcimageWindow;

%compare with actual flash only image

%form sensor images
%adjust for pbrt constant 
gtFlashOnlyOi = oiAdjustIlluminance(gtFlashOnlyOi,.1);
% oi = oiSet(oi, 'photons', oiGet(irradianceOi,'photons') * 10^13);  %some normalization issues
optics = oiGet(gtFlashOnlyOi, 'optics');  %to create proper crop at sensor
optics = opticsSet(optics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
gtFlashOnlyOi = oiSet(gtFlashOnlyOi, 'optics', optics);
gtFlashOnlyOi = oiSet(gtFlashOnlyOi,'fov', 39.60); %fov of the scene
gtFlashOnlySensor = s3dProcessSensor(gtFlashOnlyOi, 0, [], 0);   
gtFlashOnlySensor = sensorSet(gtFlashOnlySensor, 'name', oiGet(gtFlashOnlyOi, 'name'));
% ambientExposureTime = sensorGet(gtFlashOnlySensor, 'exposuretime');
vcAddAndSelectObject(gtFlashOnlySensor); sensorImageWindow;

%process sensor images
[gtFlashOnlyImage, transformMatrices] = s3dProcessImage(gtFlashOnlySensor, []);
gtFlashOnlyImage = imageSet(gtFlashOnlyImage, 'name', sensorGet(gtFlashOnlySensor, 'name'));
vcAddAndSelectObject(gtFlashOnlyImage); vcimageWindow;
gtFlashOnlyImageResult = imageGet(gtFlashOnlyImage, 'result') ;

calculatedFlashOnlyImageResult = imageGet(calculatedFlashOnlyImage, 'result');
calculatedFlashOnlyImageResult = mean(gtFlashOnlyImageResult(:))/mean(calculatedFlashOnlyImageResult(:)) .* calculatedFlashOnlyImageResult;
% calculatedFlashOnlyImageResult = 1./equalizationRatio .* calculatedFlashOnlyImageResult;

error = gtFlashOnlyImageResult - calculatedFlashOnlyImageResult;
figure; imshow(error + .5);
squaredError = (error).^2;
MSE = mean(squaredError(:))