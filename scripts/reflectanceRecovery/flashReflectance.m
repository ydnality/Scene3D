%% Example of flash surface reflectance estimation from PBRT rendered image
% Andy Lin
% 
% This script will go through the data flow of estimating surface
% reflectance, using flash only. The spectrum of the flash is known.  SVD
% is used to produced basis functions.  Currently, only 3 basis functions
% are used.  The basis functions are obtained using the ground truth
% reflectances of this scene.

%% Generate the ground truth reflectance
% This is performed by loading a "graycard" image, which has a reflectance
% of 1 for the entire spectrum, thus simulating a graycard for every single
% point of the scene.  The "graycard" image represents the illuminant 
% present at each point of the scene.  An irradiance image is then 
% produced, and divided by the "graycard" image to produce the ground truth
% reflectance.

%%% load the radiance
%load('compPhotography/reflectanceRecovery/indObjRadianceOi.mat'); 
% load('compPhotography/reflectanceRecovery/indObjSimpleRadiance2Oi.mat'); 
load([s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjPaintWhiteFlash.mat']); 
% load([s3dRootPath '/compPhotography/reflectanceRecovery/imageData/deskPaintWhiteFlash.mat']); 
irradianceOi = opticalimage;
irradianceOi = oiSet(irradianceOi,'name','Irradiance Image');
vcAddAndSelectObject(irradianceOi); oiWindow;

%%% load the "graycard" image
% load('compPhotography/reflectanceRecovery/indObjIlluminantOi.mat'); 
load([s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjGraycardFlash.mat']); 
% load([s3dRootPath '/compPhotography/reflectanceRecovery/imageData/deskGraycardFlash.mat']); 
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

%% Find the basis functions of this reflectance
% The SVD is used to compute basis functions.  In the future, we will use a
% database of reflectances, instead of computing basis functions from the
% same scene.

samples = reshape(reflectance, [size(reflectance,3) size(reflectance,1) * size(reflectance,2) ]);    %dimensions:  31 x numSamples
[bases,sValues,mn] = lmComputeBases(samples,0,3);
% bases'*bases %should be identity matrix
figure; plot(400:10:700, bases); grid on; title('Basis Functions'); xlabel('Wavelength(nm)'); ylabel('response');



%% Find the basis functions of a reflectance database
% The SVD is used to compute basis functions.  In the future, we will use a
% database of reflectances, instead of computing basis functions from the
% same scene.

% importedData = load([isetRootPath '/data/surfaces/reflectances/' 'DupontPaintChip_Vhrel.mat']);
% figure; plot(importedData.data); title('DupontPaintChip Vhrel.mat reflectances');

desiredWave = 400:10:700;
SPD = ieReadSpectra('/data/surfaces/reflectances/DupontPaintChip_Vhrel', desiredWave);

SPD2 = ieReadSpectra('/data/surfaces/reflectances/gretagDigitalColorSG', desiredWave);


SPD3 = ieReadSpectra('/data/surfaces/reflectances/Nature_Vhrel', desiredWave);
samples = cat(2, SPD, SPD2);

figure; plot(desiredWave, samples);

% samples = reshape(reflectance, [size(reflectance,3) size(reflectance,1) * size(reflectance,2) ]);    %dimensions:  31 x numSamples

[bases,sValues,mn] = lmComputeBases(samples,0,3);
% bases'*bases %should be identity matrix
figure; plot(400:10:700, bases); grid on; title('Basis Functions'); xlabel('wavelength'); ylabel('response');

%% Process the optical image
%%% make a sensor image from the irradiance image

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

%%% make an image from the sensor image (image processing)

[image, transformMatrices] = s3dProcessImage(sensor, []);
image = imageSet(image, 'name', sensorGet(sensor, 'name'));
vcAddAndSelectObject(image); vcimageWindow;
% imageData = imadjust(imageGet(image, 'result'), [],[],imageGet(image,
% 'gamma'));   %gamma is not needed in this case!

imageData = imageGet(image, 'result'); % / pixelGet(pixel, 'conversionGain');

%% calculate surface reflectance
illuminant = .255 * ones(31, 1);
pixel = sensorGet(sensor, 'pixel');
% %overall sensor response is spectral QE * fill factor 
% %.7 is a fudge factor for now since the unit conversion from PBRT to Iset
sensorResponse = sensorGet(sensor, 'spectral QE') * pixelGet(pixel, 'fill factor') * .7;  
surfaceReflectanceCalc = calculateReflectance(imageData, bases, illuminant, sensorResponse);
 
% 
% [imageData,r,c] = RGB2XWFormat(imageData);
% 
% % old code that RGB2XWFormat fixed
% % reshapedImageData = reshape(imageData, [size(imageData,1) * size(imageData,2) size(imageData, 3) ]);
% % reshapedImageData = reshapedImageData';
% 
% %% Reflectance calculation
% % Reflectance is calculated using the same method as the one used in
% % DiCarlo et al.  The sensor response function can be written in terms of
% % the surface reflectance, the ambient illuminant, and spectral response of
% % the sensor.  We then assume that the surface reflectance can be written
% % in terms of a weighted sum of basis functions.  The equation is inverted
% % to obtain the solution of the weights for the basis functions.  
% 
% % photodetectorResponse = pixelGet(pixel, 'pdspectralsr');
% % sensorResponse = sensorGet(sensor, 'colorfilters') .* repmat(photodetectorResponse, [1 3]);
% % colorTransform = imageGet(image, 'combinedTransform');
% sensorResponse = sensorGet(sensor, 'spectral QE') * pixelGet(pixel, 'fill factor') * .7;  
% %overall sensor response is spectral QE * fill factor 
% %.7 is a fudge factor for now since the unit conversion from PBRT to Iset
% %is a mystery.  This will change for different exposures... must adjust for
% %it
% 
% %relative spectral response of illuminant is difficult to determine
% % illuminant = .255 * ones(31, 1);
% illuminant = diag(ones(31,1)) * overallIlluminantMean;% * sensorGet(sensor, 'gain') * sensorGet(sensor, 'analog gain');
% surfaceReflectanceCalcXW = (bases * (((sensorResponse)' * illuminant * bases)^-1 * imageData'))';
% surfaceReflectanceCalc = XW2RGBFormat(surfaceReflectanceCalcXW,r,c);

% surfaceReflectanceCalc = reshape(surfaceReflectanceCalc', [size(imageData,1) size(imageData,2) size(surfaceReflectanceCalc, 1)]);
% surfaceReflectance(surfaceReflectance < 0) = 0;
% surfaceReflectance(surfaceReflectance > 1) = 1;
% surfaceReflectanceCalc = imresize(surfaceReflectanceCalc, [size(reflectance,1) size(reflectance,2)]);
%% Visualize Result

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',surfaceReflectanceCalc);
oi = oiSet(oi, 'name', 'Calculated Reflectance');
vcAddAndSelectObject(oi); oiWindow;

% plot some example reflectance curves
figure; plot(400:10:700, reshape(surfaceReflectanceCalc(10,10,:), [1 31]), 400:10:700, reshape(reflectance(10,10,:), [1 31]));

% reflectanceCalcOi = irradianceOi;
% reflectanceCalcOi = oiSet(reflectanceCalcOi, 'cphotons', double(surfaceReflectanceCalc));
% vcAddAndSelectObject(reflectanceCalcOi); oiWindow;