%% Example of surface reflectance estimation from PBRT rendered image
% Andy Lin


%from Joyce
% [csBasis,sValues] = lmComputeBases(RefMatrix,0,nBases);
% plot(csBasis); grid on;



%% first an extremely simple example
samples = [ 1 2 3 4; 2 3 4 5; 3 4 5 6]
samples = samples';
[bases,sValues,mn] = lmComputeBases(samples,0,1);
figure; plot(bases); grid on;

%there is only 1 distinct signal, so the basis should be a straight line 

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

%% image processing
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
