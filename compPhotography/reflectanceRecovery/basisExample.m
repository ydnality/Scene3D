
%from Joyce
% [csBasis,sValues] = lmComputeBases(RefMatrix,0,nBases);
% plot(csBasis); grid on;



%% first an extremely simple example

samples = [ 1 2 3 4; 2 3 4 5; 3 4 5 6]
samples = samples';
[bases,sValues,mn] = lmComputeBases(samples,0,1);
figure; plot(bases); grid on;

%there is only 1 distinct signal, so the basis should be a straight line 

%% a more complex example

%%generate the ground truth reflectance
 %load the radiance
% load('compPhotography/reflectanceRecovery/indObjRadianceOi.mat'); 
load('compPhotography/reflectanceRecovery/indObjSimpleRadianceOi.mat'); 

radianceOi = opticalimage;
vcAddAndSelectObject(radianceOi); oiWindow;

 %load the "graycard" image
load('compPhotography/reflectanceRecovery/indObjIlluminantOi.mat'); 
illuminantOi = opticalimage;
vcAddAndSelectObject(illuminantOi); oiWindow;

radianceValues = oiGet(radianceOi , 'photons');
illuminantValues = oiGet(illuminantOi, 'photons');

 %calculate the ground truth reflectance by dividing the radiance by the graycard image 
 %assumptions: all lambertian surfaces, only 1 bounce allowed
reflectance = radianceValues./illuminantValues;
reflectance(isnan(reflectance)) = 0;
reflectance(isinf(reflectance)) = 0;

 %show as an optical image of the reflectance
reflectanceOi = radianceOi;
reflectanceOi = oiSet(reflectanceOi, 'cphotons', double(reflectance * 10^20));
vcAddAndSelectObject(reflectanceOi); oiWindow;

%%find the basis functions of this reflectance
samples = reshape(reflectance, [size(reflectance,3) size(reflectance,1) * size(reflectance,2) ]);    %dimensions:  31 x numSamples
[bases,sValues,mn] = lmComputeBases(samples,0,3);
figure; plot(bases); grid on;


%%process image 
oi = oiSet(oi, 'photons', oiGet(radianceOi,'photons') * 10^13);  %some normalization issues
myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
myOptics = opticsSet(myOptics,'focalLength', .050);  %distance from lens to sensor - sensor will inherit this value!
oi = oiSet(oi, 'optics', myOptics);
oi = oiSet(oi,'fov', 39.60); %fov of the scene

 %sensor processing
sensor = s3dProcessSensor(oi, 0, [], [], .32);   %test scene - back 100m flash - multiplication factor - 8
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

 %image processing
[image, transformMatrices] = s3dProcessImage(sensor, []);
vcAddAndSelectObject(image); vcimageWindow;
imageData = imageGet(image, 'result');
reshapedImageData = reshape(imageData, [size(imageData,1) * size(imageData,2) size(imageData, 3) ]);

reshapedImageData = reshapedImageData';

%%reflectance calculation
sensorResponse = sensorGet(sensor, 'colorfilters');
colorTransform = imageGet(image, 'combinedTransform');
surfaceReflectanceCalc = bases * ((sensorResponse)' * bases)^-1 * colorTransform^-1 * reshapedImageData;

surfaceReflectanceCalc = reshape(surfaceReflectanceCalc', [size(imageData,1) size(imageData,2) size(surfaceReflectanceCalc, 1)]);
% surfaceReflectance(surfaceReflectance < 0) = 0;
% surfaceReflectance(surfaceReflectance > 1) = 1;
surfaceReflectanceCalc = imresize(surfaceReflectanceCalc, [size(reflectance,1) size(reflectance,2)]);

reflectanceCalcOi = radianceOi;
reflectanceCalcOi = oiSet(reflectanceCalcOi, 'cphotons', double(surfaceReflectanceCalc * 10^20));
vcAddAndSelectObject(reflectanceCalcOi); oiWindow;
