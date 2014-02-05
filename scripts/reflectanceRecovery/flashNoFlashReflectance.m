%% Example of flash/no-flash surface reflectance estimation from PBRT rendered image
% Andy Lin
% 
% This script will go through the data flow of estimating surface
% reflectance, using flash/no-flash. Currently, it is only calculating the
% flash-only image from a flash/no-flash pair.

%set this flag to false if you want to use the pre-rendered scenes
reRenderScenes = false;
ambientFile = 'indestructibleObject/paintReflectance-downAmbient.pbrt';
flashFile = 'indestructibleObject/paintReflectance-downAmbientFlash.pbrt';

ambientFile = 'desk/paint.pbrt';
flashFile = 'desk/paintFlash.pbrt';


ambientOiFile = [s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjPaintAmbientOnly.mat']
flashOiFile = [s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjPaintAmbientAndFlash.mat']
ambientFlashExposureTime = 472.44*10^-3;
ambientExposureTime = 1.8;


% ambientOiFile = [s3dRootPath '/compPhotography/reflectanceRecovery/imageData/deskPaintAmbientOnly.mat']
% flashOiFile = [s3dRootPath '/compPhotography/reflectanceRecovery/imageData/deskPaintAmbientAndFlash.mat']
% ambientFlashExposureTime = 1.91;
% ambientExposureTime = 7.21;

%% Render/load scenes
if (reRenderScenes)
    %%% Scene rendering
    %%% Ambient only scene generation
    sceneName = 'Ambient Irradiance Image';
    oiAmbient = s3dRenderOI(ambientFile, .050);
    oiAmbient = oiSet(oiAmbient, 'name', sceneName);
    % strip the file name from the path and assign that as the name of the
    % object  ... vcSaveObject(oi,);
    vcAddAndSelectObject(oiAmbient);
    oiWindow;

    %%% Ambient and flash scene generation
    sceneName = 'Flash + Ambient Irradiance Image';
    oiFlash = s3dRenderOI(flashFile, .050);
    oiFlash = oiSet(oiFlash, 'name', sceneName);
    % strip the file name from the path and assign that as the name of the
    % object  ... vcSaveObject(oi,);
    vcAddAndSelectObject(oiFlash);
    oiWindow;
else
    %%% Load pre-rendered scenes
    %%% load the ambient only image
    load(ambientOiFile); 
    ambientOi = opticalimage;
    ambientOi = oiSet(ambientOi,'name','Ambient Irradiance Image');
    ambientOi = s3dFixOi(ambientOi, .050); %this may or may not be necessary, depending on the scene (if it was processed under old or new code)
    vcAddAndSelectObject(ambientOi); oiWindow;
    oiAmbient = ambientOi;

    %%% first load the flash + ambient image
    load(flashOiFile); 
    ambientFlashOi = opticalimage;
    ambientFlashOi = oiSet(ambientFlashOi,'name','Flash + Ambient Irradiance Image');
    vcAddAndSelectObject(ambientFlashOi); oiWindow;
    oiFlash = ambientFlashOi;
end

%% Process images

%%% process ambient image
% sensor processing
sensorAmbient = s3dProcessSensor(oiAmbient, 0, oiGet(oiAmbient, 'size'),  .0, 'analog');  %we might need to play with exposure
vcAddAndSelectObject('sensor',sensorAmbient); sensorImageWindow;
%image processing
[imageAmbient, transformMatrices] = s3dProcessImage(sensorAmbient, []);
vcAddAndSelectObject(imageAmbient); vcimageWindow;
ambientImage = imageAmbient;

%%% process ambient and flash image
% sensor processing
sensorFlash = s3dProcessSensor(oiFlash, 0, oiGet(oiFlash, 'size'),  .0, 'analog');  %we might need to play with exposure
vcAddAndSelectObject('sensor',sensorFlash); sensorImageWindow;
%image processing
[imageFlash, transformMatrices] = s3dProcessImage(sensorFlash, []);
vcAddAndSelectObject(imageFlash); vcimageWindow;
ambientFlashImage = imageFlash;


%%% process pre-computed images ( run this if you have already rendered the scenes and would like to load them from file

%  %% calculate the flash only image from the existing 2
% %note this is just a concept - we will do this again with the actual
% %processed images
% calculatedFlashOnlyOi = opticalimage;
% subtractionValue = double(oiGet(ambientFlashOi, 'photons') - oiGet(ambientOi, 'photons'));
% calculatedFlashOnlyOi = oiSet(calculatedFlashOnlyOi, 'photons', subtractionValue);
% vcAddAndSelectObject(calculatedFlashOnlyOi); oiWindow;

%  %% compare the computed flash only image with the ground truth
% load([s3dRootPath '/compPhotography/reflectanceRecovery/indObjSimpleRadiance2FlashOi.mat']); 
% gtFlashOnlyOi = opticalimage;
% gtFlashOnlyOi = oiSet(gtFlashOnlyOi,'name','Ground Truth Flash Irradiance Image');
% vcAddAndSelectObject(gtFlashOnlyOi); oiWindow;

%%% process ambient + flash image

% %fix Oi and form sensor image
% ambientFlashOi = s3dFixOi(ambientFlashOi, .050); %this may or may not be necessary, depending on the scene (if it was processed under old or new code)
% ambientFlashSensor = s3dProcessSensor(ambientFlashOi, 0, [], 0);   
% ambientFlashExposureTime = sensorGet(ambientFlashSensor, 'exposuretime');
% vcAddAndSelectObject(ambientFlashSensor); sensorImageWindow;
% %process image
% [ambientFlashImage, transformMatrices] = s3dProcessImage(ambientFlashSensor, []);
% vcAddAndSelectObject(ambientFlashImage); vcimageWindow;
% 
% %%% process ambient image
% 
% %fix Oi and form sensor image
% ambientSensor = s3dProcessSensor(ambientOi, 0, [], 0);   
% ambientExposureTime = sensorGet(ambientSensor, 'exposuretime');
% vcAddAndSelectObject(ambientSensor); sensorImageWindow;
% %process image
% [ambientImage, transformMatrices] = s3dProcessImage(ambientSensor, []);
% vcAddAndSelectObject(ambientImage); vcimageWindow;


%% Compute flash only image


equalizationRatio = ambientFlashExposureTime/ambientExposureTime; 
%this is the ratio of the ambient + flash exposure to the ambient exposure.
%To convert the ambient only image to the ambient + flash exposure, we must multiply
%the ambient image by this ratio to equalize it.  
calculatedFlashOnlyImage = ambientImage;
calculatedFlashOnlyResult = imageGet(ambientFlashImage, 'result') - imageGet(ambientImage, 'result') * equalizationRatio;
calculatedFlashOnlyImage = imageSet(calculatedFlashOnlyImage, 'result', calculatedFlashOnlyResult);
calculatedFlashOnlyImage = imageSet(calculatedFlashOnlyImage, 'name', 'Calculated Flash Only Image');
vcAddAndSelectObject(calculatedFlashOnlyImage); vcimageWindow;

%%% Compare calculated with actual flash only image

% %form sensor images
% %adjust for pbrt constant 
% gtFlashOnlyOi = oiAdjustIlluminance(gtFlashOnlyOi,.1);
% gtFlashOnlyOi = s3dFixOi(gtFlashOnlyOi, .050);
% gtFlashOnlySensor = s3dProcessSensor(gtFlashOnlyOi, 0, [], 0);   
% vcAddAndSelectObject(gtFlashOnlySensor); %sensorImageWindow;
% 
% %process sensor images
% [gtFlashOnlyImage, transformMatrices] = s3dProcessImage(gtFlashOnlySensor, []);
% vcAddAndSelectObject(gtFlashOnlyImage); vcimageWindow;
% gtFlashOnlyImageResult = imageGet(gtFlashOnlyImage, 'result') ;
% 
% calculatedFlashOnlyImageResult = imageGet(calculatedFlashOnlyImage, 'result');
% calculatedFlashOnlyImageResult = mean(gtFlashOnlyImageResult(:))/mean(calculatedFlashOnlyImageResult(:)) .* calculatedFlashOnlyImageResult;
% % calculatedFlashOnlyImageResult = 1./equalizationRatio .* calculatedFlashOnlyImageResult;
% 
% error = gtFlashOnlyImageResult - calculatedFlashOnlyImageResult;
% figure; imshow(error + .5); title ('Calculated Flash Image Error');
% squaredError = (error).^2;
% MSE = mean(squaredError(:))

%% Estimate the surface reflectance

% load basis functions
% imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/bases/simpleIndObj.mat']);
% imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/bases/paintIndObj.mat']);
% imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/bases/paintAllIndObj.mat']);
imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/bases/dupontGretagSG.mat']);
% imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/bases/dupontGretagSGNature.mat']);

bases = imported.bases;

numWave = 31; % import this somehow

% illuminant = .255 * ones(31, 1);
importedData = load([isetRootPath '/data/lights/' 'VivitarFlash.mat']);
desiredWave = 400:10:700;
SPD = ieReadSpectra('data/lights/VivitarFlash', desiredWave);
illuminant = SPD * 1.5;

pixel = sensorGet(sensorAmbient, 'pixel');
% %overall sensor response is spectral QE * fill factor 
% %.7 is a fudge factor for now since the unit conversion from PBRT to Iset
sensorResponse = sensorGet(sensorAmbient, 'spectral QE') * pixelGet(pixel, 'fill factor') * .7;  
imageData = imageGet(calculatedFlashOnlyImage, 'result');
surfaceReflectanceCalc = calculateReflectance(imageData, bases, illuminant, sensorResponse);

oi = oiCreate;
oi = initDefaultSpectrum(oi);
oi = oiSet(oi,'photons',surfaceReflectanceCalc);
oi = oiSet(oi, 'name', 'Calculated Reflectance');
vcAddAndSelectObject(oi); oiWindow;

% plot some example reflectance curves
imported = load([s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjReflectance.mat']);
reflectance = oiGet(imported.opticalimage, 'photons');

%plot red
plotPoints = surfaceReflectanceCalc(10:5:160,10:5:160,:);
figure; hold on;
pCalculated = plot(400:10:700, reshape(plotPoints(:), [size(plotPoints, 1) * size(plotPoints,2) 31]));
set(pCalculated,'Color','blue','LineWidth',1);

pGT = plot(400:10:700, reshape(reflectance(10,10,:), [1 31]));
set(pGT,'Color','green','LineWidth',3, 'LineStyle','--');
xlabel('Wavelength (nm)');
ylabel('Reflectance');
hold off;

%plot gray
plotPoints = surfaceReflectanceCalc(220:5:290, 10:5:220,:);
figure; hold on;
pCalculated = plot(400:10:700, reshape(plotPoints(:), [size(plotPoints, 1) * size(plotPoints,2) 31]));
set(pCalculated,'Color','blue','LineWidth',1);

pGT = plot(400:10:700, reshape(reflectance(240, 130,:), [1 31]));
set(pGT,'Color','green','LineWidth',3, 'LineStyle','--');
xlabel('Wavelength (nm)');
ylabel('Reflectance');
hold off;


%plot blue
plotPoints = surfaceReflectanceCalc(190:5:231,258:5:313,:);
figure; hold on;
pCalculated = plot(400:10:700, reshape(plotPoints(:), [size(plotPoints, 1) * size(plotPoints,2) 31]));
set(pCalculated,'Color','blue','LineWidth',1);

pGT = plot(400:10:700, reshape(reflectance(200,300,:), [1 31]));
set(pGT,'Color','green','LineWidth',3, 'LineStyle','--');
xlabel('Wavelength (nm)');
ylabel('Reflectance');
hold off;

%plot green
plotPoints = surfaceReflectanceCalc(45:8:165,360:8:432,:);
figure; hold on;
pCalculated = plot(400:10:700, reshape(plotPoints(:), [size(plotPoints, 1) * size(plotPoints,2) 31]));
set(pCalculated,'Color','blue','LineWidth',1);

pGT = plot(400:10:700, reshape(reflectance(120,400,:), [1 31]));
set(pGT,'Color','green','LineWidth',3, 'LineStyle','--');
xlabel('Wavelength (nm)');
ylabel('Reflectance');
hold off;

%% Estimate the ambient illuminant


numIlluminants = 3;
illuminantBank = cell(numIlluminants, 1);

illuminantBank{1} = illuminantCreate('tungsten', 400:10:700);
illuminantBank{2} = illuminantCreate('fluorescent', 400:10:700);
illuminantBank{3} = illuminantCreate('D65', 400:10:700);
estimatedIlluminantImage = zeros(size(imageData,1), size(imageData,2));
ambientImageRGB = imageGet(ambientImage, 'result');


for i = 1:1: size(surfaceReflectanceCalc,2)
    i
    for j = 1:1: size(surfaceReflectanceCalc,1)
        tempDistance = zeros(1, numIlluminants);
        for index = 1:numIlluminants
            tempRGB = sensorResponse' * diag(reshape(surfaceReflectanceCalc(j, i,: ), [numWave 1])) * ...
                illuminantGet(illuminantBank{index}, 'photons');
            illuminantBank{index}.rgbTraj(j,i,:) = reshape(tempRGB  , [1 1 3]) ;
            
            % using the distance from a point to line formula:  see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
            point0 = reshape(ambientImageRGB(j, i , :), [1 3]);
            point1 = [0 0 0];
            point2 = tempRGB' * 10^-10;
            
            vector1 = point0 - point1;
            vector2 = point0 - point2;
            vector3 = point2 - point1;
            
            tempDistance(index) =  norm(cross(vector1, vector2))/ norm(vector3);
            illuminantBank{index}.distance(j,i) = tempDistance(index); 
        end
        
        [minDistance minIndex] = min(tempDistance);
        estimatedIlluminantImage(j, i) = minIndex;
    end
end

% Pool the results for better results
pooledImage = ones(size(estimatedIlluminantImage));
blockSize = 8;
for i = 1:blockSize:size(estimatedIlluminantImage,2)-blockSize-1
    i
    for j = 1:blockSize:size(estimatedIlluminantImage,1)-blockSize-1
        tempBlock = estimatedIlluminantImage(j:j+blockSize-1, i:i+blockSize-1);
        pooledImage(j:j+blockSize-1,i:i+blockSize-1) = mode(tempBlock(:)) .* ones(blockSize,blockSize);
    end
end
figure; imagesc(pooledImage);
colorbar;

figure; imagesc(estimatedIlluminantImage);
colorbar;

% Calculate error
error = sum(not(estimatedIlluminantImage(:)==  1))

blockError = sum(not(pooledImage(:)==  1))

% figure; imagesc(illuminantBank{1}.distance);
% figure; imagesc(illuminantBank{2}.distance);
% figure; imagesc(illuminantBank{3}.distance);
% figure; imshow(double((illuminantBank{2}.rgbTraj))/max(illuminantBank{1}.rgbTraj(:)));
%% 
