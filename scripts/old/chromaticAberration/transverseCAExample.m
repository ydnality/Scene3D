%% Example of flash/no-flash surface reflectance estimation from PBRT rendered image
% Andy Lin
% 
% This script will go through the data flow of estimating surface
% reflectance, using flash/no-flash. Currently, it is only calculating the
% flash-only image from a flash/no-flash pair.

%set this flag to false if you want to use the pre-rendered scenes
reRenderScenes = true;
inputFile = 'testTargetCA/realistic2PinkTest.pbrt';

preloadedInputFile = '';

% ambientOiFile = [s3dRootPath '/compPhotography/reflectanceRecovery/imageData/indObjPaintAmbientOnly.mat']


%% Render/load scenes
if (reRenderScenes)
    %%% Scene rendering
    %%% Ambient only scene generation
    sceneName = 'Chromatic Aberration Image';
    oi = s3dRenderScene(inputFile, .050);
    oi = oiSet(oi, 'name', sceneName);
    % strip the file name from the path and assign that as the name of the
    % object  ... vcSaveObject(oi,);
    vcAddAndSelectObject(oi);
    oiWindow;

else
    %%% Load pre-rendered scenes
    %%% load the ambient only image
    load(ambientOiFile); 
    ambientOi = opticalimage;
    ambientOi = oiSet(ambientOi,'name','Ambient Irradiance Image');
    ambientOi = s3dFixOi(ambientOi, .050); %this may or may not be necessary, depending on the scene (if it was processed under old or new code)
    vcAddAndSelectObject(ambientOi); oiWindow;

end


%oiAdjustIlluminance can be used to scale the irradiance values
% oi = oiSet(oi, 'photons', oiGet(oi, 'photons') * 10^-12);
% vcAddAndSelectObject(oi); oiWindow;

% shorterDistance = oi;
data = load('noApertureOffset091713.mat');
shorterDistance = data.opticalimage;
irradianceShort = oiGet(shorterDistance, 'photons');
spreadShort = mean(irradianceShort(:, 160:190, 1),2)/mean(irradianceShort(:));

% longerDistance = oi;
data = load('apertureOffset1p5Above091713.mat');
longerDistance = data.opticalimage;
irradianceLong = oiGet(longerDistance, 'photons');
spreadLong = mean(irradianceLong(:, 160:190, 1), 2)/mean(irradianceLong(:));
% figure; plot(1:400,spreadShort);
% figure; plot(1:400,spreadLong);
figure; plot(1:400,spreadShort,'g', 1:400,spreadLong , '--m');





% 
% irradianceShort = oiGet(shorterDistance, 'photons');
% spreadShort = mean(irradiance(:, 160:190, 1),2);
% 
% % longerDistance = oi;
% irradianceLong = oiGet(longerDistance, 'photons');
% spreadLong = mean(irradiance(:, 160:190, 1), 31);
% % figure; plot(1:400,spreadShort);
% % figure; plot(1:400,spreadLong);
% figure; plot(1:400,spreadShort,'g', 1:400,spreadLong , '--m');