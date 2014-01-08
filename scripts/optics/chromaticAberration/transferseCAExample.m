%% Example of flash/no-flash surface reflectance estimation from PBRT rendered image
% Andy Lin
% 
% This script will go through the data flow of estimating surface
% reflectance, using flash/no-flash. Currently, it is only calculating the
% flash-only image from a flash/no-flash pair.

%set this flag to false if you want to use the pre-rendered scenes
reRenderScenes = true;
inputFile = 'testTargetCA/realisticPointTest.pbrt';

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
