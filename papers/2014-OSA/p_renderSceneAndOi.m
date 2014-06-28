%% 2014 OSA Conference 
%
% This script renders several versions of the indestructible Object scene
% using the reverse calculation.
%
% The first half renders these scenes using our modified version of PBRT
% (required for it to run properly).
%
% The second half renders the processed images (through the sensor and image
% processor), given the optical image .mat file and displays the results.
%
% AL Vistalab 2014

%% TOFIX more comments

%% Render optical images of various versions of indObj using PBRT 
    %% Render the defocused oi
    sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainDefocused.pbrt');
    oi = s3dRenderOi(sceneName, .050, 'indObj');

    %% Render the pinhole scene

    %the pinhole scene can be treated as the scene radiance
    sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinhole.pbrt');
    oi = s3dRenderOi(sceneName, .050, 'indObj');

    %% Render the depth map

    %the pinhole scene can be treated as the scene radiance
    sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainDepthMap.pbrt');
    % oi = s3dRenderOi(sceneName, .050, 'indObj');

    depthMap = s3dRenderDepthMap(sceneName, 1);
    figure; imagesc(depthMap);
    title('Depth-map For Indestructible Object (units mm)');
    colorbar;
    colormap('gray');

%% Create images for figures in paper with image slices of the previously rendered optical images
    %% create figure images for defocused oi

    % load pre-rendered oi's
    blurredOiFName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indObj50mm10mmAp.mat');
    load(blurredOiFName);
    blurredOi = opticalimage;
    vcAddObject(blurredOi);
    oiWindow;

    %get the photons and set some parameters
    blurredPhotons = oiGet(blurredOi, 'photons');
    sliceChannels = round(linspace(1, 31, 4));  %channels we wish to use
    wave = oiGet(blurredOi, 'wave');
    colorCode = [1 0 1;   %this determines the color of the slice images
                 0 0 1;
                 0 1 0
                 1 0 0];
    %get the 400nm, 500n, 600nm, and 700nm bands of the oi
    for i = 1:4
        imageSlice = blurredPhotons(:,:, sliceChannels(i));
        coloredSlice = imageSlice/max(imageSlice(:) * .65);
        coloredSlice = repmat(coloredSlice, [1 1 3]) .* repmat(reshape(colorCode(i,:), [1 1 3]), size(imageSlice));
        figure; imshow(coloredSlice);
        imwrite(coloredSlice, [blurredOiFName num2str(wave(sliceChannels(i))) 'nm.png']);
    end

    %% create figure images for pinhole image
    % load pre-rendered oi's
    pinholeOiFName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indObj50mmPinhole.mat');
    load(pinholeOiFName);
    pinholeOi = opticalimage;
    vcAddObject(pinholeOi);
    oiWindow;

    %get the photons and set some parameters
    pinholePhotons = oiGet(pinholeOi, 'photons');
    sliceChannels = round(linspace(1, 31, 4));   %channels we wish to use
    wave = oiGet(pinholeOi, 'wave');
    colorCode = [1 0 1;   %this determines the color of the slice images
                 0 0 1;
                 0 1 0
                 1 0 0];
    %get the 400nm, 500n, 600nm, and 700nm bands of the oi
    for i = 1:4
        imageSlice = pinholePhotons(:,:, sliceChannels(i));
        coloredSlice = imageSlice/max(imageSlice(:) * .65);
        coloredSlice = repmat(coloredSlice, [1 1 3]) .* repmat(reshape(colorCode(i,:), [1 1 3]), size(imageSlice));
        figure; imshow(coloredSlice);
        imwrite(coloredSlice, [blurredOiFName num2str(wave(sliceChannels(i))) 'nm.png']);
    end

%% Render the processed images using the ISET pipeline, given a saved optical image .mat file
    %% create sensor images
    readNoise = 0;
    oi = pinholeOi;  
    % oi = blurredOi;
    sensor = s3dProcessSensor(oi, readNoise, [500 500],[], 'analog');    %low noise, auto exposure
    vcAddObject(sensor); sensorWindow;

    %% create processed image
    image = s3dProcessImage(sensor);
    image = imageSet(image, 'gamma', .6);
    image = imageSet(image, 'internalcolorspace', 'XYZ');
    image = vcimageCompute(image,sensor);

    vcAddAndSelectObject(image); vcimageWindow;

    %save images
    finalImage = imadjust(imageGet(image, 'results'), [], [], imageGet(image, 'gamma'));
    figure; imshow(finalImage);
    imwrite(finalImage, 'processedOutput.png');

