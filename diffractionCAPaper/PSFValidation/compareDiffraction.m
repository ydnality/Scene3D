% Use the point test file, blur it in 2 different ways.  The first way is
% to use the new ray tracing method which uses Heisenburg Uncertainty Ray
% Bending (HURB).  The second way is the classical way, using theoretical PSF's.  

%% HURB ray tracing results
%for this example only

chdir(PSFValidationPath);

sampleArray = cell(1, 1);

sampleArray{2}.rayTraceFile = 'rayTraceScratch.mat' 
sampleArray{2}.focalLength = 50
sampleArray{2}.filmDistance = 51.281394156842644
sampleArray{2}.apertureDiameter = 2.2727
sampleArray{2}.targetDistance = 2


% sampleArray{1}.rayTraceFile = 'rayTrace25mm32res.mat' 
% sampleArray{1}.focalLength = 25
% sampleArray{1}.apertureDiameter = 1.1364
% sampleArray{1}.filmDistance = 25.3163
% sampleArray{1}.targetDistance = 2
% 
% sampleArray{2}.rayTraceFile = 'rayTrace50mm32res.mat' 
% sampleArray{2}.focalLength = 50
% sampleArray{2}.filmDistance = 51.281394156842644
% sampleArray{2}.apertureDiameter = 2.2727
% sampleArray{2}.targetDistance = 2
% 
% sampleArray{3}.rayTraceFile = 'rayTrace100mm32res.mat' 
% sampleArray{3}.focalLength = 100
% sampleArray{3}.apertureDiameter = 4.5455
% sampleArray{3}.filmDistance = 105.2604
% sampleArray{3}.targetDistance = 2
% 
% sampleArray{4}.rayTraceFile = 'rayTrace25mm4m32res.mat' 
% sampleArray{4}.focalLength = 25
% sampleArray{4}.apertureDiameter = 1.1364
% sampleArray{4}.filmDistance = 25.1572
% sampleArray{4}.targetDistance = 4
% 
% sampleArray{5}.rayTraceFile = 'rayTrace50mm4m32res.mat' 
% sampleArray{5}.focalLength = 50
% sampleArray{5}.filmDistance = 50.6329
% sampleArray{5}.apertureDiameter = 2.2727
% sampleArray{5}.targetDistance = 4
% 
% sampleArray{6}.rayTraceFile = 'rayTrace100mm4m32res.mat' 
% sampleArray{6}.focalLength = 100
% sampleArray{6}.apertureDiameter = 4.5455
% sampleArray{6}.filmDistance = 102.5641
% sampleArray{6}.targetDistance = 4
% 
% sampleArray{7}.rayTraceFile = 'rayTrace25mm1m32res.mat' 
% sampleArray{7}.focalLength = 25
% sampleArray{7}.apertureDiameter = 1.1364
% sampleArray{7}.filmDistance = 25.6410
% sampleArray{7}.targetDistance = 1
% 
% sampleArray{8}.rayTraceFile = 'rayTrace50mm1m32res.mat' 
% sampleArray{8}.focalLength = 50
% sampleArray{8}.filmDistance = 52.6316
% sampleArray{8}.apertureDiameter = 2.2727
% sampleArray{8}.targetDistance = 1
% 
% sampleArray{9}.rayTraceFile = 'rayTrace100mm1m32res.mat' 
% sampleArray{9}.focalLength = 100
% sampleArray{9}.apertureDiameter = 4.5455
% sampleArray{9}.filmDistance = 111.1111
% sampleArray{9}.targetDistance = 1

for index = 2:2 %length(sampleArray)
    
    %initialize variables
    load(sampleArray{index}.rayTraceFile);
    filmDistance = sampleArray{index}.filmDistance;
    sensorWidth = .2/sqrt(2)
    horFieldofView = 2 * atan(sensorWidth/(2 * filmDistance)) * 180/pi
    focalLength = sampleArray{index}.focalLength;

    %plot raytrace PSF
    oi = opticalimage;
    oiIlluminance = oiGet(oi, 'illuminance');
    PSFLine = oiIlluminance(size(oiIlluminance,1)/2, :);
    oi = oiSet (oi, 'horizontalfieldofview', horFieldofView);
    vcAddAndSelectObject(oi);
    position = linspace(-sensorWidth/2 *1000 , sensorWidth/2 *1000 , length(PSFLine));  % the range is different from theoretical because the theoretical result is a square image
%     figure;
%     plot(position, PSFLine);
%     title('Raytracing PSF');
%     xlabel('um');
%     ylabel('Illuminance');


    %% Theoretical results
    %scene = sceneCreate('point array',256,128);
    %scene = sceneSet(scene,'fov',8);
    %vcAddAndSelectObject(scene); sceneWindow;

    %load scene file
    scene = sceneFromFile('pointTest.png', 'rgb');
    %scene = sceneFromFile('usairforce.png', 'rgb');

    scene = sceneSet(scene,'fov',horFieldofView);
    scene = sceneSet(scene, 'distance', 2001);
    vcAddAndSelectObject(scene); sceneWindow;

    %create optical image
    oiT = oiCreate;
    optics = oiGet(oiT,'optics');           
    optics = opticsSet(optics,'fnumber',22);
    % In this example we set the properties of the optics to include cos4th
    % falloff for the off axis vignetting of the imaging lens
    optics = opticsSet(optics,'offaxis','cos4th');
    optics = opticsSet(optics,'focallength',focalLength * 10^-3);    
    oiT = oiSet(oiT,'optics',optics);
    oiT = oiCompute(scene,oiT);
    vcAddAndSelectObject(oiT); oiWindow;

    %plot both PSFs on 1 figure
    oiIlluminanceT = oiGet(oiT, 'illuminance');
    PSFLineT = oiIlluminanceT(size(oiIlluminanceT,1)/2, :);
%     PSFLineTS = PSFLineT * max(PSFLine(:))/max(PSFLineT(:));
    PSFLineTS = PSFLineT /max(PSFLineT(:));
    PSFLineS = PSFLine / max(PSFLine(:));
    positionT = linspace(-sensorWidth/2 *1000, sensorWidth/2 *1000, length(PSFLineT));
    figure;
    plot(position, PSFLineS, positionT, PSFLineTS);
    title(['PSF Comparison;' num2str(sampleArray{index}.focalLength) 'mm;f/' ...
        num2str(sampleArray{index}.focalLength/(sampleArray{index}.apertureDiameter), 2) ';' ...
        num2str(sampleArray{index}.targetDistance) 'm Target Distance' ]);
    xlabel('um');
    ylabel('Relative illuminance');
    legend Ray-tracing Theoretical
    
    fileName = ['PSFC_' num2str(sampleArray{index}.focalLength) 'mm_f' ...
        num2str(sampleArray{index}.focalLength/(sampleArray{index}.apertureDiameter), 2) '_' ...
        num2str(sampleArray{index}.targetDistance) 'mTargDis' ];
    hgexport(gcf, [fileName '.tif'], hgexport('factorystyle'), 'Format', 'tiff');
end