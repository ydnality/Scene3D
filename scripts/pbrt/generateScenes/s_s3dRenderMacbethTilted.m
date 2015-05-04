%% Render the macbeth color chiecker as a 3D object, tilted
%
% AL Vistasoft, 2015

%%
ieInit

%% Rendering parameters
addLight     = true;
% theta = 30 * pi/180;  % Angle of tilt
theta1 = 1;
theta2 = 0.8;
theta3 = 1.2;

%% initialization
clear curPbrt;
curPbrt = pbrtObject();

% Candidate for a function
curPbrt.removeMaterial();
curPbrt.removeGeometry();
curPbrt.removeLight();

%% Add a light source

% light properties
wave = 400:10:700;
light = ieReadSpectra('D65',wave, []);
spectrum = Energy2Quanta(wave, light); %convert to photons
spectrum = reshape([wave(:) spectrum(:)]',1,2*length(wave));
spectrumObject = pbrtSpectrumObject('spectrum L', spectrum(:));
lightFront = pbrtLightDistantObject('light',spectrumObject, [0 0 80], [0 0 79]);
% curPbrt.addLightSource(lightFront);

% Alternative light
lightFrom = [ 0 0 0];    % Position of source
lightTo =   [0 0 -1 ];   % Direction of principal ray
coneAngle      = 180;    % Angle of rays from light source
coneDeltaAngle = 180;    % Drop off of illumination???

lightSource = pbrtLightSpotObject('light', spectrumObject, coneAngle, coneDeltaAngle, lightFrom, lightTo);
curPbrt.addLightSource(lightSource);

%% Add spot light and infinite light sources

if addLight
    % Make the spot light
    disp('Making spot light')
    spotLight = pbrtLightSpotObject();
    spotLight.setName('spot');
    spotLight.setSpectrum(pbrtSpectrumObject('rgb I', [1000 1000 1000]));
    spotLight.setAngle(90);      % 180
    spotLight.setDeltaAngle(90); %180
    spotLight.setFrom([-142.3855 -286.2024  13.0082]);
    spotLight.setTo([ -141.9582 -285.2984   13.0082]);
    curPbrt.addLightSource(spotLight);
    
    %infinite light (for diffuse lighting)
    infiniteLight = pbrtLightInfiniteObject();
    curPbrt.addLightSource(infiniteLight);
    
end

%% Standard camera properties

%uncomment to use a 2ElLens instead
from = [ 0 0 0];
to =   [0 0 -1];
position = [from; to; 0 1 0 ];
%  Need to test the rank or something of the positon matrix
% rank(position) cannot be 1
% the third row determines the up vector for the camera

curPbrt.camera.setPosition(position);

% Two element lens example
% curPbrt.camera.setLens(pbrtLensRealisticObject());
% curPbrt.camera.lens.filmDistance = 133.33; %90; % 70;  % 133.33;
% curPbrt.camera.lens.filmDiag = 70;
% curPbrt.camera.lens.specFile = '2ElLens.dat';
% curPbrt.camera.lens.apertureDiameter = 16; % in mm
% curPbrt.camera.lens.curveRadius = 0;       % Experimental
% curPbrt.camera.setResolution(450, 300);

% This is a pinhole camera (the default)
% Not assigning a specFile makes it a pinhole camera.
newCamPos =    [0  0 0;
    0   0 -1;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);
curPbrt.camera.lens.filmDistance = 133.33;
curPbrt.camera.lens.filmDiag = 70;

%% Backdrop
backDropDepth = -160;
matRGB = [400 1 500 1 600 1 700 1 ]';
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('spectrum Kd', matRGB));
curPbrt.addMaterial(newMaterial);

%add a backdrop
backDropTransform = ...
    [50 0 0 0;
    0 50 0 0 ;
    0 0 1 0;
    0 0 backDropDepth  1];
backDrop = pbrtGeometryObject('backdrop', 'grayMat', [], [], backDropTransform);
curPbrt.addGeometry(backDrop);

%% MCC

macbethSpectrumFile = fullfile(isetRootPath, 'data', 'surfaces', 'macbethChart.mat');
[reflectances,wave,comment,fName]  = ieReadSpectra(macbethSpectrumFile, 400:10:700, []);

scaleFactor = 4;

foregroundDepth = -80;
xValues = linspace(-2.5*scaleFactor, 2.5*scaleFactor, 6);
yValues = linspace(-1.5*scaleFactor, 1.5*scaleFactor, 4);
[xOffset, yOffset] = meshgrid(xValues, yValues);

for index = 1:24
    spectrum= [wave;
        reflectances(:, index)'];
    spectrumObject = pbrtPropertyObject('spectrum Kd', spectrum(:));
    newMaterial = pbrtMaterialObject(['macbeth' int2str(index)], 'matte', spectrumObject);
    curPbrt.addMaterial(newMaterial);
end

for ii = 1:6
    for jj = 1:4
        %add a foreground target
        
        foregroundTransform = ...
            [2 0 0 0;
            0 2 0 0 ;
            0 0 1 0;
            xOffset(jj,ii) yOffset(jj,ii) foregroundDepth  1];
        
        % Think about clarying this, like XY, ZY, YZ rotations and so
        % forth. Maybe a function?
        rotationTransform1 = ...
            [1 0 0 0;
            0 cos(theta1) -sin(theta1) 0;
            0 sin(theta1) cos(theta1) 0;
            0 0 0 1]';
        
        rotationTransform2 = ...
            [cos(theta2) -sin(theta2) 0 0;
             sin(theta2) cos(theta2)  0 0; 
               0            0       1 0;
            0 0 0 1]';
        
        rotationTransform3 = ...
            [cos(theta3)  0 -sin(theta3)  0; 
            0      1      0       0;
            sin(theta3) 0 cos(theta3)   0;
            0 0 0 1]';
        
        totalTransform = rotationTransform1*rotationTransform2 * rotationTransform3 * foregroundTransform;
        
        frontSquare = pbrtGeometryObject(['checker' int2str(ii) int2str(jj)],...
            ['macbeth' int2str(3 -(jj -1) + (ii-1) * 4 + 1)], [], [], ...
            totalTransform);
        curPbrt.addGeometry(frontSquare);
        
        %frontSquare = pbrtGeometryObject(['checker' int2str(ii) int2str(jj)], ['macbeth'], [], [], totalTransform);
        
        
    end
end

%% set sampler
curPbrt.sampler.removeProperty();
nSamples = 256;                      % 512;
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', nSamples));

%% Render the oi

dockerFlag = true;
oi = s3dRenderOIAndDepthMap(curPbrt, 'Macbeth', dockerFlag);
vcAddAndSelectObject(oi); oiWindow;

%% OLD CODE

%% Render Scene Radiance Using pbrtObjects (make all checkers completely white)
% for i = 1:14
%     
%     clear curPbrt;
%     curPbrt = pbrtObject();
%     
%     %camera position
%     newCamPos =    [0  0 0;
%         0   0 -1;
%         0 1.00000 0];
%     curPbrt.camera.setPosition(newCamPos);
%     curPbrt.camera.lens.filmDistance = 133.33;
%     curPbrt.camera.lens.filmDiag = 70;
%     
%     scaleFactor = 4;
%     % curPbrt.camera.setResolution(100, 100);    %LQ mode
%     
%     %uncomment to use a 2 element lens instead of a pinhole
%     % curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));
%     
%     %sampler
%     sampler = curPbrt.sampler.removeProperty();
%     sampler.value = 256;
%     curPbrt.sampler.addProperty(sampler);
%     
%     %backdrop Depth
%     % backDropDepth = -100 * scaleFactor;  %backdrop distance increases with depth of spheres
%     backDropDepth = -160;
%     foregroundDepth = -80;
%     foregroundDepth2 = -70;
%     foregroundDepth3= -90;
%     
%     %calculate sphere offsets
%     xValues = linspace(-2.5*scaleFactor, 2.5*scaleFactor, 6);
%     yValues = linspace(-1.5*scaleFactor, 1.5*scaleFactor, 4);
%     [xOffset yOffset] = meshgrid(xValues, yValues);
%     
%     % lightSpectrumFile = fullfile(s3dRootPath, 'papers', 'ReflectanceAndDepth', 'Illuminant.mat');
%     [lights,wave,comment,fName]  = ieReadSpectra('D65', 400:10:700, []);
%     lights = Energy2Quanta(wave, lights); %convert to photons
%     %res is returned as a 31 x 14 matrix, where it's rows: wavelength and cols:
%     
%     %light source
%     % lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
%     curPbrt.removeLight();
%     spectrum = lights(:, i)';
%     %spectrum = ones(size(spectrum));  %temp debug
%     tempMatrix = [400:10:700; spectrum];  %helps put data in 400 1 500 .5 600 .5 700 1 format
%     spectrumObject = pbrtSpectrumObject('spectrum I', tempMatrix(:));
%     %%for infinite light source
%     %spectrumObject = pbrtSpectrumObject('spectrum I', tempMatrix(:));     %for finite light source
%     lightFront = pbrtLightSpotObject(['light' int2str(i)], spectrumObject, [], [], [0 0 0], [0 0 -1]);
%     %lightFront = pbrtLightDistantObject(['light' int2str(i)],spectrumObject, [0 0 80], [0 0 79]);
%     curPbrt.addLightSource(lightFront);
%     
%     %add a new material
%     matRGB= [400 1 500 1 600 1 700 1 ]';
%     newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('spectrum Kd', matRGB));
%     curPbrt.addMaterial(newMaterial);
%     
%     %add material file
%     curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'materials', 'simpleTarget-mat.pbrt'));
%     
%     % remove default geometry
%     curPbrt.removeGeometry();
%     %add a backdrop
%     backDropTransform = ...
%         [50 0 0 0;
%         0 50 0 0 ;
%         0 0 1 0;
%         0 0 backDropDepth  1];
%     
%     backDrop = pbrtGeometryObject('backdrop', 'Material', [], [], backDropTransform);
%     curPbrt.addGeometry(backDrop);
%     
%     %read macbeth color checker refletance values
%     macbethSpectrumFile = fullfile(isetRootPath, 'data', 'surfaces', 'macbethChart.mat');
%     [reflectances,wave,comment,fName]  = ieReadSpectra(macbethSpectrumFile, 400:10:700, []);
%     
%     %add new material for macbeth color checker reflectances
%     %     for index = 1:24
%     %         spectrum= [wave;
%     %                    reflectances(:, index)'];
%     %         spectrumObject = pbrtPropertyObject('spectrum Kd', spectrum(:));
%     %         newMaterial = pbrtMaterialObject(['macbeth' int2str(index)], 'matte', spectrumObject);
%     %         curPbrt.addMaterial(newMaterial);
%     %     end
%     
%     % for ii = 1:6
%     %     for jj = 1:4
%     %add a foreground target
%     foregroundTransform = ...
%         [10 0 0 0;
%         0 10 0 0 ;
%         0 0 1 0;
%         0 0 foregroundDepth  1];
%     
%     rotationTransform = ...
%         [1 0 0 0;
%         0 cos(theta) -sin(theta) 0;
%         0 sin(theta) cos(theta) 0;
%         0 0 0 1]';
%     
%     totalTransform = rotationTransform * foregroundTransform;
%     %totalTransform = rotationTransform;
%     %frontSquare = pbrtGeometryObject(['checker' int2str(ii) int2str(jj)], ['macbeth'], [], [], totalTransform);
%     
%     %uncomment to make all targets white
%     frontSquare = pbrtGeometryObject(['checker'], ['grayMat'], [], [], totalTransform);
%     
%     curPbrt.addGeometry(frontSquare);
%     
%     % noScale = true;
%     dockerFlag = false;
%     focalLength = 0.050;
%     oi = s3dRenderOIAndDepthMap( curPbrt, focalLength, 'simpleScene', dockerFlag);
%     
%     vcAddObject(oi); oiWindow;
%     %fullName = vcSaveObject(scene, fullfile(s3dRootPath, 'papers', 'ReflectanceAndDepth', 'Data', '03052015_finiteLight', ['whitelight' int2str(i) '.mat']));
%     % fullName = vcSaveObject(scene, fullfile(s3dRootPath, 'papers', 'ReflectanceAndDepth', 'Data', '02062015_scenes', ['whitelight' int2str(i) '.mat']));
% end