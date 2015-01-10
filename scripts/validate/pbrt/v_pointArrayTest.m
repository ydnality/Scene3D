%%  Creates a pointArray using point light sources
%
%  PBRT related
%
% The point light sources are created using tiny spherical area lights 
% diffraction is turned off in this example!
%
% AL VISTASOFT 2014

%initialization and directory organization
% filePath = [datapath '/validate/pbrtObject/']
% chdir(filePath);
% mkdir('batchPbrtFiles');
% unix('rm batchPbrtFiles/*');
% unix('cp * ./batchPbrtFiles/');
% chdir('batchPbrtFiles');

%make a new pbrt object
clear curPbrt;
curPbrt = pbrtObject();
curPbrt.camera.setLens(fullfile(dataPath, 'validate', 'pbrtObject', 'idealLensDiffraction-50mm.pbrt'));
curPbrt.removeLight();   %remove existing lights and geometry
curPbrt.removeGeometry();
curPbrt.removeMaterial();
curPbrt.sampler.removeProperty();
curPbrt.sampler.addProperty(pbrtPropertyObject('integer pixelsamples', 64));

% curPbrt.sampler.setPixelSamples(64);  %curPbrt.sampler.setPixelSamples(131070);
curPbrt.camera.setResolution(201, 201);

%add a grid of point light sources
xPos = -8:2:8;
yPos = -8:2:8;
for i = xPos
    for j = yPos
        clear pointSource;
        pointSource = pbrtLightAreaObject();
        pointSource.addTransform(pbrtTransformObject('Translate', [4.5 + i 0 7 + j]));
        pointSource.removeShape();
        pointSource.addShape(pbrtShapeObject('sphere', 'radius', .02));   %add a tiny spherical area light
        curPbrt.addLightSource(pointSource);
    end
end

%loop through different depths
for depth = 50:10:130

    curPbrt.camera.setPosition([4.5 -depth 7; % Starting up
                    4.5 -depth + 1 7; % Ending up
                    0 0 1]);

    %run pbrt
    tmpFileName = fullfile(dataPath, 'generatedPbrtFiles', ['deleteMe' int2str(depth) '.pbrt']);
    curPbrt.writeFile(tmpFileName);
%     oi = s3dRenderOI(tmpFileName, 50, [filePath '/batchPbrtFiles/']);
    oi = s3dRenderOI(tmpFileName, .050, tmpFileName);
end


% chdir('..');
%the result should be a grid of points



