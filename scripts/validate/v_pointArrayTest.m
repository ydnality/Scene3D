%%  Creates a pointArray using point light sources.  
% The point light sources are created using tiny spherical area lights 
% diffraction is turned off in this example!

%initialization and directory organization
filePath = [datapath '/validate/pbrtObject/']
chdir(filePath);
mkdir('batchPbrtFiles');
unix('rm batchPbrtFiles/*');
unix('cp * ./batchPbrtFiles/');
chdir('batchPbrtFiles');

%make a new pbrt object
clear curPbrt;
curPbrt = pbrtObject();
curPbrt.camera.setLens('idealLensDiffraction-50mm.pbrt');
curPbrt.removeLight();   %remove existing lights and geometry
curPbrt.removeGeometry();
curPbrt.removeMaterial();
curPbrt.sampler.setPixelSamples(131070);
curPbrt.camera.setResolution(201, 201);

%add a grid of point light sources
xPos = -8:2:8;
yPos = -8:2:8;
for i = xPos
    for j = yPos
        clear pointSource;
        pointSource = lightAreaObject();
        pointSource.addTransform(transformObject('Translate', [4.5 + i 0 7 + j]));
        pointSource.removeShape();
        pointSource.addShape(shapeObject('sphere', 'radius', .02));   %add a tiny spherical area light
        curPbrt.addLightSource(pointSource);
    end
end

%run pbrt
tmpFileName = ['deleteMePointArray' '.pbrt'];
curPbrt.writeFile(tmpFileName);
oi = s3dRenderScene(tmpFileName, 50, [filePath '/batchPbrtFiles/'], tmpFileName);

chdir('..');

%the result should be a grid of points



