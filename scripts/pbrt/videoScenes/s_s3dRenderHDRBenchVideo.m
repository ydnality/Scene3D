% Set camera position
goodCamPos = [179.701691 -170.845230 36.446182 
       178.975494 -170.161819 36.371479 
       -0.046965 0.059093 0.997147];
   
v1 = [179.701691 -170.845230 36.446182 ];
v2 = [ 178.975494 -170.161819 36.371479];

direction = v2 - v1;  
direction = direction./norm(direction, 2); %normalize this

orthoDirection = [-direction(2) direction(1) 0];


directionMat = [direction;
                direction;
                0 0 0];

orthoDirectionMat = [orthoDirection;
                     orthoDirection;
                     0 0 0];
 
%these positions are for reference only
initialCamPos = goodCamPos - directionMat * 200;
intermediateCamPos = goodCamPos + directionMat * 100;
finalCamPos = goodCamPos + directionMat * 100 + orthoDirectionMat * 150;

numStepsZoom = 150;
camPosZoom = zeros(3,3, numStepsZoom);
%linspace over the whole matrix to create intermediate matrices for zoom
%steps
for i = 1:3
    for j = 1:3
        camPosZoom(i, j, :) = linspace(initialCamPos(i,j), intermediateCamPos(i,j), numStepsZoom);
    end
end

numStepsPan = 150;
camPosPan = zeros(3,3, numStepsPan);
%linspace over whole matrix to create intermediate matrices for pan steps
for i = 1:3
    for j = 1:3
        camPosPan(i, j, :) = linspace(intermediateCamPos(i,j), finalCamPos(i,j), numStepsPan);
    end
end

%concatenate camPos matrix
camPos = cat(3, camPosZoom, camPosPan);

for i = 1:size(camPos, 3)
    %% Render Scene Radiance Using pbrtObjects
    clear curPbrt;
    curPbrt = pbrtObject();
    
    curPbrt.camera.setPosition(camPos(:,:,i));
    curPbrt.camera.lens.filmDistance = 70; %133.33;
    curPbrt.camera.lens.filmDiag = 70;
    curPbrt.camera.setResolution(300, 300);    %LQ mode
    
    %uncomment to use a 2 element lens instead of a pinhole
    % curPbrt.camera.setLens(fullfile(s3dRootPath, 'data', 'lens', '2ElLens50mm.pbrt'));
    
    % Sampler
    sampler = curPbrt.sampler.removeProperty();
    sampler.value = 1024;
    curPbrt.sampler.addProperty(sampler);
    
    % Backdrop Depth
    backDropDepth = -100;
    foregroundDepth = -65;
    
    % Light source
    curPbrt.addLightSource(fullfile(s3dRootPath,'data','pbrtScenes','benchScene', 'sunsetEnvironmentLight.pbrt'));
    
    % Add material file
    curPbrt.addMaterial(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'benchScene', 'default-mat.pbrt'));
    
    % Remove default geometry
    curPbrt.removeGeometry();
    
    % Add geometry
    curPbrt.addGeometry(fullfile(s3dRootPath, 'data', 'pbrtScenes', 'benchScene','default-geom-big-bigfloor.pbrt'));
    
    % Render scene and depth map
    scene = s3dRenderSceneAndDepthMap(curPbrt, 'simpleScene', true);
    
    %vcAddObject(scene); sceneWindow;
    %instead of visualizing scenes, we will save it
    fullName = vcSaveObject(scene, fullfile(dataPath, 'pbrtScenes', 'benchScene', 'HDRVideo', ['frame' int2str(i) '.mat']));
end