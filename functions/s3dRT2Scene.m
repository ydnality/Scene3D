function scene = s3dRT2Scene(fName,dMapFile,wave,dsample)
% Converts a renderToolbox radiance image into an ISET scene
% 
%  scene = s3dRT2Scene(fName,dMapFile,[wave=400:10:700],[dsample=1])
% 
% The RTB scene is calculated by???
% 
% Inputs:
%  fName:  Full path file name of the rendered radiance
%  dMapFile:  Corresponding depth map file - there may be problems if depth map dimensions are not square.
%  wave:   Sample wavelengths in the radiance file
%  dsample: Spatially downsample amount for scene.  Can be useful for speed/memory.
% 
% Return
%  scene: Returned ISET scene.  
% 
% Example:
%  fName = fullfile(s3dRootPath,'picMat_orange_rad.mat');
%  scene = s3dRT2Scene(fName,[],[],4);
%  vcAddAndSelectObject(scene); sceneWindow
% 
%  scene = s3dRT2Scene(fName,[],[],2);
%  vcAddAndSelectObject(scene); sceneWindow
%
% Copyright, Stanford, 2011

%% Check arguments
if ieNotDefined('fName'), fName = vcSelectDataFile; end
if ieNotDefined('wave'), wave = 400:10:700; end
if ieNotDefined('dsample'), dsample = 1;   end % Don't down sample

% Use the default.  But this only works with the picMat_orange_rad file.
if ieNotDefined('dMapFile')
    if strcmp(fName,fullfile(s3dRootPath,'picMat_orange_rad.mat'))
        dMapFile = fullfile(s3dRootPath,'depthMap.zbf'); 
    else
        error('No dMapFile specified');
    end
end


%% Load the .mat file produced by the RTB
tmp   = load(fName);
[r,c] = size(tmp.picMat{1});
w = length(tmp.picMat);

% Spatially downsample if requested.  This is needed for speed and memory
% in some cases.
if dsample == 1
    img = zeros(r,c,w);
    for ii = 1:w, img(:,:,ii) = tmp.picMat{ii}; end
else
    test = tmp.picMat{1}(1:dsample:r,1:dsample:c);
    img = zeros(size(test,1),size(test,2),w);
    for ii = 1:w, 
        img(:,:,ii) = tmp.picMat{ii}(1:dsample:r,1:dsample:c); 
    end
end


%% Make a default scene structure.
scene = sceneCreate;
scene = sceneSet(scene,'name',fName);
scene = sceneSet(scene,'wave',wave);

% We must check whether the data are in photons or energy
% scene = sceneSet(scene,'photons',img);
% Andy:  the units are in energy, determined by trial and error
imgphotons = Energy2Quanta(wave,img);
scene = sceneSet(scene,'photons',imgphotons);

% scene = sceneSet(scene,'knownReflectance',[img(10,10,10),10,10,10]);
% scene = sceneSet(scene,'illuminantEnergy',ones(nWave,1));
d65 = vcReadSpectra('D65',wave);  % D65 in units of energy  
% plot(wave,d65)
scene = sceneSet(scene, 'illuminantEnergy', d65);

% I don't understand this code.  So I took it out.  Probably we need to do
% something to make sure the reflectance and illuminant values scale
% properly.  We don't, so heavens knows what is happening now.
%
%  load('knownReflectance.mat');
%  scene = sceneSet(scene, 'knownReflectance', knownReflectance);
%  scene.data.knownReflectance = [];  
% Andy:for some reason, the pre-made scenes have a known reflectance set
%already, which are incorrect.  This was put into place so that
%sceneIlluminantScale would work properly and scale the data correctly.  In
%the future, to be more precise, we can go through each scene, and specify
%the known reflectance.  One more way to be more precise with illumination
%is to specify the illuminant on a pixel-by-pixel basis.
scene = sceneIlluminantScale(scene);  % Make reflectances reasonable.
scene = sceneAdjustLuminance(scene,100);

%% The depth map file needs to be created, specified, and then read.
% This only works for the one special case in the scene3D repository.

%tmp = load(dMapFile);
% What are the units on this depth?
% scene = sceneSet(scene,'depth map',tmp.test);
scene.depthMap = s3dReadDepthMapFile(dMapFile);
scene.depthMap = imresize(scene.depthMap, [sceneGet(scene, 'rows') sceneGet(scene, 'cols')]);

return;

        