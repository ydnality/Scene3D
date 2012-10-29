% Script requirements: latest version of ISET.  Please contact imageEVAL to
% obtain this version of ISET.
% This script demonstrates non-uniform illumination functionality.  First,
% it loads a PBRT scene.  The PBRT illuminant image is then loaded in as a
% point-by-point varying illuminant.  We then plot the reflectance of the
% 550 nm band.  Compared to the uniform illuminant, the nonuniform
% illuminant gives superior estimation of reflectances, which are invariant
% to lighting changes, something which was not present previously.


%% Input Scene Illumination From PBRT
fName = 'desert.dat';
fNameIlluminant = 'desert-white.dat';

%% load the scene

fID = fopen(fName,'r','l');
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 Inf]);
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
photons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

%% Make a default scene structure.
fileWave = 400:10:700;   %hardcoded for now - will change this later
scene = sceneCreate;
scene = sceneSet(scene,'name',fName);
scene = sceneSet(scene, 'wave', fileWave);

%imgphotons = Energy2Quanta(filewave,img(:, :, :));
scene = sceneSet(scene,'photons',photons);

%% load the space-varying illuminant
fID = fopen(fNameIlluminant,'r','l');
[A, cnt] = fscanf(fID,'%d %d %d\n',[3 Inf]);
photons = fread(fID,prod(A),'double');

%A(2) = rows, A(1) = columns
illuminantPhotons = reshape(photons,A(2),A(1),A(3));
fclose(fID);

scene = sceneSet(scene,'illuminant photons',illuminantPhotons);
scene = sceneIlluminantScale(scene);  

vcAddAndSelectObject(scene); sceneWindow;

%retrieve the reflectance
reflectance = sceneGet(scene, 'reflectance');
%plot the 500 channel of the reflectance
figure; imagesc(reflectance(:,:,16));
title('Reflectance Calculated Using Nonuniform Illuminant');

%% demonstrate what the reflectance would have looked like with a uniform illuminant
d65 = vcReadSpectra('D65',fileWave);  % D65 in units of energy
scene = sceneSet(scene, 'illuminantEnergy', d65);
scene = sceneIlluminantScale(scene);  

vcAddAndSelectObject(scene); sceneWindow;

%retrieve the reflectance
reflectance = sceneGet(scene, 'reflectance');
%plot the 500 channel of the reflectance
figure; imagesc(reflectance(:,:,16));
title('Reflectance Calculated Using Uniform Illuminant');
