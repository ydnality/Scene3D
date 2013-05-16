%load flash image
% fullName = 'benchFlash.mat';
fullName = 'indObjFlash.mat';
% fullName = 'indObjFlashLessStrong.mat';
% fullName = 'deskFlash.mat';
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

%load no-flash image
% fullName = 'benchNoFlash.mat';
% fullName = 'indObjNoFlash.mat';
fullName = 'indObjNoFlashLessNoise.mat';
% fullName = 'indObjFlashLessStrong.mat';
% fullName = 'deskNoFlash.mat';
%fullName = 'deskTest.mat';   %testing exposure things
load(fullName,'vci');
vciNoFlash = vci;
vcAddAndSelectObject('vcimage',vciNoFlash);
vcimageWindow;


%% obtain the ratio image

flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
HSVFlashImage = rgb2hsv(flashImage);
intensityImage = HSVFlashImage(:,:,3);

HSVNoFlashImage = rgb2hsv(noFlashImage);
intensityImageNF = HSVNoFlashImage(:,:,3);

ratioImage = intensityImage./intensityImageNF;
figure; imagesc(ratioImage);

%% obtain flash only image - use this to help with finding depth

flashExposure = 1.5020e-35;
noFlashExposure = 0.0138;

%flash only image can help build reflectance image - problem is that depth
%should be used because flash intensity diminishes with distance - so it is
%better to know depth in order to better estimate the reflectance.



%% use 2 flash distances to estimate distance
% fullName = 'indObjFlashLessStrong.mat';       %autoExpTime =  0.2230
fullName = 'indObjFlashLambertianNoAmbient.mat'; 
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

fullName = 'indObjFlashBackLambertianNoAmbient.mat';     %autoExpTime =0.2361
% fullName = 'indObjFlashLessStrongBackLambertian.mat'; 
load(fullName,'vci');
vciFlashBack = vci;
vcAddAndSelectObject('vcimage',vciFlashBack);
vcimageWindow;

flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
flashImageBack = imadjust(imageGet(vciFlashBack, 'results'), [], [], imageGet(vciFlash, 'gamma'));

%ratio image
HSVFlashImage = rgb2hsv(flashImage);
HSVFlashImageBack = rgb2hsv(flashImageBack);

linearIntensityFlash = imadjust(HSVFlashImage(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));
linearIntensityFlashBack = imadjust(HSVFlashImageBack(:,:,3), [], [], 1/imageGet(vciFlash, 'gamma'));

ratioImage = linearIntensityFlash./linearIntensityFlashBack;
figure; imagesc(ratioImage);

%normalizeRatioImage
% ratioImage = ratioImage./.15;
% figure; imagesc(ratioImage);

%experiment with converting to depth
% sqrtRatio = sqrt(ratioImage);
% depth = 1 ./ (sqrtRatio - 1.01);
% figure; imagesc(depth); 
% figure; mesh(depth);

% found solution for depth, given phi (angle), d1/d2 = 1/sqrt(ratio), and (1/ratio),

%first calculate phi at each part of the image.  assume a horizontal field
%of view of 39.60 as used in the rendering

fieldOfView = 39.60;
theta = linspace(-fieldOfView/2, fieldOfView/2, size(ratioImage, 2));
phi = 90 - abs(theta);
phi = repmat(phi, [size(ratioImage,1) 1]);
py = 200;
phiRadians = phi * pi/180;
d1overd2 = 1./sqrt(ratioImage);
d1overd22 = 1./ratioImage;
overd2solution1 = (2 .* sin(phiRadians) + sqrt(4 .* d1overd22 .* sin(phiRadians) .* sin(phiRadians) - 4 .* py^2 * (d1overd22 - 1)))./(2 .* py^2);
overd2solution2 = (2 .* sin(phiRadians) - sqrt(4 .* d1overd22 .* sin(phiRadians) .* sin(phiRadians) - 4 .* py^2 * (d1overd22 - 1)))./(2 .* py^2);

figure; imagesc(abs(1./overd2solution1));
figure; imagesc(abs(1./overd2solution2));
