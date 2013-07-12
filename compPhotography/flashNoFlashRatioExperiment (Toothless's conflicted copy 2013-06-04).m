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
fullName = 'indObjNoFlash.mat';
% fullName = 'deskNoFlash.mat';
%fullName = 'deskTest.mat';   %testing exposure things
load(fullName,'vci');
vciNoFlash = vci;
vcAddAndSelectObject('vcimage',vciNoFlash);
vcimageWindow;


%% 

flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));
noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
HSVFlashImage = rgb2hsv(flashImage);
intensityImage = HSVFlashImage(:,:,3);

HSVNoFlashImage = rgb2hsv(noFlashImage);
intensityImageNF = HSVNoFlashImage(:,:,3);

ratioImage = intensityImage./intensityImageNF;
figure; imagesc(ratioImage);