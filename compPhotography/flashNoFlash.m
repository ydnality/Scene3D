% Flash/noflash experiment - basedly loosely off of "Flash Photography
% Enhancement via Instrinsic Relighting" by Eisenmann, and Durand.

%load flash image
%fullName = 'benchFlash.mat';
%fullName = 'indObjFlash.mat';
fullName = 'deskFlash.mat';
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

%load no-flash image
%fullName = 'benchNoFlash.mat';
% fullName = 'indObjNoFlash.mat';
fullName = 'deskNoFlashLessNoise.mat';
load(fullName,'vci');
vciNoFlash = vci;
vcAddAndSelectObject('vcimage',vciNoFlash);
vcimageWindow;


%% flash image processing

%separate into color and intensity

% this commented out section is implementing a modification mentioned in
% the appendix.  However, they don't specify how to get the color data, so
% we won't use this for now until we figure this out.

% tempSum = redChannel + greenChannel + blueChannel;
% intensityImage = redChannel./tempSum .* redChannel + ...
%                  greenChannel./tempSum .* greenChannel + ...
%                  blueChannel./tempSum .* blueChannel;

% use YCbCr colorspace to decouple color and intensity
flashImage = imadjust(imageGet(vciFlash, 'results'), [], [], imageGet(vciFlash, 'gamma'));

YCbCrFlashImage = rgb2ycbcr(flashImage);
intensityImage = YCbCrFlashImage(:,:,1);
figure; imshow(intensityImage); title('Flash Intensity Image');

%separate into details and large scale

%large scale
%smallImage = imresize(intensityImage, .5);  %put this line back in if you
%want to debug on small images
largeScale = bilateralFilter(intensityImage, 5, .025, 5);  %play with parameters if needed
figure; imshow(largeScale); title('Large Scale Image');

%details
details = intensityImage./largeScale;  %intensityImage = largeScale * details
figure; imshow(details - .5); title('Details');

%verification - commented out for now
% intensityImageRecovery = largeScale .* details;
% figure; imshow(intensityImageRecovery); title('Intensity Recovery');
% sum(sum(intensityImage - intensityImageRecovery))   %checks to see if near 0
%% no flash image processing

noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
YCbCrNoFlashImage = rgb2ycbcr(noFlashImage);
intensityImageNF = YCbCrNoFlashImage(:,:,1);
figure; imshow(intensityImageNF); title('No Flash Intensity Image');

%separate into details and large scale

%large scale
%smallImage = imresize(intensityImage, .5);  %put this line back in if you
%want to debug on small images
largeScaleNF = bilateralFilter(intensityImageNF, 5, .025, 5);  %play with parameters if needed
figure; imshow(largeScaleNF); title('Large Scale Image');

%contrast enhancement

logLargeScale = log(largeScaleNF);

maxI = max(largeScaleNF(:))
minI = min(largeScaleNF(:))
range = maxI - minI;
wantedRange = .78;
wantedRange = range;

compressFactor = log(wantedRange)/(log(range))
logOffset = -log(maxI) * compressFactor
logIntensityOut = logLargeScale * compressFactor;
scaledIntensity = 10.^(logIntensityOut);
figure; imshow(scaledIntensity);

%details
detailsNF = intensityImageNF./largeScaleNF;  %intensityImage = largeScale * details
figure; imshow(detailsNF - .5); title('Details');

%% shadow treatment - todo
%% combination of image decomposition

finalIntensityImage = largeScaleNF .* details;
%finalIntensityImage = scaledIntensity .* details;

figure; imshow(finalIntensityImage);

    %max color experiment
    %duoColor = cat(4, YCbCrFlashImage, YCbCrNoFlashImage);
    %maxImage = max(duoColor,[], 4);
    %finalImage = maxImage;.9

%white balance tricks
NoFlashWhiteBalance = mean(mean(noFlashImage,1),2)
FlashWhiteBalance = mean(mean(flashImage,1),2)

FinalWhiteBalance = NoFlashWhiteBalance.^.0001
FinalWhiteBalance = FinalWhiteBalance./(max(FinalWhiteBalance))
FinalWhiteBalance = [1 1 1];  %disable this white balance business for now... 


finalImage = YCbCrFlashImage; 
finalImage(:,:,1) = finalIntensityImage;
finalImageRGB = ycbcr2rgb(finalImage);
finalImageRGBColorCast = cat(3, finalImageRGB(:,:,1).* FinalWhiteBalance(1), finalImageRGB(:,:,2).* FinalWhiteBalance(2), finalImageRGB(:,:,3).* FinalWhiteBalance(3));

figure; imshow(finalImageRGBColorCast);
