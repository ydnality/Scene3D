% Flash/noflash experiment - basedly loosely off of "Flash Photography
% Enhancement via Instrinsic Relighting" by Eisenmann, and Durand.

%% load images
%load flash image
% fullName = 'benchFlash.mat';
% fullName = 'indObjFlash.mat';
fullName = 'deskFlash.mat';
load(fullName,'vci');
vciFlash = vci;
vcAddAndSelectObject('vcimage',vciFlash);
vcimageWindow;

%load no-flash image
% fullName = 'benchNoFlash.mat';
% fullName = 'indObjNoFlash.mat';
fullName = 'deskNoFlash.mat';
%fullName = 'deskTest.mat';   %testing exposure things
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
figure; imshow(flashImage); title('Original Flash Image');
imwrite(flashImage, 'flashImage.tiff', 'TIFF', 'compression', 'none');

HSVFlashImage = rgb2hsv(flashImage);
intensityImage = HSVFlashImage(:,:,3);
figure; imshow(intensityImage); title('Flash Intensity Image');

%separate into details and large scale

%large scale
%smallImage = imresize(intensityImage, .5);  %put this line back in if you
%want to debug on small images
largeScale = bilateralFilter(intensityImage, 5, .1, 5);  %play with parameters if needed
figure; imshow(largeScale); title('Large Scale Image Flash');

%details
details = intensityImage./largeScale;  %intensityImage = largeScale * details
figure; imshow(details - .5); title('Details Flash');

%verification - commented out for now
% intensityImageRecovery = largeScale .* details;
% figure; imshow(intensityImageRecovery); title('Intensity Recovery');
% sum(sum(intensityImage - intensityImageRecovery))   %checks to see if near 0
%% no flash image processing

noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], imageGet(vciNoFlash, 'gamma'));
% noFlashImage = imadjust(imageGet(vciNoFlash, 'results'), [], [], .6);
figure; imshow(noFlashImage);  title('Original No Flash Image');
imwrite(noFlashImage, 'noFlashImage.tiff', 'TIFF', 'compression', 'none');

HSVNoFlashImage = rgb2hsv(noFlashImage);
intensityImageNF = HSVNoFlashImage(:,:,3);
% figure; imshow(intensityImageNF); title('No Flash Intensity Image');

%separate into details and large scale

%large scale
%smallImage = imresize(intensityImage, .5);  %put this line back in if you
%want to debug on small images

%classic bilateral filter
% largeScaleNF = bilateralFilter(intensityImageNF, 5, .1, 5);  %play with parameters if needed   %works for indObj
largeScaleNF = bilateralFilter(intensityImageNF, 5, .007, 5);  %play with parameters if needed 
figure; imshow(largeScaleNF); title('Large Scale Image (bilateral No Flash)');

%use cross-bilateralFilter
% largeScaleNF = crossBilateralFilter(intensityImageNF, intensityImage, 5, .4, 5);  %play with parameters if needed
% 
% 
% figure; imshow(largeScaleNF); title('Large Scale Image (cross-bilateral)');

%contrast enhancement
logLargeScale = log(largeScaleNF);

maxI = max(largeScaleNF(:))
minI = min(largeScaleNF(:))
range = maxI - minI;
wantedRange = 1;
% wantedRange = range;

compressFactor = log(wantedRange)/(log(range))
logOffset = -log(maxI) * compressFactor
logIntensityOut = logLargeScale * compressFactor;
scaledIntensity = 10.^(logIntensityOut);
% figure; imshow(scaledIntensity);

%details
detailsNF = intensityImageNF./scaledIntensity; 
% detailsNF = intensityImageNF./largeScaleNF;  %intensityImage = largeScale * details
figure; imshow(detailsNF - .5); title('Details No Flash');


%% experiment with no-flash image - using the shadow weighted bilateral filter... 

% meanIntensityNF = mean(intensityImageNF(:));
% meanIntensityF = mean(intensityImage(:));
% 
% normalizedIntensityNF = intensityImageNF/meanIntensityNF;
% normalizedIntensity = intensityImage/meanIntensityF;
% 
% differenceImage = normalizedIntensity - normalizedIntensityNF;
% differenceImage(differenceImage < 0) = 0;
% figure; imagesc(differenceImage);
% 
% largeScaleNFShadow = shadowBilateralFilter(intensityImageNF, differenceImage, 5, .015,.1, 5);  %play with parameters if needed. 
% %The shadow bilateral filter is our new technique that prevents the
% %bilateral filter from smoothing through shadow boundaries in the no-flash
% %image
% 
% figure; imshow(largeScaleNFShadow);  title('large scale shadow retention')

%% shadow treatment - todo
%% combination of image decomposition

finalIntensityImage = largeScaleNF .* details;
%finalIntensityImage = largeScaleNFShadow .* details;   % with shadow correction of ambient image

%finalIntensityImage = scaledIntensity .* details;
% figure; imshow(finalIntensityImage);

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

finalImage = HSVFlashImage; 
finalImage(:,:,3) = finalIntensityImage;
finalImageRGB = hsv2rgb(finalImage);
finalImageRGBColorCast = cat(3, finalImageRGB(:,:,1).* FinalWhiteBalance(1), finalImageRGB(:,:,2).* FinalWhiteBalance(2), finalImageRGB(:,:,3).* FinalWhiteBalance(3));

figure; imshow(finalImageRGB); title('Final Image');
imwrite(finalImageRGB, 'finalImageRGB.tiff', 'TIFF', 'compression', 'none');
% figure; imshow(finalImageRGBColorCast);
%% experiment with saturation - this shows some strange things happening with the Ycbcr colorspace
% testImage = YCbCrNoFlashImage;
% testImage(:,:,1) = YCbCrFlashImage(:,:,1);
% testImage = ycbcr2rgb(testImage);
% figure; imshow(testImage);
% 
% testImage = YCbCrFlashImage;
% testImage(:,:,1) = YCbCrNoFlashImage(:,:,1);
% figure; imshow(testImage)
% testImage = ycbcr2rgb(testImage);
% figure; imshow(testImage);

%% experiment2  with saturation - this shows using the hsv colorspace is much better for some reason... perhaps rounding/precision error!
% hsvFlashImage = rgb2hsv(flashImage);
% figure; imshow(hsvFlashImage);  title('hsv flash image');
% 
% hsvNoFlashImage = rgb2hsv(noFlashImage);
% figure; imshow(hsvNoFlashImage);  title('hsv noflash image');
% 
% 
% testImage = hsvFlashImage;
% testImage(:,:,3) = hsvNoFlashImage(:,:,3);
% testImageRgb = hsv2rgb(testImage);
% figure; imshow(testImageRgb);
