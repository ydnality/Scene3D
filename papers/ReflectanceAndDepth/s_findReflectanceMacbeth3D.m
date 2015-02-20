close all;
clear variables;
clc;

s_initISET;
%%
nBasis = 12;
testFileName = 'MacbethImage';
forceSelection = 0;
run(sprintf('%s/Camera Model/modelFlea3.m',ReflDepthRootPath));


% Load the illuminant (NO NORMALIZATION!)
fName = sprintf('%s/Illuminant.mat',ReflDepthRootPath);
illuminant = ieReadSpectra(fName,standardWave);
illuminant = Energy2Quanta(standardWave,illuminant);

% Load the reference reflectances
fileName = fullfile(isetRootPath,'data','surfaces','macbethChart');
macbethReflectances = ieReadSpectra(fileName,standardWave);

% Create Macbeth basis functions
[basisFcns, ~, score] = pca(macbethReflectances','Centered',false);


%% Make a prediction of the sensor responsivity on an ideally white,
%  uniform surface.
load(sprintf('%s/Data/02062015_scenes/UniformImage.mat',ReflDepthRootPath));

prediction = deltaL*cameraGain.*((cameraFilters')*illuminant); % + cameraOffset;


% Generate the local gain map

Background = im2double(Img); clear Img; 
h = size(Background,1);
w = size(Background,2);
Background = reshape(Background,[h, w,  nChannels, nFilters]);
% Subtract the offset
Background = max(Background - repmat(shiftdim(cameraOffset',-2),[h,w,1,1]),0);

% Display the data
for i=1:1
    figure;
    for j=1:nChannels
        subplot(4,4,j); imshow(Background(:,:,j,i).^(1/2.2));
        title(sprintf('Filter %i, LED %i',i,j));
    end
end


%% The illumination non-uniformity is best derived from the channel 1 data 
%  This is the channel without any additional low-pass filter in the
%  optical path.

% Apply low pass filtering
kernel = fspecial('gaussian',11,11);
Background = imfilter(Background,kernel,'same');


for f=1:1
    predictionMap = repmat(shiftdim(prediction(f,:),-1),[h,w]);
    Background(:,:,:,f) = Background(:,:,:,f)./predictionMap;
end

gainMap = Background; clear Background; clear predictionMap;


for f=1:1
figure;
for c=1:nChannels
    subplot(4,4,c); imagesc(gainMap(:,:,c,f)); colorbar
end
end



%% Extract data from a Macbeth image

load(sprintf('%s/Data/02062015_scenes/%s.mat',ReflDepthRootPath,testFileName));
Img = im2double(Img);
Img = reshape(Img,[h, w,  nChannels, nFilters]);

% Display images

for f=1:1
figure; 
for i=1:nChannels
    subplot(4,4,i); imagesc(Img(:,:,i,f).^(1/2.2)); colorbar;
    title(sprintf('Filter %i, LED %i',f,i));
end
end


% Remove the camera offsets
Img = max(Img - repmat(shiftdim(cameraOffset',-2),[h,w,1,1]),0);

% Correct for the illuminant non-uniformity
Img = Img./gainMap;

measVals = zeros(nFilters,nChannels,24);
cornerPoints = [1    98
   128    96
   128     9
     1     9];

for f=1:1
    
    sensor = createFlea3Camera(f);
    
    for c=1:nChannels
   
        sensor = sensorSet(sensor,'voltage',Img(:,:,c,f)*(pixelGet(sensorGet(sensor,'pixel'),'voltageswing')));
        sensor = sensorSet(sensor,'dv',Img(:,:,c,f)*(2^sensorGet(sensor,'nbits')));
        vcAddObject(sensor);
        
        [vals] = macbethSelect(sensor,[],1,cornerPoints);
        measVals(f,c,:) = cellfun(@nanmean,vals);
        
    end
end

measVals = measVals/(2^sensorGet(sensor,'nbits'));


%% Use only the monochromatic channel
channelIDs = [1 2 3 4 5 6 7 8 13 14];
filterIDs = 1;

lambdaSet = linspace(0,1,1001);

[optLambda, rmsError, est, ~, ~, pred] = xValidateLambdaV2(measVals(filterIDs,channelIDs,:),...
    cameraFilters(:,filterIDs),...
    cameraGain(filterIDs,channelIDs)*deltaL,...
    zeros(length(filterIDs),length(channelIDs)),...
    illuminant(:,channelIDs),...
    basisFcns(:,1:nBasis),...
    lambdaSet,...
    macbethReflectances);

figure; 
plot(lambdaSet,rmsError);
xlabel('Lambda');
ylabel('RMS error');
title(sprintf('Optimal lambda %f',optLambda));

figure;
grid on; box on; hold on;
plot(est,macbethReflectances,'x');
plot(est,est,'r');
xlabel('Estimated');
ylabel('True');
title('Prediction accuracy');
xlim([0 1]);
ylim([0 1]);

figure;
for tx=1:6
    for ty=1:4
        
        indx = ty + (tx-1)*4;
        subfigIndx = (ty-1)*6 + tx;
        
        subplot(4,6,subfigIndx);
        hold on; box on; grid on;
        plot(standardWave,est(:,indx),'rx');
        plot(standardWave,macbethReflectances(:,indx));
        title(sprintf('RMS error %f\n',rms(est(:,indx) - macbethReflectances(:,indx))));
        ylim([0 1]);
        xlim([min(standardWave), max(standardWave)]);
        
    end
end

tmp = measVals(filterIDs,channelIDs,:);

figure;
hold on; grid on; box on; axis equal;
xlim([0 1]);
ylim([0 1]);
plot(tmp(:),pred(:),'x');
xlabel('True');
ylabel('Estimated');
title(sprintf('RMS %f',rms(tmp(:) - pred(:))));


%% Use all channels
%{
channelIDs = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
filterIDs = [1 2 3 4 5 6 7 8];
lambdaSet = linspace(0,1,1000);

[optLambda, meas2basisCoeffs, est, lambdaSet, rmsError, pred] = xValidateLambdaV2(measVals(filterIDs,channelIDs,:),...
                                                                                cameraFilters(:,filterIDs),...
                                                                                cameraGain(filterIDs,channelIDs)*deltaL,...
                                                                                cameraOffset(filterIDs,channelIDs),...
                                                                                illuminant(:,channelIDs),...
                                                                                basisFcns(:,1:nBasis),...
                                                                                lambdaSet,...
                                                                                macbethReflectances);

figure;
plot(lambdaSet,rmsError);
xlabel('Lambda');
ylabel('RMS error');
title('Smoothnes parameter x-val rms error.');

figure;
grid on; box on; hold on;
plot(est,macbethReflectances,'x');
plot(est,est,'r');
xlabel('Estimated');
ylabel('True');
title('Prediction accuracy');
xlim([0 1]);
ylim([0 1]);

figure;
for tx=1:6
    for ty=1:4
        
        indx = ty + (tx-1)*4;
    
        subfigIndx = (ty-1)*4 + tx;
        
        subplot(4,6,subfigIndx);
        hold on; box on; grid on;
        plot(standardWave,est(:,indx),'rx');
        plot(standardWave,macbethReflectances(:,indx));
        title(sprintf('RMS error %f\n',rms(est(:,indx) - macbethReflectances(:,indx))));
        ylim([0 1]);
        xlim([min(standardWave), max(standardWave)]);
        
    end
end

tmp = measVals(:,filterIDs,channelIDs);


figure;
hold on; grid on; box on; axis equal;
xlim([0 1]);
ylim([0 1]);
plot(tmp(:),pred(:),'x');
xlabel('True');
ylabel('Estimated');
title(sprintf('RMS %f',rms(tmp(:) - pred(:))));



for i=1:8
    figure;
    hold on; grid on; box on;
    plot(0:0.1:1,0:0.1:1,'r');
    xlabel('True');
    ylabel('Estimated');
    title(sprintf('Filter %i',i));
    plot(squeeze(measVals(:,i,:)),squeeze(pred(:,i,:)),'x');
end


for i=1:14
    figure;
    hold on; grid on; box on;
    plot(squeeze(measVals(:,:,i)),squeeze(pred(:,:,i)),'x');
    xlabel('True');
    ylabel('Estimated');
    title(sprintf('LED %i',i));
end
%}
