close all;
clear variables;
clc;

%% Define script parameters

run(sprintf('%s/Camera Model/modelFlea3.m',ReflDepthRootPath));
macbethReflectances = ieReadSpectra(fullfile(isetRootPath,'data','surfaces','macbethChart'),standardWave);

illuminantEnergy = ieReadSpectra(sprintf('%s/Illuminant.mat',ReflDepthRootPath),standardWave);
illuminant = Energy2Quanta(standardWave,illuminantEnergy);

%% Run simulations
Img = zeros([simRows, simCols, nChannels*nFilters]);

cameraGain = zeros(nFilters,nChannels);
cameraOffset = zeros(nFilters,nChannels);
cameraExposure = zeros(nFilters,nChannels);

cameraPixelVals = zeros(nFilters,nChannels,24);
modelPixelVals = zeros(nFilters,nChannels,24);

cornerPoints = [1    98
   128    96
   128     9
     1     9];

 
for f=1:1
    tic
    [sensor,optics] = createFlea3Camera(f);
    
    for ch = 1:nChannels
        
        fprintf('Simulating filter %i channel %i\n',f,ch);
        
        load(sprintf('%s/Data/1208Scenes/macbeth%i.mat',ReflDepthRootPath,ch));
        depthMap = sceneGet(scene,'depthMap');
        scene = sceneSet(scene,'photons',sceneGet(scene,'photons')*1e6);
        centerDist = depthMap(round(size(depthMap,1)/2),round(size(depthMap,2)/2));
        centerDist = centerDist / 1000;
        scene = sceneSet(scene,'distance',centerDist);
        scene = sceneSet(scene,'fov',0.08);
        vcAddObject(scene);
        
        
        %{
        spectrum.wave = standardWave;
        scene = sceneCreate('macbethd65',20,spectrum);
        scene = sceneAdjustIlluminant(scene,7.66*illuminantEnergy(:,ch),0);
        scene = sceneSet(scene,'distance',1);
        scene = sceneSet(scene,'fov',0.48);
        scene = sceneSet(scene,'name',sprintf('Filter %i, channel %i',f,ch));
        vcAddObject(scene);
        %}
        
        % Compute the optical image
        oi = oiCreate();
        oi = oiSet(oi,'optics',optics);
        oi = oiSet(oi,'Name',sprintf('Filter %i, channel %i',f,ch));
        oi = oiCompute(scene,oi);
        vcAddObject(oi);
        
        % Compute the sensor image
        cameraExposure(f,ch) = autoExposure(oi,sensor,0.95,'luminance');
        sensor = sensorSet(sensor,'exposureTime',cameraExposure(f,ch));
        sensor = sensorSet(sensor,'Name',sprintf('Filter %i, channel %i',f,ch));
        sensor = sensorCompute(sensor,oi);
        vcAddObject(sensor);
        
        Img(:,:,(f-1)*nChannels + ch) = sensorGet(sensor,'dv')/(2^sensorGet(sensor,'nbits'));
        [cameraGain(f,ch), cameraOffset(f,ch)] = getSensorTotalGain(scene,oi,sensor);
        
        
        [vals, ~, ~, ~] = macbethSelect(sensor,[],1,cornerPoints);
        cameraPixelVals(f,ch,:) = cellfun(@nanmean,vals);
        modelPixelVals(f,ch,:) = cameraGain(f,ch)*cameraFilters(:,f)'*diag(illuminant(:,ch))*macbethReflectances*deltaL;
        
    end
    toc
end

cameraPixelVals = cameraPixelVals/255;

for i=1:1
    figure;
    hold on; grid on; box on;
    for j=1:nChannels
    plot(squeeze(modelPixelVals(i,j,:))',squeeze(cameraPixelVals(i,j,:))','.');
    pause
    end
    title(sprintf('Filter %i',i))
    ylabel('ISET model values');
    xlabel('Linear model values');
end

figure;
hold on; grid on; box on;
for ch=1:nChannels
    A = [squeeze(modelPixelVals(1,ch,:))];
    b = squeeze(cameraPixelVals(1,ch,:)) - cameraOffset(1,ch);
    
    coeff = A\b;
    % cameraOffset(1,ch) = coeff(1)*cameraOffset(1,ch) + coeff(2);
    % cameraGain(1,ch) = cameraGain(1,ch)*coeff(1);
    fprintf('Channel %i, gain %f\n',ch,coeff(1));
    plot(b,A*coeff,'.');
    
end

if ~exist(sprintf('%s/Data/',ReflDepthRootPath),'dir')
   mkdir(sprintf('%s/Data/',ReflDepthRootPath));
end
save(sprintf('%s/Data/MacbethImage.mat',ReflDepthRootPath),'Img','cameraGain','cameraOffset','cameraExposure');

%%
sceneWindow;
oiWindow;
sensorWindow;
