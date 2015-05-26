
%% Render left pinhole oi
%the pinhole scene can be treated as the scene radiance
%sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinhole.pbrt');
sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinholeWhite.pbrt');
pinholeOiLeft = s3dRenderOIAndDepthMap(sceneName, 'indObj');
vcAddObject(pinholeOiLeft); oiWindow;

%save depthMap as image
depthMap = oiGet(pinholeOiLeft, 'depth map');
imwrite(depthMap./255, 'indObjGTDepth.png');

%% Render right pinhole oi
%the pinhole scene can be treated as the scene radiance
%sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinholeRight.pbrt');
sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinholeRightWhite.pbrt');
pinholeOiRight = s3dRenderOI(sceneName, 'indObj');
vcAddObject(pinholeOiRight); oiWindow;

%% Process left pinhole oi through sensor and image processor

% first set the oi mean illuminance
meanIlluminance = 10; %.5;  %10

% specify which oi you want to use
oi = pinholeOiRight;
%oi = pinholeOiLeft;

    %% create sensor images
    readNoise = 0;
    oi = oiSet(oi,'optics focal length',0.034);  % I set these variables by hand just to avoid the zero condition
    oi = oiSet(oi,'optics f number', 2);    
    oi = oiAdjustIlluminance(oi, meanIlluminance);
    
    %% process sensor and image processing data
    %
    % Create a sensor in which each pixel is aligned with a single sample in
    % the OI.  Then produce the sensor data (which will include color filters)
    %
    ss = oiGet(oi,'sample spacing','m');
    sensor = sensorCreate;
    sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
    sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
    sensor = sensorSet(sensor,'exp time',0.0001);
    %sensor = sensorSet(sensor,'exp time',0.25386);
    %sensor = sensorSet(sensor, 'autoexposure', 1);

    % pixel = sensorGet(sensor, 'pixel')
    %sensor = sensorSet(sensor, 'prnu level', .001);
    sensor = sensorSet(sensor, 'dsnu level', .001);

    % Describe
    sensorGet(sensor,'pixel size','um')
    sensorGet(sensor,'size')
    sensorGet(sensor,'fov',[],oi)

    %% Compute the sensor response
    sensor = sensorCompute(sensor,oi);
    vcAddObject(sensor); sensorWindow('scale',1);

    %% Interpolate the color filter data to produce a full sensor
    %
    ip = ipCreate;
    ip = ipCompute(ip,sensor);
    vcAddObject(ip); ipWindow;

    % Show in a separate window
    %rgb = ipGet(ip,'result');
    
%% code needed to calculate direction vectors
%v1 = [-39.8207  -69.2222   26.4228]
% v2 = [-39.3933  -68.3181   26.2080]
% 
% dirVector= v2-v1
% 
% sideVector = [dirVector(2) -dirVector(1) 0];
% sideVector = sideVector./sum(sideVector .*sideVector)
% 
% v1New = v1 + 5 * sideVector
% v2New = v2 + 5 * sideVector
