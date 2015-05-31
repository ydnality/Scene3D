%% Render the pinhole scene and Depth Map

%the pinhole scene can be treated as the scene radiance
sceneName = fullfile(s3dRootPath, 'papers', '2014-OSA', 'indestructibleObject', 'mainPinhole.pbrt');
pinholeOi = s3dRenderOIAndDepthMap(sceneName, 'indObj');

%%  Create a sensor in which each pixel is aligned with a single sample in
% the OI.  Then produce the sensor data (which will include color filters)
%
ss = oiGet(pinholeOi,'sample spacing','m');
sensor = sensorCreate;
sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
sensor = sensorSet(sensor,'exp time',0.0005);
%sensor = sensorSet(sensor,'exp time',0.25386);
%sensor = sensorSet(sensor, 'autoexposure', 1);

% pixel = sensorGet(sensor, 'pixel')
%sensor = sensorSet(sensor, 'prnu level', .001);
sensor = sensorSet(sensor, 'dsnu level', .001);

% Describe
sensorGet(sensor,'pixel size','um')
sensorGet(sensor,'size')
sensorGet(sensor,'fov',[],pinholeOi)

%
vcReplaceObject(pinholeOi); oiWindow;
oiGet(pinholeOi,'fov')

%%

sensor = sensorCompute(sensor,pinholeOi);
vcAddObject(sensor); sensorWindow('scale',1);

%% Interpolate the color filter data to produce a full sensor
%
ip = ipCreate;
ip = ipCompute(ip,sensor);
% vcAddObject(ip); ipWindow;

% Show in a separate window
rgb = ipGet(ip,'result');