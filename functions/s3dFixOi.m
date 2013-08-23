%% oi = s3dFixOi(oi, focalLength)
%
% focalLength: the theoretical focal length of the lens.  This becomes the
% sensor distance.
% Puts in the proper optics and FOV and returns the oi
function oi = s3dFixOi(oi, focalLength)
    oi = oiSet(oi, 'photons', oiGet(oi,'photons') * 10^13);  %some normalization issues
    myOptics = oiGet(oi, 'optics');  %to create proper crop at sensor
    myOptics = opticsSet(myOptics,'focalLength', focalLength);  %distance from lens to sensor - sensor will inherit this value!
    oi = oiSet(oi, 'optics', myOptics);
    sensor = s3dProcessSensor(oi, 0, [], 0); %we pre-compute this to find the sensor FOV
    fov = sensorGet(sensor, 'fov');
    oi = oiSet(oi,'fov', fov);  
end