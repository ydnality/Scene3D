%% sensor = s3dProcessSensor(oi, readNoise, size, exposureTime, bitDepth)
% 
% Processes the sensor using assigned readNoise, resolution, exposureTime,
% and bitDepth.  
%
% size: 1 x 2 int specifying sensor size. For example, [400 600]
% exposureTime: exposure duration in seconds.
% bitDepth: the bit depth of the sensor.  For float, use 'analog.'
% Otherwise, '8 bit', '12 bit' etc. are possible options. '12 bit' is the
% max.
function sensor = s3dProcessSensor(oi, readNoise, size, exposureTime, bitDepth)

    if (ieNotDefined('readNoise'))
       readNoise = .00096; %default low noise 
    end
    
    if (ieNotDefined('size'))
        size = [400 600];
    end
    
    if (ieNotDefined('exposureTime'))
       exposureTime = -1; 
    end
    
    
    if (ieNotDefined('bitDepth'))
       bitDepth = 'analog'; 
    end
    
    % oi = oiSet (oi, 'horizontalfieldofview', 8 * 200/150 );
    vcAddAndSelectObject(oi);
    oiWindow;
    m = oiGet(oi, 'mean illuminance')

    % sensor processing - may want to consider putting this into an external
    % helper function
    sensor = sensorCreate('bayer (gbrg)');
    sensor = sensorSet(sensor, 'size', size);   %*make this a parameter
    
    %using desired sensor size, set pixel size accordingly
    sensorWidth = 36;  %*make this a parameter later
    sensorHeight = 24;
    pixelSize = sensorWidth * .001/size(2);
    
    % We set the sensor properties using sensorSet and sensorGet routines.
    %
    % Just as the optical irradiance gives a special status to the optics, the
    % sensor gives a special status to the pixel.  In this section we define
    % the key pixel and sensor properties, and we then put the sensor and pixel
    % back together.

    % To get the pixel structure from the sensor we use:
    pixel =  sensorGet(sensor,'pixel');

    % Here are some of the key pixel properties
    voltageSwing   = 1.15;  % Volts
    wellCapacity   = 90000;  % Electrons  %flash image
    conversiongain = voltageSwing/wellCapacity;   
    fillfactor     = 0.45;       % A fraction of the pixel area
    %pixelSize      = 2.2*1e-6;   % Meters
    darkvoltage = 0;
    %darkvoltage    = 1e-005;     % Volts/sec
    readnoise      = readNoise;    % Volts   

    % We set these properties here
    pixel = pixelSet(pixel,'size',[pixelSize pixelSize]);   
    pixel = pixelSet(pixel,'conversiongain', conversiongain);        
    pixel = pixelSet(pixel,'voltageswing',voltageSwing);                                             
    pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;               
    pixel = pixelSet(pixel,'readnoisevolts',readnoise);  

    %  Now we set some general sensor properties
    %     exposureDuration = 03.0; % commented because we set autoexposure
    %     dsnu =  0.0010;           % Volts (dark signal non-uniformity)     no flash image
    dsnu = 0;    %flash image
    %prnu = 0; 
    prnu = 0.2218;            % Percent (ranging between 0 and 100) photodetector response non-uniformity
    %prnu = 0.4218; 
    analogGain   = 1;         % Used to adjust ISO speed
    analogOffset = 0;         % Used to account for sensor black level

    %      autoExpTime = autoExposure(oi,sensor, 1, 'default')
     if (exposureTime <= 0)
         sensor = sensorSet(sensor,'auto exposure', 1);
     else
        sensor = sensorSet(sensor,'exptime', exposureTime);
     end

    sensor = sensorSet(sensor,'dsnulevel',dsnu);  
    sensor = sensorSet(sensor,'prnulevel',prnu); 
    sensor = sensorSet(sensor,'analogGain',analogGain);     
    sensor = sensorSet(sensor,'analogOffset',analogOffset);   

    % Stuff the pixel back into the sensor structure
    sensor = sensorSet(sensor,'pixel',pixel);
    
    sensor = pixelCenterFillPD(sensor,fillfactor);

    % It is also possible to replace the spectral quantum efficiency curves of
    % the sensor with those from a calibrated camera.  We include the
    % calibration data from a very nice Nikon D100 camera as part of ISET.
    % To load those data we first determine the wavelength samples for this sensor.
    wave = sensorGet(sensor,'wave');

    % Then we load the calibration data and attach them to the sensor structure
    fullFileName = fullfile(isetRootPath,'data','sensor','colorfilters','NikonD100.mat');
    [data,filterNames] = ieReadColorFilter(wave,fullFileName); 
    sensor = sensorSet(sensor,'filterspectra',data);
    sensor = sensorSet(sensor,'filternames',filterNames);
    sensor = sensorSet(sensor, 'quantization', bitDepth);
    
    % We are now ready to compute the sensor image
    sensor = sensorCompute(sensor,oi);

    % Add name to correspond to oi name
    sensor = sensorSet(sensor, 'name', oiGet(oi, 'name'));
    % We can view sensor image in the GUI.  Note that the image that comes up
    % shows the color of each pixel in the sensor mosaic. Also, please be aware
    % that the Matlab rendering algorithm often introduces unwanted artifacts
    % into the display window.  You can resize the window to eliminate these.
    % You can also set the display gamma function to brighten the appearance in
    % the edit box at the lower left of the window.
end