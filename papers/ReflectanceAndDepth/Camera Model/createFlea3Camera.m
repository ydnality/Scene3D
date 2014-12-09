function [sensor, optics] = createFlea3Camera(filterID)

run(sprintf('%s/Camera Model/modelFlea3.m',MsVideoRootPath));

sensor = sensorCreate('Monochrome');
pixel =  sensorGet(sensor,'pixel');

% We set these properties here
pixel = pixelSet(pixel,'sizesamefillfactor',[pixelSize pixelSize]);   
pixel = pixelSet(pixel,'conversiongain', conversiongain);        
pixel = pixelSet(pixel,'voltageswing',voltageSwing);                                             
pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;               
pixel = pixelSet(pixel,'readnoisevolts',readnoise);  
pixel = pixelSet(pixel,'wave',standardWave);
pixel = pixelSet(pixel,'pixelspectralqe',qe); 


% Set these sensor properties
sensor = sensorSet(sensor,'wave',standardWave);
sensor = sensorSet(sensor,'exposuretime',exposureDuration); % commented because we set autoexposure
sensor = sensorSet(sensor,'autoExposure','off'); 
sensor = sensorSet(sensor,'noiseFlag',noiseFlag);
sensor = sensorSet(sensor,'rows',simRows);
sensor = sensorSet(sensor,'cols',simCols);
% sensor = sensorSet(sensor,'rows',rows);
% sensor = sensorSet(sensor,'cols',cols);
sensor = sensorSet(sensor,'dsnulevel',dsnu);  
sensor = sensorSet(sensor,'prnulevel',prnu); 
sensor = sensorSet(sensor,'analogGain',analogGain);     
sensor = sensorSet(sensor,'analogOffset',analogOffset);
sensor = sensorSet(sensor,'quantizationmethod',quantizationMethod);

% Stuff the pixel back into the sensor structure
sensor = sensorSet(sensor,'pixel',pixel);
sensor = pixelCenterFillPD(sensor,fillfactor);
% Then we load the calibration data and attach them to the sensor structure
sensor = sensorSet(sensor,'filterspectra',ones(length(standardWave),1));
sensor = sensorSet(sensor,'infraredfilter',ones(length(standardWave),1));

% Select the bandpass filter.
if ~isempty(filterCenterWaves{filterID})
    fName = fullfile(MsVideoRootPath,'Camera Model',sprintf('Bandpass%inm',filterCenterWaves{filterID}));
    filterData = ieReadColorFilter(standardWave,fName);
    sensor = sensorSet(sensor,'filterspectra',filterData);
end

%% Optics

oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'model','DiffractionLimited');
optics = opticsSet(optics,'off axis method','Skip');
optics = opticsSet(optics,'focallength',focalLength);

end

