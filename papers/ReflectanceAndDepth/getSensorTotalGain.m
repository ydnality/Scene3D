function [ISETGain, ISETOffset] = getSensorTotalGain(scene,oi,sensor)

% [ISETGain, ISETOffset] = getSensorTotalGain(scene,oi,sensor)
%
% This function takes into account all the scene, optics and sensor
% parameters to compute the total gain and offest applied to a camera
% pixel.
%
% In a linear image formation model the pixel intensity is given by
%
%    m = ISETGain*sum(reflectance.*illuminant.*sensitivity) + ISETOffset
%
%
%   Copuright, Henryk Blasinski 2014


%% GAIN

% 1. The first part of the gain has to do with the radiance to irradiance conversion
sDist = sceneGet(scene,'distance');
fN    = opticsGet(oiGet(oi,'optics'),'fNumber');     
m     = opticsGet(oiGet(oi,'optics'),'magnification',sDist);
rad2irr = pi /(1 + 4*fN^2*(1+abs(m))^2);


% 2. The second part of the gain is due to scene surface area to pixel area
% mappings
q = vcConstants('q');
gridSpacing = 1/sensorGet(sensor,'nSamplesPerPixel');
gridSpacing = 1/round(1/gridSpacing);
nGridSamples = 1/gridSpacing;

if nGridSamples == 1, pdArray = pixelGet(sensorGet(sensor,'pixel'),'fillfactor');
else                  pdArray = sensorPDArray(sensor,gridSpacing);
end

area = pixelGet(sensorGet(sensor,'pixel'),'area')/(nGridSamples^2);


% 3. Finally we have to take into account the electrn to voltage
% conversions
pixel = sensorGet(sensor,'pixel');
cur2volt = sensorGet(sensor,'integrationTime')*pixelGet(pixel,'conversionGain') / q;
volt2dv = 1/pixelGet(pixel,'voltageswing');
analogGain = sensorGet(sensor,'analoggain');

% We combine the effects of 1, 2 and 3 into ISETGain
ISETGain = q*volt2dv*cur2volt*pdArray*area*rad2irr/analogGain;



%% OFFSET
% The offset parameter is due to an analog offest introduced on purpose and
% the dark voltage that is a function of the integration time.

ISETOffset = sensorGet(sensor,'analogoffset')/sensorGet(sensor,'analoggain') + pixelGet(pixel,'dark Voltage')*sensorGet(sensor,'integrationTime');

end