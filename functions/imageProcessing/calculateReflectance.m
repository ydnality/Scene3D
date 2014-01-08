%%% surfaceReflectance = calculateReflectance(imageData, bases, illuminant, sensorResponse)

% imageData: the rgb image data in h x w x n format
% bases: bases used for surface reflectance estimation.  nWave x nBases
% illuminat: assumed illuminant of the scene n x 1
% sensorResponse: overall sensor response.  nWave x nBases

% Reflectance is calculated using the same method as the one used in
% DiCarlo et al.  The sensor response function can be written in terms of
% the surface reflectance, the ambient illuminant, and spectral response of
% the sensor.  We then assume that the surface reflectance can be written
% in terms of a weighted sum of basis functions.  The equation is inverted
% to obtain the solution of the weights for the basis functions. 
function surfaceReflectanceCalc = calculateReflectance(imageData, bases, illuminant, sensorResponse)
 
    [imageData,r,c] = RGB2XWFormat(imageData);
    %relative spectral response of illuminant is difficult to determine
    illuminant = diag(illuminant);
    surfaceReflectanceCalcXW = (bases * (((sensorResponse)' * illuminant * bases)^-1 * imageData'))';
    surfaceReflectanceCalc = XW2RGBFormat(surfaceReflectanceCalcXW,r,c);

end