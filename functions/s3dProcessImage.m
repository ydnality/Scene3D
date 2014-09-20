function [vci, transformMatrices ] = s3dProcessImage(sensor, wantedTransformationMatrices)

    if (ieNotDefined('wantedTransformationMatrices'))
        disp ('transformation Matrix not supplied');
        %default transforms (identity for all assumed)
        wantedTransformationMatrices = {eye(3), eye(3), eye(3)};
    end
    
    vci = vcimageCreate;

    % The routines for setting and getting image processing parameters are
    % imageGet and imageSet.
    %
    vci = imageSet(vci,'name','Unbalanced');
    vci = imageSet(vci,'scaledisplay',1);
    vci = imageSet(vci,'renderGamma',0.6);

    % The default properties use bilinear demosaicking, no color conversion or
    % balancing.  The sensor RGB values are simply set to the display RGB
    % values.
    
    vci = imageSet(vci,'illuminant correction method','None');
    vci = imageSet(vci,'illuminant correction matrix','None');
    vci = imageSet(vci,'internal color space','Sensor');
    vci = imageSet(vci, 'sensorconversionmethod', 'None');

%     vci = vcimageCompute(vci,sensor);;

    vci = imageSet(vci, 'transforms', wantedTransformationMatrices);  
    
    vci = vcimageCompute(vci,sensor);
    vci = imageSet(vci, 'transforms', wantedTransformationMatrices);
    
    
    vci = vcimageCompute(vci,sensor);
    transformMatrices = imageGet(vci, 'transforms');
    
    % Add name to correspond to same name as sensor
    vci = imageSet(vci, 'name', sensorGet(sensor, 'name'));
    
    % As in the other cases, we can bring up a window to view the processed
    % data, this time a full RGB image.
    
%     vcAddAndSelectObject(vci); vcimageWindow

%     % You can experiment by changing the processing parameters in many ways,
%     % such as:
%     vci2 = imageSet(vci,'name','More Balanced');
%     vci2 = imageSet(vci2,'internalCS','XYZ');
%     vci2 = imageSet(vci2,'colorConversionMethod','MCC Optimized');
%     vci2 = imageSet(vci2,'colorBalanceMethod','Gray World');
% 
%     % With these parameters, the colors will appear to be more accurate
%      image = vcimageCompute(vci,sensor);

end