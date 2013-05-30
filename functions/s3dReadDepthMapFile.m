function A3 = s3dReadDepthMapFile(depthMapFile, imageSize)
% Reads in a zbf file with binary depth data
% 
%  A3 = s3dReaddepthMapFile(depthMapFile)
% 
% depthMapFile is path to depth map file, for example, 'depthmap.zbf'
% How are these zbf files created?
%
% Andy's Response: zbf files are created within the RenderToolBox pipeline.
%  I added an additional command in the pipeline that tells Radiance to
%  render a zbf file, which is a binary file containing depth map
%  information.  
% 
% Inputs: 
% depthMapFile: string of the .zbf depth map file produced by Radiance.
% 
% Return values:
%   A3: This function returns the UNNORMALIZED 2D depth map image. 
%   This function also writes, depthMapView.tif, the resulting NORMALIZED depth map 
%   preview image. For the normalized depth map preview image, the value 0 corresponds to 
%   the lower limit and the value 1 corresponds to the upper limit of the unnormalized range.  
% 
% Example: 
%  depthMap = s3dReaddepthMapFile('depthmap.zbf');
% 
% (c) Stanford VISTA Team

%%

% There should be error checking here.



% fid=fopen('depthmap.zbf', 'r', 'l');
fid = fopen(depthMapFile, 'r', 'l');

%A = fread(fid, 'float32');
A = fread(fid, 'double');

% if size is not defined, assume a square
if (ieNotDefined('imageSize'))
    imageSize = [round(sqrt(size(A,1))) round(sqrt(size(A,1)))]
end


% Please fix ...
% assumes a square image, need to change to proportion of image dimensions
% if not a square image

% test = round(sqrt(size(A,1)));

% Why is there garbage?  Please put in comments.
%condition the data, throw away garbage
%A2 = A .* (abs(A) < 1000);

% A2 = A;
% A2 = A2 (1:test^2);

%new stuff
A2 = A(1:imageSize(1)*imageSize(2));
%reshape data
A3 = reshape(A2, [imageSize(2) imageSize(1)])';

%make the range shown appropriate - only for debugging purposes
% lowerLimit = 900;
% upperLimit = 1000;

% This used to print out the lowerLimit and upperLimit
% Either they should be returned, or we should calculate them quietly.
% viewImage = abs(A3 - 900)/100;
% viewImage = abs(A3 - lowerLimit)/(upperLimit-lowerLimit);
%lowerLimit = min(A2);    %will this cause scaling problems later with inconsistent scaling?  No - this is only for viewing purposes!  
%upperLimit = max(A2);

% Scale the data to between 0 and 1
%viewImage = abs((A3 - lowerLimit)/(upperLimit-lowerLimit));
%figure; imshow(viewImage);
% 
% % Write out a tiff file.  This provides a basic debugging output to see if
% % the depthmap was properly read.  
% imwrite(viewImage, 'depthMapView.tif');

return

