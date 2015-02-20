%% s_s3dChessboard
%
% Andy's lightfield calculation idea
%
%  We will create a light field image from a computer graphics scene.
%  We will create the image using PBRT.
%  The output data will get put into the LF format below.
%  We will render it using the LF toolbox utilities
%
%  Process:
%    Create an ordinary camera model in PBRT.
%
%    Add an additional small aperture at some position in front of the
%    sensor. The sensor samples are (s1 , s2).  For example, (17,17). These
%    are behind the pinhole.  The pin holes will be at positions (p1,p2)
%
%    Run PBRT to calculate the intensity at each wavelength at each of the
%    (s1,s2) samples. Put the samples in the (1,1) position of each of the
%    (p1,p2) pinholes into the LF structure as follows.
%
%    Suppose wavelength = 1.  
%    Then the (s1,s2) sample from the (p1,p2) pinhole should be placed into
%    the light field entry LF(s1,s2,p1,p2,w)
%
%
% This is a contributed LF toolbox
% http://www.mathworks.com/matlabcentral/fileexchange/49683-light-field-toolbox-v0-4
%
% A few demos here.  We can add more examples.

%addpath(genpath((fullfile(isetRootPath,'../external','LFToolbox0.4'))))

%% Example data set
cd(fullfile(dataPath,'lightfields'))

% Read the light field data
lightF = LFReadGantryArray('Chess/rectified',struct('UVLimit',256));

% Here is one picture
% img = squeeze(LF(10,10,:,:,2)); 
% vcNewGraphWin; imagesc(img); colormap(gray)

%%
LFDispMousePan(lightF)

%%
LFDispVidCirc(lightF)


