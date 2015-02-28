%  Experiments with light fields and ISET
%
% In this case, we start with data Andy produced to perform various
% manipulations.
%
% (BW) Vistasoft Team, 2015

% The superpixel scene files are in:
% /home/ydna/Scene3D/data/pbrtScenes/benchScene/LF/ 
% They are named superpixel'i'_'j'.mat
% 
% On black, I made a link to Scene3d/data/LF
%
% So, we need to be able to read these in and convert them into one
% irradiance matrix suitable for an OI.
%
% There are 6400 superpixels (80 x 80).
%   dDir = fullfile(s3dRootPath,'data','LF');
%   files = dir(fullfile(dDir,'superpixel*.mat'));
%
% The 5D lightfield (i,j, rows, col, wavelengths) is located here:
% /home/ydna/Scene3D/data/pbrtScenes/benchScene/LF/LF.mat
% 
load(fullfile(dDir,'LF.mat'),'lightField');

%% Some views
% lightField(:,:, row, col, :) gives us a pinhole view at a specific
% position on the aperture. 
%  Notice that the pixels at the edges don't really get any rays or if they
%  do they get very little late (are noisier).
vcNewGraphWin;
cnt = 1;
row = 9; col = 9;
for rr=1:row
    for cc=1:col
        img = squeeze(lightField(:,:,rr,cc,:));
        img = imageTranspose(img);
        subplot(9,9,cnt), imageSPD(img,400:10:700);
        cnt = cnt + 1;
    end
end


%%

% lightField(i,j, :, :, :) gives us the view of
% the superpixel at position (i,j).
% 
% 
% The last sections of s_s3dRenderHDrBenchLF.m has some examples where I
% manipulate the lightField to create some simple cool examples.
%