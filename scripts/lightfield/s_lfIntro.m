%  Experiments with light fields and ISET
%

%
%  See also: s_s3dISETLF.m, LFToolbox0.4
%
% (AL, BW) Vistasoft Team, 2015

%%
ieInit
clx
ieInit

%% load a lightfield as an oi object

%  To retrieve some of the lightfield data use
% urlBase = 'http://scarlet.stanford.edu/validation/SCIEN/LIGHTFIELD';
% fname = 'indObjLFOiDirect.mat';

urlBase = 'http://scarlet.stanford.edu/validation/SCIEN/LIGHTFIELD/scenes';
fname = 'benchLF.mat';

urlwrite(fullfile(urlBase,fname),fname)
in = load(fname);
oi = in.oi;

% This should be stored in the file.  But that file isn't right, so we set
% parameters here until the data are fixed.
oi = oiSet(oi,'fov',20);
oi = oiSet(oi,'optics focal length',3.5e-3);
oi = oiAdjustIlluminance(oi,10);

% This OI was created using PBRT. The way in which this file was created
% should be up included with the data somewhere, perhaps on the remote data
% path (at scarlet).
% in = load(fullfile(dataPath, 'lightfields', 'indObjLFOiDirect.mat'));
% oi = in.opticalimage;

%% Check the oi and show it.

oiGet(oi,'sample spacing','um')
oiGet(oi,'fov')
oiGet(oi,'mean illuminance')

vcAddObject(oi); oiWindow;

%% process sensor and image processing data
%
% Create a sensor in which each pixel is aligned with a single sample in
% the OI.  Then produce the sensor data (which will include color filters)
%
ss = oiGet(oi,'sample spacing','m');
sensor = sensorCreate;
sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
sensor = sensorSet(sensor,'exp time',0.010);

% Describe
sensorGet(sensor,'pixel size','um')
sensorGet(sensor,'size')
sensorGet(sensor,'fov',[],oi)

%% Compute the sensor response
sensor = sensorCompute(sensor,oi);
vcAddObject(sensor); sensorWindow('scale',1);

%% Interpolate the color filter data to produce a full sensor
%
ip = ipCreate;
ip = ipCompute(ip,sensor);
vcAddObject(ip); ipWindow;

% Show in a separate window
rgb = ipGet(ip,'result');
size(rgb)

%% Pack the samples of the rgb image into the lightfield structure

% If we had a lightfield structure, lf, this could become
%    rgb2lf(rgb,lf)

superPixelH = size(rgb,1)/in.numPinholesH;
superPixelW = size(rgb,2)/in.numPinholesW;

% Allocate space I AM WORRIED ABOUT THE ORDERING OF 3 and 4.  I put in what
% sweems right, but the original code had H and W reversed. (BW).
lightfield = zeros(superPixelH, superPixelW, in.numPinholesH, in.numPinholesW, 3);

% For numerical calculations, we would use this
for i = 1:in.numPinholesW
    for j = 1:in.numPinholesH
        lightfield(:,:, j, i, :) = ...
            rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end

% For visualization, this might be a good idea - use lrgb2srgb
% Allocate space I AM WORRIED ABOUT THE ORDERING OF 3 and 4.  I put in what
% sweems right, but the original code had H and W reversed. (BW).
LF = zeros(superPixelH, superPixelW, in.numPinholesH, in.numPinholesW, 3);
rgb = lrgb2srgb(double(rgb));
for i = 1:in.numPinholesW
    for j = 1:in.numPinholesH
        LF(:,:, j, i, :) = ...
            rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end


%% Some views of the light field data

% lightField(:,:, row, col, :) gives us a view from corresponding pixels in
% each of the pinhole (microlens array) data sets.
%
% The pixels at the edges don't really get any rays or if they do they get
% very little late (are noisier).

vcNewGraphWin;
cnt = 1;
row = superPixelH; col = superPixelW;
rList = 1:2:row;
cList = 1:2:col;
for rr=rList
    for cc=cList
        img = squeeze(lightfield(rr,cc,:,:,:));
        img = lrgb2srgb(img);
        subplot(length(rList),length(cList),cnt), imshow(img);
        cnt = cnt + 1;
    end
end

%% Compare the leftmost and rightmost images in the middle
vcNewGraphWin([],'wide');

img = squeeze(lightfield(3,2,:,:,:));
img = lrgb2srgb(img);
subplot(1,2,1), imshow(img);

img = squeeze(lightfield(3,8,:,:,:));
img = lrgb2srgb(img);
subplot(1,2,2), imshow(img);

%% render some example images

% If we sume all the r,g and b pixels within each aperture we get a single
% RGB image corresponding to the mean.  This is the image at the microlens
% itself
tmp = squeeze(sum(sum(lightfield,2),1));
vcNewGraphWin;
tmp = tmp/max(tmp(:));
imagescRGB(lrgb2srgb(tmp));

%% Now what would the image have been if we move the sensor forward?

vcNewGraphWin
for Slope = -0.5:0.1:0.5
    ShiftImg = LFFiltShiftSum(lightfield, Slope );
    imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
    axis image; truesize
    title(sprintf('Parameter %0.2f',Slope))
    pause(0.1)
end

%% The white image

% What is this fourth plane?  I think it is an overall intensity estimate.
% We need to calculate this for our simulation.  At present, it is just
% arbitrary.
wImage = ShiftImg(:,:,4);
vcNewGraphWin;
imagesc(wImage);

%%  Interact with the lightfield using the toolbox

% In this case, we use the srgb representation because we are just
% visualizing
LFDispMousePan(LF)

%% 
LFDispVidCirc(LF)

%% END
