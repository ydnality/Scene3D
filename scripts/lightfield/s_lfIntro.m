% s_lfIntro
% 
% Introduction to light field camera simulation using ISET
%
% See also: s_s3dISETLF.m, LFToolbox0.4
%
% (AL, BW) Vistasoft Team, 2015

% Set up ISET
ieInit

%% load a lightfield as an oi object from the remote data site

% Another light field we calculated.
% fname = 'metronomeLF.mat';
fname = 'benchLF.mat';
rd  = ieRdata('create');
val = ieRdata('load data',rd,fname);

%% Unpack the data 
oi  = val.oi;
oi.lightfield.pinholes = [val.numPinholesH, val.numPinholesW];
oi.lightfield.pbrt     = val.curPbrt;

% We set parameters here.
oi = oiSet(oi,'fov',20);
oi = oiSet(oi,'optics focal length',3.5e-3);
oi = oiAdjustIlluminance(oi,10);

% Show the oi as an image
ieAddObject(oi); oiWindow;

% Show a table of the oi parameters
iePTable(oi);

%% Create a sensor 
%
% Each pixel is aligned with a single sample in the OI.  Then produce the
% sensor data (which will include color filters)
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

sensor = sensorCompute(sensor,oi);
ieAddObject(sensor); sensorWindow('scale',1);

%% Interpolate (bilinear) the color filter data to produce an image
%
ip = ipCreate;
ip = ipCompute(ip,sensor);
vcAddObject(ip); ipWindow;

% Show in a separate window
rgb = ipGet(ip,'result');

%% Pack the samples of the rgb image into the lightfield structure
% This is the structure used by the light field toolbox

% The calculations in this section could become
%    rgb2lf(rgb,lf)

superPixelH = size(rgb,1)/val.numPinholesH;
superPixelW = size(rgb,2)/val.numPinholesW;
lightfield = zeros(superPixelH, superPixelW, ...
    val.numPinholesH, val.numPinholesW, 3);

% For lightfield calculations, we use this rgb format
for i = 1:val.numPinholesW
    for j = 1:val.numPinholesH
        lightfield(:,:, j, i, :) = ...
            rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end

% For visualization we display the rgb as srgb 
% This makes it easier to see.
% So, LF is like lightfield but with srgb values
LF = zeros(superPixelH, superPixelW, val.numPinholesH, val.numPinholesW, 3);
rgb = lrgb2srgb(double(rgb));
for i = 1:val.numPinholesW
    for j = 1:val.numPinholesH
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
        subplot(length(rList),length(cList),cnt), imagescRGB(img);
        cnt = cnt + 1;
    end
end

%% Compare the leftmost and rightmost images in the middle
vcNewGraphWin([],'wide');

img = squeeze(lightfield(3,2,:,:,:));
img = lrgb2srgb(img);
subplot(1,2,1), imagescRGB(img);

img = squeeze(lightfield(3,8,:,:,:));
img = lrgb2srgb(img);
subplot(1,2,2), imagescRGB(img);

%% render some example images

% If we sume all the r,g and b pixels within each aperture we get a single
% RGB image corresponding to the mean.  This is the image at the microlens
% itself
tmp = squeeze(sum(sum(lightfield,2),1));
tmp = tmp/max(tmp(:));

vcNewGraphWin;
imagescRGB(lrgb2srgb(tmp));
title('Image at microlens plane')

%% Now what would the image have been if we move the sensor forward?

vcNewGraphWin([],'wide');

% Use different Slopes for the benchLF and metronomeLF pictures
if strcmp(fname,'metronomeLF.mat'),     Slope = -0.5:0.2:1.0;
elseif strcmp(fname,'benchLF.mat'),     Slope = -0.5:0.5:3.0;
end

for ii = 1:length(Slope)
    ShiftImg = LFFiltShiftSum(lightfield, Slope(ii) );
    subplot(1,length(Slope),ii);
    imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
    axis image;
    title(sprintf('Parameter %0.2f',Slope(ii)))
end

%% The white image

% What is this fourth plane?  I think it is an overall intensity estimate.
% We need to calculate this for our simulation.  At present, it is just
% arbitrary.
% wImage = ShiftImg(:,:,4);
% vcNewGraphWin;
% imagesc(wImage);

%%  Interact with the lightfield using the toolbox

% In this case, we use the srgb representation because we are just
% visualizing
% LFDispMousePan(LF)

%% 
LFDispVidCirc(LF)

%% END
