%  Experiments with light fields and ISET for the slanted bar
%
%  See also: s_s3dISETLF.m, LFToolbox0.4
%
%  To retrieve some of the lightfield data use
%     urlBase = 'http://scarlet.stanford.edu/validation/SCIEN/LIGHTFIELD/scenes';
%     fname = 'metronomeLF.mat';
%     fname = 'benchLF.mat';
%     urlwrite(fullfile(urlBase,fname),fname)
%
% This will work some day (sort of works now, but more checking needed)
%    s3dGetRdata(fullfile(urlBase,fname),fname);
%
% (AL) Vistasoft Team, 2015

%%
ieInit

%% inputs

% load a lightfield as an oi object.
% These should be available on the scarlet/validation web site.
%in = load('benchLFHQ.mat');
in = load('slantedBarMultLF160Diffract.mat');
% in = load('metronomeLF.mat');

oi = in.oi;

numPinholesW = in.numPinholesW
numPinholesH = in.numPinholesH
curPbrt = in.curPbrt

% The optics should be reasonable.  We set a focal length of 35mm here.
%oi = oiSet(oi,'optics focal length',3.5e-2);
oi = oiSet(oi,'optics focal length',6.8e-2);

% We say the field of view is 10 deg.  This produces a reasonable sample
% spacing, corresponding to about one pretty big pixel.
oi = oiSet(oi,'fov',20);

% Set the mean illuminance in lux to 10.  Reasonably bright (bright room)
oi = oiAdjustIlluminance(oi,10);

vcAddObject(oi); oiWindow;

% Description - Just checking
oiGet(oi,'sample spacing','um')
oiGet(oi,'fov')
oiGet(oi,'mean illuminance')

%convert to RGB (we are skipping the sensor processing step for simplicity.
%That will come later)
% rgb = oiGet(oi, 'rgb');
% vcNewGraphWin;
% imagescRGB(rgb);

%% process sensor and image processing data
%
% Create a sensor in which each pixel is aligned with a single sample in
% the OI.  Then produce the sensor data (which will include color filters)
%
ss = oiGet(oi,'sample spacing','m');
sensor = sensorCreate;
sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
sensor = sensorSet(sensor,'exp time',0.0005);
%sensor = sensorSet(sensor,'exp time',0.25386);
%sensor = sensorSet(sensor, 'autoexposure', 1);

% pixel = sensorGet(sensor, 'pixel')
%sensor = sensorSet(sensor, 'prnu level', .001);
sensor = sensorSet(sensor, 'dsnu level', .001);

% Describe
sensorGet(sensor,'pixel size','um')
sensorGet(sensor,'size')
sensorGet(sensor,'fov',[],oi)

%
vcReplaceObject(oi); oiWindow;
oiGet(oi,'fov')


%% create a larger sensor to see the difference from diffraction

ss = oiGet(oi,'sample spacing','m');
sensor = sensorCreate;
sensor = sensorSet(sensor,'pixel size same fill factor',ss(1)/2);
sensor = sensorSet(sensor,'size',oiGet(oi,'size') * 2);
sensor = sensorSet(sensor,'exp time',0.0005/4);
%sensor = sensorSet(sensor,'exp time',0.25386);
%sensor = sensorSet(sensor, 'autoexposure', 1);

% pixel = sensorGet(sensor, 'pixel')
%sensor = sensorSet(sensor, 'prnu level', .001);
sensor = sensorSet(sensor, 'dsnu level', .001);

% Describe
sensorGet(sensor,'pixel size','um')
sensorGet(sensor,'size')
sensorGet(sensor,'fov',[],oi)

%
vcReplaceObject(oi); oiWindow;
oiGet(oi,'fov')
%% Compute the sensor response

sensor = sensorCompute(sensor,oi);
vcAddObject(sensor); sensorWindow('scale',1);

%% Interpolate the color filter data to produce a full sensor
%
ip = ipCreate;
ip = ipCompute(ip,sensor);
% vcAddObject(ip); ipWindow;

% Show in a separate window
rgb = ipGet(ip,'result');
% vcNewGraphWin; image(lrgb2srgb(rgb)); axis image

%% Pack the samples of the rgb image into the lightfield structure

% If we had a lightfield structure, lf, this could become
%    rgb2lf(rgb,lf)

% These parameters should always be part of an oi lightfield description.
%lightField(i,j, :,:, :) = photons(1:9, 1:9, :); 
sz = oiGet(oi,'size');
superPixelW = sz(2)/in.numPinholesW;
superPixelH = sz(1)/in.numPinholesH;

% This is the array size of pinholes (or microlens)
% The reason for floor() is ... well rounding or something.  Shouldn't
% really be needed.
numSuperPixW = floor(size(rgb, 2)/superPixelW);
numSuperPixH = floor(size(rgb, 1)/superPixelH);

% Allocate space
% lightfield = zeros(numSuperPixW, numSuperPixH, superPixelW, superPixelH, 3);
lightfield = zeros(superPixelH, superPixelW, in.numPinholesW, in.numPinholesH, 3);

% For numerical calculations, we would use this
for i = 1:numSuperPixW
    for j = 1:numSuperPixH
        lightfield(:,:, j, i, :) = ...
            rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end

% For visualization, thismight be a good idea - use lrgb2srgb
% LF = zeros(superPixelH, superPixelW, numSuperPixW, numSuperPixH, 3);
LF = zeros(superPixelH, superPixelW, in.numPinholesW, in.numPinholesH, 3);
rgb = lrgb2srgb(double(rgb));
for i = 1:numSuperPixW
    for j = 1:numSuperPixH
        LF(:,:, j, i, :) = ...
            rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end


%% Some views of the light field data

% lightField(:,:, row, col, :) gives us a view from corresponding pixels in
% each of the pinhole (microlens array) data sets.

% The pixels at the edges don't really get any rays or if they do they get
% very little late (are noisier).

% vcNewGraphWin;
% cnt = 1;
% row = superPixelH; col = superPixelW;
% rList = 1:2:row;
% cList = 1:2:col;
% for rr=rList
%     for cc=cList
%         img = squeeze(LF(rr,cc,:,:,:));
%         subplot(length(rList),length(cList),cnt), imshow(img);
%         cnt = cnt + 1;
%     end
% end

%% create the full tiling of the pinhole images and put it in an image


vcNewGraphWin;
cnt = 1;
row = superPixelH; col = superPixelW;
rList = 1:1:row;
cList = 1:1:col;
tiledImage = zeros(superPixelH * numPinholesH, superPixelW * numPinholesW, 3);

for rr=rList
    for cc=cList
        img = squeeze(LF(rr,cc,:,:,:));
        
        tiledImage((rr-1) * numPinholesH + 1: rr * numPinholesH , (cc-1) * numPinholesW +1: cc * numPinholesW , :) = img;
    end
end


imshow(tiledImage);
%% Autofocus as a function call

slopeRange = [-5 -4];  %change this 
slopeRange = [-1.3 -.9];  %for 160 left
slopeRange = [0 .4]; % for 160 right
slopeRange = [1.4 2.8];
%slopeRange = [-.6 -.3]
stepSize = .05;
bestFocusImage = s3dLFAutofocus(lightfield, [], slopeRange, stepSize, []); 

%In focus Slopes 
%middle: 0
%left: -.35 (80)  -1.15 (160)   -2.68 (240) -4.45 (320)
%right: .15 or . .1(80)     .35(160)   .81 (240)         1.5(320)   

% fig = 

Slope = -1.15
vcNewGraphWin;
ShiftImg = LFFiltShiftSum(lightfield, Slope );
imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
bestFocusImage = ShiftImg(:,:,1:3);
axis image; truesize
title(sprintf('Parameter %0.2f',Slope))


%% perform slanted bar analysis (automatically) - eventually put this in another script

% fig = vcNewGraphWin
% 
% Slope = 0;
% ShiftImg = LFFiltShiftSum(lightfield, Slope );
% imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
% axis image; truesize
% title(sprintf('Parameter %0.2f',Slope))

%masterRect = round(getrect(fig));

% Run the computation for the monochrome sensor
% sensor = sensorCompute(sensorC,oi);
% vcReplaceObject(sensor);   
% vci = ipCompute(vci,sensor);
% vcReplaceObject(vci); ipWindow;
% 
% % Find a good rectangle
% masterRect = ISOFindSlantedBar(vci);
% h = ieDrawShape(vci,'rectangle',masterRect);



% these are for 80 x 80 x 18 x 18 LF's.  Multiply for different spatial
% resolutions...

% We want to try... 80 x 80 x 36 x 36
%                   120 x 120 x 24 x 24 
%                   160 x 160 x 18 x 18  
%                   240 x 240 x 12 x 12
%                   320 x 320 x 9 x 9
%

scaleFactor = numPinholesW/80;   %depending on number of pinholes, 
middleRect = [37 33 7 15] * scaleFactor;   
rightRect = [60 35 5 9] * scaleFactor;
leftRect = [15 31 8 19] * scaleFactor;
masterRect = leftRect;

ip = ipSet(ip, 'result', bestFocusImage);

barImage = vcGetROIData(ip,masterRect,'result');
c = masterRect(3)+1;
r = masterRect(4)+1;
barImage = reshape(barImage,r,c,3);
% vcNewGraphWin; imagesc(barImage(:,:,1)); axis image; colormap(gray);

% Run the ISO 12233 code.  The results are stored in the window.


%ss = oiGet(oi,'sample spacing','m');
originalSensorRes = oiGet(oi, 'size');
ss = oiGet(oi,'sample spacing','m');
ss = ss .* originalSensorRes(1)/numPinholesW;  %sample size accounts for integration of light field
sensor = sensorCreate;
sensor = sensorSet(sensor,'pixel size same fill factor',ss(1));
sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
sensor = sensorSet(sensor,'exp time',0.010);

% pixel = sensorGet(sensor, 'pixel')
sensor = sensorSet(sensor, 'dsnu level', .004);

% Describe
sensorGet(sensor,'pixel size','um')
sensorGet(sensor,'size')
sensorGet(sensor,'fov',[],oi)


sensor = sensorSet(sensor, 'rows', numPinholesH);
sensor = sensorSet(sensor, 'cols', numPinholesW);
pixel = sensorGet(sensor,'pixel');
dx = pixelGet(pixel,'width','mm');
[results, fitme, esf, h] = ISO12233(barImage, dx, [] , 'luminance') 
results.mtf50



%%  Interact with the lightfield using the toolbox

% In this case, we use the srgb representation because we are just
% visualizing
% LFDispMousePan(LF)

%
LFDispVidCirc(LF)

%% Results
% Spatial Samples	Left	Middle	Right
% 80	1.2	1.2	1
% 160	1.2	2.6	2.2
% 240	1	3.2	2.4
% 320	0.8	3.6	2
spatialSamples = [80 160 240 320];
leftMTF = [1.2 1.2 1 .8];
middleMTF = [1.2 2.6 3.2 3.6];
rightMTF = [1 2.2 2.4 2];

vcNewGraphWin;
plot(spatialSamples, leftMTF, spatialSamples, middleMTF, spatialSamples, rightMTF);
legend('Left MTF50', 'Middle MTF50', 'Right MTF50');
xlabel('Number of Spatial Samples');
ylabel('MTF50');
title('Spatial vs. Angular Resolution Trade-offs');


%% compare diffraction vs. no diffraction

diffractRes = load('MTFMiddle320.mat');
noDiffractRes = load('MTFMiddle320NoDiffract.mat');

freq = diffractRes.results.freq;
vcNewGraphWin;
plot(freq, diffractRes.results.mtf(:,1), 'r', freq, noDiffractRes.results.mtf(:,1), 'b');

legend('Diffraction', 'No Diffraction');

vcNewGraphWin;
plot(freq, noDiffractRes.results.mtf(:,1));



%% compare mtf plots for different spatial/angular tradeoffs
mtfLeft320 = load('MTFLeft320.mat');
%mtfLeft240 = load('MTFLeft240.mat');
mtfLeft160 = load('MTFLeft160.mat');
mtfLeft80 = load('MTFLeft80.mat');


vcNewGraphWin;
plot(mtfLeft320.results.freq, mtfLeft320.results.mtf(:,1), mtfLeft160.results.freq, mtfLeft160.results.mtf(:,1) ,mtfLeft80.results.freq,  mtfLeft80.results.mtf(:,1));
legend('320 Pixel Cells', '160 Pixel Cells','80 Pixel Cells');
xlabel('Spatial frequency (cy/mm on sensor)');
ylabel('Contrast reduction (SFR)');
xlim([0 7])

