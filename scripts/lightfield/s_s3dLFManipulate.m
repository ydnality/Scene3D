%  Experiments with light fields and ISET
%
%  See also: s_s3dISETLF.m, LFToolbox0.4
%
%  To retrieve some of the lightfield data use
%     urlBase = 'http://scarlet.stanford.edu/validation/SCIEN/Lightfield';
%     fname = 'indObjLFOiDirect.mat';
%     fname = 'benchLFScene.mat';
%     urlwrite(fullfile(urlBase,fname),fname)
% 
% (AL) Vistasoft Team, 2015

%%
ieInit

%% inputs

%load a lightfield as an oi object
%in = load(fullfile(dataPath, 'lightfields', 'benchLFSceneDirect.mat'))

% How did this get created?  Through PBRT.
% The way in which got created should be up at scarlet
% This is a script named XXX
in = load(fullfile(dataPath, 'lightfields', 'indObjLFOiDirect.mat'));
oi = in.opticalimage;

% The optics should be reasonable.  We set a focal length of 3.5mm here.
oi = oiSet(oi,'optics focal length',3.5e-3);

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

% These parameters should always be part of an oi lightfield description.
%lightField(i,j, :,:, :) = photons(1:9, 1:9, :); 
superPixelW = 9;
superPixelH = 9;

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
sensor = sensorSet(sensor,'size',oiGet(oi,'size')+[superPixelH superPixelW]);
sensor = sensorSet(sensor,'exp time',0.010);

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
vcAddObject(ip); ipWindow;

% Show in a separate window
rgb = ipGet(ip,'result');
vcNewGraphWin; image(lrgb2srgb(rgb));

%% Convert to a 5D matrix form

numSuperPixW = floor(size(rgb, 2)/superPixelW);
numSuperPixH = floor(size(rgb, 1)/superPixelH);

lightfield = zeros(numSuperPixW, numSuperPixH, superPixelW, superPixelH, 3);

for i = 1:numSuperPixW
    for j = 1:numSuperPixH
        lightfield(i,j, :, :, :) = rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end

%% Some views
% Use the light field toolbox to calculate stuff?

% lightField(:,:, row, col, :) gives us a pinhole view at a specific
% position on the aperture. 
%  Notice that the pixels at the edges don't really get any rays or if they
%  do they get very little late (are noisier).
vcNewGraphWin;
cnt = 1;
row = superPixelH; col = superPixelW;
rList = 1:3:row;
cList = 1:3:col;
for rr=rList
    for cc=cList
        img = squeeze(lightfield(:,:,rr,cc,:));
        img = imageTranspose(img);
        img = lrgb2srgb(img);
        subplot(length(rList),length(cList),cnt), imshow(img);
        cnt = cnt + 1;
    end
end


%% render some example images

%sum all the sub aperture views (3rd and 4th dimensions)  
summedimage = sum(sum(lightfield, 3), 4);
summedimage = reshape(summedimage, [80 80 3]);
summedimage = permute(summedimage, [2 1 3]);
%summedimage = summedimage(:,end:-1:1, :);
summedimage = summedimage./(superPixelW * superPixelH); %normalize image by the number of summed images

vcNewGraphWin; imshow(summedimage);


%% LF workshop stuff

%change into LF workshop format
LF = permute(lightfield, [4 3 1 2 5]);
ShiftSumSlope1 = .2;
ShiftSumSlope2 = 0;

%---Demonstrate shift sum filter---
for( Slope = -.3:.1:.3 )
	fprintf('Applying shift sum filter');
	[ShiftImg, FiltOptionsOut] = LFFiltShiftSum( LF, Slope );
	fprintf(' Done\n');
	FiltOptionsOut
	
	%LFFigure(CurFigure); 
	%CurFigure = CurFigure + 1;
    figure;
    ShiftImg = permute(ShiftImg, [2 1 3]);
    imshow(ShiftImg(:,:,1:3));
	%LFDisp(ShiftImg);
	axis image off
 	truesize
	title(sprintf('Shift sum filter, slope %.3g', Slope));
	drawnow
end



%% Experiment with "hyperfan" filter.  I initially thought this could render all in focus images, but I appera to be incorrect...

%---Demonstrate 4D Hyperfan filter---
LFSize = size(LF);
HyperfanSlope1 = 0; HyperfanSlope2 = .3;
HyperfanBW = 0.035;  % What does this mean?!?!

fprintf('Building 4D frequency hyperfan... ');
[H, FiltOptionsOut] = LFBuild4DFreqHyperfan( LFSize, HyperfanSlope1, HyperfanSlope2, HyperfanBW );
fprintf('Applying filter');
[LFFilt, FiltOptionsOut] = LFFilt4DFFT( LF, H, FiltOptionsOut );
FiltOptionsOut

% LFFigure(CurFigure);
% CurFigure = CurFigure + 1;
%figure; 
%LFFilt = permute(LFFilt, [2 1 3]);
%imshow(LFFilt);
figure;
LFDisp(LFFilt);

%sum all the sub aperture views (3rd and 4th dimensions)  
% summedimage = sum(sum(LFFilt(:,:,:,:, 1:3), 1), 2);
% summedimage = reshape(summedimage, [80 80 3]);
% summedimage = permute(summedimage, [2 1 3]);
% %summedimage = summedimage(:,end:-1:1, :);
% summedimage = summedimage./(superPixelW * superPixelH); %normalize image by the number of summed images
% vcNewGraphWin; imshow(summedimage);


axis image off
truesize
title(sprintf('Frequency hyperfan filter, slopes %.3g, %.3g, HyperfanBW %.3g', HyperfanSlope1, HyperfanSlope2, HyperfanBW));
drawnow


