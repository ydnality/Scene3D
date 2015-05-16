    %  Experiments with light fields and ISET
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
in = load('slantedBarMultLF320Diffract.mat');
% in = load('metronomeLF.mat');

oi = in.oi;

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
sensor = sensorSet(sensor,'exp time',0.010);

% pixel = sensorGet(sensor, 'pixel')
sensor = sensorSet(sensor, 'dsnu level', .004);

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

vcNewGraphWin;
cnt = 1;
row = superPixelH; col = superPixelW;
rList = 1:2:row;
cList = 1:2:col;
for rr=rList
    for cc=cList
        img = squeeze(LF(rr,cc,:,:,:));
        subplot(length(rList),length(cList),cnt), imshow(img);
        cnt = cnt + 1;
    end
end

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

%% Compare the leftmost and rightmost images in the middle
% vcNewGraphWin([],'wide');
% 
% img = squeeze(LF(3,2,:,:,:));
% subplot(1,2,1), imagescRGB(img);
% 
% img = squeeze(LF(3,8,:,:,:));
% subplot(1,2,2), imagescRGB(img);

%% render some example images

% If we sume all the r,g and b pixels within each aperture we get a single
% RGB image corresponding to the mean.  This is the image at the microlens
% itself
% tmp = squeeze(sum(sum(LF,2),1));
% vcNewGraphWin;
% tmp = tmp/max(tmp(:));
% imagescRGB(tmp);

%% Now what would the image have been if we move the sensor forward?

% -0.8:.2:0.6 for indestructible object.
% -.6 to 3.4 for bench light field
fig = vcNewGraphWin
for Slope = -0.5:.25:.5  %1.5 
%for Slope = [0 .25 1.75]   
    ShiftImg = LFFiltShiftSum(lightfield, Slope );
    imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
    axis image; truesize
    title(sprintf('Parameter %0.2f',Slope))
    pause(0.2)
    tmpImg = ShiftImg(:,:,1:3);
  %  imwrite(tmpImg./max(tmpImg(:)), sprintf('BenchLFImageDSNUSlope%.2f.png', Slope))

end


%.175 is in focus for far slanted bar



%% Autofocus as a function call

slopeRange = [-5 -4]  %change this 
stepSize = .05
bestFocusImage = s3dLFAutofocus(lightfield, [], slopeRange, stepSize, []); 

%In focus Slopes 
%middle: 0
%left: -.6 (80)      -4.45 (320)
%right: .15

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
middleRect = [37 33 7 15];   
rightRect = [60 35 5 9];
leftRect = [15 31 8 19] * 4;
masterRect = leftRect;

ip = ipSet(ip, 'result', ShiftImg(:,:,1:3));

barImage = vcGetROIData(ip,masterRect,'result');
c = masterRect(3)+1;
r = masterRect(4)+1;
barImage = reshape(barImage,r,c,3);
% vcNewGraphWin; imagesc(barImage(:,:,1)); axis image; colormap(gray);

% Run the ISO 12233 code.  The results are stored in the window.
pixel = sensorGet(sensor,'pixel');
dx = pixelGet(pixel,'width','mm');
[results, fitme, esf, h] = ISO12233(barImage, dx, [] , 'luminance') 
results.mtf50



%% Autofocus a section of the image (original experimental script)

% rect = getrect(fig)
% roundRect = round(rect);
% %select a rectangle on the last figure
% index = 1;
% SlopeVec = -0.7:.05:.5;
% varianceVec = zeros(size(SlopeVec));
% 
% for Slope = SlopeVec  %1.5 
%     ShiftImg = LFFiltShiftSum(lightfield, Slope );
%     imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
%     tmpImg = ShiftImg(:,:,1:3);
%     
%     varianceVec(index) = sum(sum(var(tmpImg(roundRect(2):roundRect(2) + roundRect(4), roundRect(1): roundRect(1) + roundRect(3), :)))); 
%     index = index + 1;
% end
% 
% [maxVar, maxInd] = max(varianceVec)
% 
% Slope = SlopeVec(maxInd);
% vcNewGraphWin;
% ShiftImg = LFFiltShiftSum(lightfield, Slope );
% imagescRGB(lrgb2srgb(ShiftImg(:,:,1:3)));
% axis image; truesize
% title(sprintf('Parameter %0.2f',Slope))
% tmpImg = ShiftImg(:,:,1:3);


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

%
LFDispVidCirc(LF)

%%
% Change into LF workshop format.
% We should do this up above, if this is true.
% ShiftSumSlope1 = .2;
% ShiftSumSlope2 = 0;

%---Demonstrate shift sum filter---
% for( Slope = -.3:.1:.3 )
% 	fprintf('Applying shift sum filter');
% 	[ShiftImg, FiltOptionsOut] = LFFiltShiftSum( LF, Slope );
% 	fprintf(' Done\n');
% 	FiltOptionsOut
% 	
% 	%LFFigure(CurFigure); 
% 	%CurFigure = CurFigure + 1;
%     figure;
%     ShiftImg = permute(ShiftImg, [2 1 3]);
%     imshow(ShiftImg(:,:,1:3));
% 	%LFDisp(ShiftImg);
% 	axis image off
%  	truesize
% 	title(sprintf('Shift sum filter, slope %.3g', Slope));
% 	drawnow
% end



%% Experiment with "hyperfan" filter.  I initially thought this could render all in focus images, but I appera to be incorrect...

%---Demonstrate 4D Hyperfan filter---
% LFSize = size(LF);
% HyperfanSlope1 = 0; HyperfanSlope2 = .3;
% HyperfanBW = 0.035;  % What does this mean?!?!
% 
% fprintf('Building 4D frequency hyperfan... ');
% [H, FiltOptionsOut] = LFBuild4DFreqHyperfan( LFSize, HyperfanSlope1, HyperfanSlope2, HyperfanBW );
% fprintf('Applying filter');
% [LFFilt, FiltOptionsOut] = LFFilt4DFFT( LF, H, FiltOptionsOut );
% FiltOptionsOut
% 
% % LFFigure(CurFigure);
% % CurFigure = CurFigure + 1;
% %figure; 
% %LFFilt = permute(LFFilt, [2 1 3]);
% %imshow(LFFilt);
% figure;
% LFDisp(LFFilt);

%sum all the sub aperture views (3rd and 4th dimensions)  
% summedimage = sum(sum(LFFilt(:,:,:,:, 1:3), 1), 2);
% summedimage = reshape(summedimage, [80 80 3]);
% summedimage = permute(summedimage, [2 1 3]);
% %summedimage = summedimage(:,end:-1:1, :);
% summedimage = summedimage./(superPixelW * superPixelH); %normalize image by the number of summed images
% vcNewGraphWin; imshow(summedimage);


% axis image off
% truesize
% title(sprintf('Frequency hyperfan filter, slopes %.3g, %.3g, HyperfanBW %.3g', HyperfanSlope1, HyperfanSlope2, HyperfanBW));
% drawnow


