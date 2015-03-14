%  Experiments with light fields and ISET

%
% (AL) Vistasoft Team, 2015

%% inputs

%load LF as oi
%in = load(fullfile(dataPath, 'lightfields', 'benchLFSceneDirect.mat'))
in = load(fullfile(dataPath, 'lightfields', 'indObjLFOiDirect.mat'));
oi = in.opticalimage;

%convert to RGB (we are skipping the sensor processing step for simplicity.
%That will come later)
rgb = oiGet(oi, 'rgb');


%% process sensor and image processing data



%% Convert to a 5D matrix form

%lightField(i,j, :,:, :) = photons(1:9, 1:9, :); 
superPixelW = 9;
superPixelH = 9;

numSuperPixW = size(rgb, 2)/superPixelW;
numSuperPixH = size(rgb, 1)/superPixelH;

lightfield = zeros(numSuperPixW, numSuperPixH, superPixelW, superPixelH, 3);

for i = 1:numSuperPixW
    for j = 1:numSuperPixH
        lightfield(i,j, :, :, :) = rgb(((j-1)*superPixelH + 1):(j*superPixelH), ...
            ((i-1) * superPixelW + 1):(i*superPixelW), :);
    end
end

%% Some views
% lightField(:,:, row, col, :) gives us a pinhole view at a specific
% position on the aperture. 
%  Notice that the pixels at the edges don't really get any rays or if they
%  do they get very little late (are noisier).
vcNewGraphWin;
cnt = 1;
row = superPixelH; col = superPixelW;
for rr=1:row
    for cc=1:col
        img = squeeze(lightfield(:,:,rr,cc,:));
        img = imageTranspose(img);
        
        subplot(9,9,cnt), imshow(img);
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


