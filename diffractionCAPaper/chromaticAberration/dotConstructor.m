testTarget = zeros(512,512,3);

% testTarget(206, 256, 1) = 1;
% testTarget(206, 256, 3) = 1;
% testTarget(306, 256, 1) = 1;
% testTarget(306, 256, 3) = 1;

testTarget(128, 256, 1) = 1;
testTarget(128, 256, 3) = 1;
testTarget(384, 256, 1) = 1;
testTarget(384, 256, 3) = 1;
figure; imshow(testTarget);
imwrite(testTarget, '2PinkDots.tif');