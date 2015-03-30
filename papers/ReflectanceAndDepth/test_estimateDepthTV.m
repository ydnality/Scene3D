% close all;
% clear variables;
% clc;

%load('/Users/hblasinski/Dropbox/reflectanceanddepth/Data/03192015_depthEstimation/a_b_c_groundTruthDepth.mat');
load([s3dRootPath '/papers/ReflectanceAndDepth/Data/03192015_depthEstimation/a_b_c_groundTruthDepth.mat']);

c(isnan(c)) = 0;

delta = 100;  %used to be 10

[est, init, hist] = estimateDepthTV(a,b,c,delta,'maxIter',150,'verbose',true);

figure; 
subplot(1,3,1);
imagesc(real(init),[0 130]);
title('Initial estimate');

subplot(1,3,2);
imagesc(real(est),[0 130]);
title(sprintf('TV %i',delta));
subplot(1,3,3);
imagesc(groundTruthDepthMap,[0 130]);
title('Ground truth');

figure;
semilogy([hist.prRes, hist.dualRes]);
legend('Primal','Dual');
title('Residuals');

figure;
plot(hist.rho);
title('Rho');

