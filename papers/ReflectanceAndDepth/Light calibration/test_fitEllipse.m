close all;
clear variables;
clc;

xc = 10;
yc = -10;
a = 5;
b = 2;
tau = 33;

angles = 1:360;
nAngles = length(angles);
rotMat = [cosd(tau) -sind(tau); sind(tau) cosd(tau)];
cVec = [xc; yc];

res = repmat(cVec,1,nAngles) + rotMat*[a*cosd(angles); b*sind(angles)] + 0.1*randn(2,nAngles);

figure;
hold on; grid on; box on;
plot(res(1,:),res(2,:),'x');
axis square; axis equal;

% Fit an ellipse
est = fitEllipse(res(1,:)',res(2,:)');

A = est(1);
B = est(2);
C = est(3);
D = est(4);
E = est(5);
F = est(6);

M0 = [F D/2 E/2; D/2 A B/2; E/2 B/2 C];
M = [A B/2; B/2 C];
eVal = sort(eig(M));

aEst = sqrt(-det(M0)/(det(M)*eVal(1)));
bEst = sqrt(-det(M0)/(det(M)*eVal(2)));
xcEst = (B*E - 2*C*D)/(4*A*C - B^2);
ycEst = (B*D - 2*A*E)/(4*A*C - B^2);
tauEst = acotd((A-C)/B)/2;

rotMatEst = [cosd(tauEst) -sind(tauEst); sind(tauEst) cosd(tauEst)];
cVecEst = [xcEst; ycEst];

angles = 1:360; nAngles = length(angles);
resEll = repmat(cVecEst,1,nAngles) + rotMatEst*[aEst*cosd(angles); bEst*sind(angles)];

plot(resEll(1,:),resEll(2,:),'g');



%% 'Reshape' the ellipse into a circle

resCirc = diag([1/aEst, 1/bEst])*inv(rotMatEst)*(resEll - repmat([xcEst; ycEst],1,size(resEll,2)));

xNew = resCirc(1,:);
yNew = resCirc(2,:);

figure;
hold on; grid on; box on;
plot(xNew,yNew,'og');
axis square;

resEllHomog = [resEll; ones(1,size(resEll,2))];

Tr = [1 0 -xcEst; 0 1 -ycEst; 0 0 1];
Rot = [inv(rotMatEst) zeros(2,1); zeros(1,2) 1];
Sc = [1/aEst 0 0; 0 1/bEst 0; 0 0 1];
M = Sc*Rot*Tr;

resCircHomog = M*resEllHomog;

figure; 
hold on; grid on; box on;
plot(resCircHomog(1,:), resCircHomog(2,:),'c^');
axis equal;

