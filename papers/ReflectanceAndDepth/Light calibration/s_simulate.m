close all;
clear variables;
clc;

load('angledPlane.mat');
vcAddObject(scene); sceneWindow();

oi = oiCreate();
oi = oiCompute(scene,oi);
vcAddObject(oi); oiWindow();

sensor = sensorCreate('monochrome');
sensor = sensorSet(sensor,'quantizationmethod','8bit');
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'fov'));
sensor = sensorCompute(sensor,oi);

vcAddObject(sensor); sensorWindow();

img = uint8(sensorGet(sensor,'dv'));

[X, Y] = meshgrid(1:size(img,2),1:size(img,1));


aEstVec = [];
bEstVec = [];
xcEstVec = [];
ycEstVec = [];
tauEstVec = [];

for i=255:-1:220

    indx = img == i;
    nPoints = sum(indx(:));
    if nPoints < 10, continue; end;
    
    x = X(indx);
    y = Y(indx);
    
    [ Cmat, aEst, bEst, xcEst, ycEst, tauEst ] = fitEllipse( x(:),y(:) );
    
    aEstVec = [aEstVec, aEst];
    bEstVec = [bEstVec, bEst];
    xcEstVec = [xcEstVec, xcEst];
    ycEstVec = [ycEstVec, ycEst];
    tauEstVec = [tauEstVec, tauEst];
    
    rotMatEst = [cosd(tauEst) -sind(tauEst); sind(tauEst) cosd(tauEst)];
    cVecEst = [xcEst; ycEst];

    angles = 1:360; nAngles = length(angles);
    resEll = repmat(cVecEst,1,nAngles) + rotMatEst*[aEst*cosd(angles); bEst*sind(angles)];
    
    resCirc = [cosd(angles); sind(angles)];
 
    figure;
    hold on; grid on; box on;
    imagesc(img == i);
    plot(resEll(1,:),resEll(2,:),'go');

end

%% Fit params
ellipseCoords = (resEll - repmat(cVecEst,1,nAngles));
circleCoords = [resCirc; ones(1,nAngles)];
cvx_begin
    variables R(2,3)
    minimize norm(ellipseCoords - R*circleCoords,'fro');
    subject to
        R(2,1) == 0
        cvx_end



tx = zeros(length(aEstVec),1);
ty = zeros(length(aEstVec),1);
for t=1:length(aEstVec)
    Tr = [1 0 -xcEstVec(t); 0 1 -ycEstVec(t); 0 0 1];
    rotMatEst = [cosd(tauEstVec(t)) -sind(tauEstVec(t)); sind(tauEstVec(t)) cosd(tauEstVec(t))];
    Rot = [inv(rotMatEst) zeros(2,1); zeros(1,2) 1];
    Sc = [1/aEstVec(t) 0 0; 0 1/bEstVec(t) 0; 0 0 1];
    M = Sc*Rot*Tr;
    
    tx(t) = M(1,3);
    ty(t) = M(2,3);
    
end

plot(xcEstVec,ycEstVec,'rs');



