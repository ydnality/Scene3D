% This script renders the basic setup of 2 side-by-side light sources and a
% simple surface in between them.  There will be 2 captures.  One with each
% light source.  The hope is that we can derive some information on the
% normal vectors using this setup.  These captures are rendered using
% pbrtObjects. 




s_initISET
%%
clear curPbrt;
curPbrt = pbrtObject();

%add a new material
matRGB= [1 1 1];
newMaterial = pbrtMaterialObject('grayMat', 'matte', pbrtPropertyObject('color Kd', matRGB));
curPbrt.addMaterial(newMaterial);
% newCamPos = [.0001 .0001 0;   %this matrix must be invertible? why?
%              .0001 .0001 -1;
%              .000001 .0000001 .9999999];
% curPbrt.camera.setPosition(newCamPos);


% newCamPos =    [0  -80.0000    0;
%     0  -79.0000    0;
%          0         0    1.0000];
     
newCamPos =    [0  0 80.0000;
    0   0 79.0000;
    0 1.00000 0];
curPbrt.camera.setPosition(newCamPos);

%add new geoemtry
translateTransform = [1 0 0 0;
    0 1 0 0 ;
    0 0 1 0;
    0 0 0  1]; %8.87306690216     %x direction is to the right, y is into the screen, z is up


theta = 45 * pi/180;
rotationTransform = [cos(theta) 0 sin(theta) 0;  %rotation about the y axis
                     0 1 0 0
                     -sin(theta) 0 cos(theta) 0;
                     0 0 0 1];
totalTransform = rotationTransform * translateTransform;


lightLeft = pbrtLightSpotObject('leftLight', [], [], [], [-50 0 80], [-50 0 79]);
lightRight = pbrtLightSpotObject('rightLight', [], [], [], [50 0 80], [50 0 79]);
% lightRight = pbrtLightSpotObject('rightLight', [], [], [], inFrom, inTo);
curPbrt.removeLight();
% curPbrt.addLightSource(lightLeft);
curPbrt.addLightSource(lightRight);

newGeometry = pbrtGeometryObject('newGeom', 'grayMat', [], [], totalTransform);

% newGeometry.setTransform(eye(4));
curPbrt.removeGeometry();
curPbrt.addGeometry(newGeometry);

tmpFileName = ['deleteMe' '.pbrt'];
curPbrt.writeFile(tmpFileName);
oi = s3dRenderOI(curPbrt, .050, tmpFileName);
% end
