%% s_pbrtExample
%
%
%
% Illustrate how to write a batch script to render a small scene multiple
% times with pbrt, using pbrtCreate, pbrtSet, pbrtWrite and so forth
%
%

%% simple case
pbrt=pbrtCreate
pbrtWrite(pbrt, 'deleteMe.pbrt');
oi = s3dRenderScene('deleteMe.pbrt', 50, [dataPath '/tmp/'])


%% example batch job using a different light source for each job

% make a cell array of pbrt objects with different lighting and execute them in batch
pbrtArray = cell(1,1)

%General spotlight
light.type = 'spot';
light.spectrum.type = 'rgb I';
light.spectrum.value = [1000 1000 1000];
light.coneangle = 180;
light.conedeltaangle = 180;
light.from = [4.5 -90 8.5];
light.to = [4.5 -89 8.5];

%green spotlight
greenLight = light;
greenLight.spectrum.value = [0 1000 0];
pbrtArray{1} = pbrt;
pbrtArray{1}.lightSource{1} = greenLight;
pbrtArray{1}.name = 'greenLight';

%red spotlight
redLight = light;
redLight.spectrum.value = [1000 0 0];
pbrtArray{2} = pbrt;
pbrtArray{2}.lightSource{1} = redLight;
pbrtArray{2}.name = 'redLight';

%blue spotlight
blueLight = light;
blueLight.spectrum.value = [0 0 1000];
pbrtArray{3} = pbrt;
pbrtArray{3}.lightSource{1} = blueLight;
pbrtArray{3}.name = 'blueLight';

%white spotlight
whiteLight= light; 
pbrtArray{4} = pbrt;
pbrtArray{4}.lightSource{1} = whiteLight;
pbrtArray{4}.name = 'whiteLight';

%% execute this array in batch
for i = 1:length(pbrtArray)
    fileName = [pbrtArray{i}.name '.pbrt'];
    pbrtWrite(pbrtArray{i}, fileName);   
    oi = s3dRenderScene(fileName, 50, [dataPath '/tmp/'], fileName);
end


%% example of using the new pbrtObject class
clear testPbrtObject;
testPbrtObject = pbrtObject();

% declare an additional spotlight
clear blueLight;
% light.type = 'light';
% light.lightType = 'spot';
% light.spectrum.type = 'rgb I';
% light.spectrum.value = [1000 1000 1000];
% light.coneangle = 180;
% light.conedeltaangle = 180;
% light.from = [4.5 -90 8.5];
% light.to = [4.5 -89 8.5];


blueLight = lightObject();
blueLight.spectrum.setValue([0 0 1000]);
%add this light source to the object
testPbrtObject.addLightSource(blueLight);

testPbrtObject.writeFile('deleteMe.pbrt');
oi = s3dRenderScene('deleteMe.pbrt', 50, [dataPath '/tmp/'], 'deleteMe.pbrt');

%% End
