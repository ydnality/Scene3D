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


%% declare some light Sources

%green spotlight
greenLight.type = 'spot';
greenLight.spectrum.type = 'rgb I';
greenLight.spectrum.value = [0 1000 0];
greenLight.coneangle = 180;
greenLight.conedeltaangle = 180;
greenLight.from = [4.5 -90 8.5];
greenLight.to = [4.5 -89 8.5];

%red spotlight
redLight.type = 'spot';
redLight.spectrum.type = 'rgb I';
redLight.spectrum.value = [1000 0 0];
redLight.coneangle = 180;
redLight.conedeltaangle = 180;
redLight.from = [4.5 -90 8.5];
redLight.to = [4.5 -89 8.5];

%blue spotlight
blueLight.type = 'spot';
blueLight.spectrum.type = 'rgb I';
blueLight.spectrum.value = [0 0 1000];
blueLight.coneangle = 180;
blueLight.conedeltaangle = 180;
blueLight.from = [4.5 -90 8.5];
blueLight.to = [4.5 -89 8.5];

%white spotlight
whiteLight.type = 'spot';
whiteLight.spectrum.type = 'rgb I';
whiteLight.spectrum.value = [1000 1000 1000];
whiteLight.coneangle = 180;
whiteLight.conedeltaangle = 180;
whiteLight.from = [4.5 -90 8.5];
whiteLight.to = [4.5 -89 8.5];


%% make a cell array of pbrt objects with different lighting and execute them in batch
pbrtArray = cell(1,1)
pbrtArray{1} = pbrt;
pbrtArray{1}.lightSource{1} = greenLight;
pbrtArray{2} = pbrt;
pbrtArray{2}.lightSource{1} = redLight;
pbrtArray{3} = pbrt;
pbrtArray{3}.lightSource{1} = blueLight;
pbrtArray{4} = pbrt;
pbrtArray{4}.lightSource{1} = whiteLight;

%execute this array in batch

for i = 1:length(pbrtArray)
    fileName = ['deleteMe' int2str(i) '.pbrt'];
    pbrtWrite(pbrtArray{i}, fileName);   
    oi = s3dRenderScene(fileName, 50, [dataPath '/tmp/']);
end

%% End
