%% still under construction 

clear pbrt;
templatePbrt = pbrtObject();

% spectrum list
lightList = [1000 0 0;
    0 1000 0;
    0 0 1000;
    1000 1000 1000];


% declare an additional spotlight
light.type = 'light';
light.lightType = 'spot';
light.spectrum.type = 'rgb I';
light.spectrum.value = [1000 1000 1000];
light.coneangle = 180;
light.conedeltaangle = 180;
light.from = [4.5 -90 8.5];
light.to = [4.5 -89 8.5];

%loop through all spectrums and render an image
for i = 1:size(lightList, 1)
   curPbrt = templatePbrt;
   curPbrt.light.spectrum.value = lightList(i, :);
   curPbrt.writeFile('deleteMe.pbrt');
   oi = s3dRenderScene('deleteMe.pbrt', 50, [dataPath '/tmp/'], 'deleteMe.pbrt');
end



% 
% 
% clear blueLight;
% blueLight = light;
% blueLight.spectrum.value = [0 0 1000];
% 
% %green spotlight
% greenLight = light;
% greenLight.spectrum.value = [0 1000 0];
% pbrtArray{1} = pbrt;
% pbrtArray{1}.lightSource{1} = greenLight;
% pbrtArray{1}.name = 'greenLight';
% 
% %red spotlight
% redLight = light;
% redLight.spectrum.value = [1000 0 0];
% pbrtArray{2} = pbrt;
% pbrtArray{2}.lightSource{1} = redLight;
% pbrtArray{2}.name = 'redLight';
% 
% %blue spotlight
% blueLight = light;
% blueLight.spectrum.value = [0 0 1000];
% pbrtArray{3} = pbrt;
% pbrtArray{3}.lightSource{1} = blueLight;
% pbrtArray{3}.name = 'blueLight';
% 
% %white spotlight
% whiteLight= light; 
% pbrtArray{4} = pbrt;
% pbrtArray{4}.lightSource{1} = whiteLight;
% pbrtArray{4}.name = 'whiteLight';
% 
% 
% 
% 
% %add this light source to the object
% testPbrtObject.addLightSource(blueLight);
% 
% testPbrtObject.writeFile('deleteMe.pbrt');
% oi = s3dRenderScene('deleteMe.pbrt', 50, [dataPath '/tmp/'], 'deleteMe.pbrt');