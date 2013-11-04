%% new depth-map method for perspective camera - consider putting this into a function

[path,name,ext] = fileparts(outfile);
dMapFile = [path '/' name '_DM.dat'];
badDepthMap = s3dReadDepthMapFile(dMapFile);
newDepthMap = badDepthMap;
numSamples = 64;

for j = 1:size(badDepthMap,1);
    for i = j+1:size(badDepthMap,2);
        newDepthMap(j, i) =  badDepthMap(j, i)/numSamples;
    end
end

figure; imagesc(badDepthMap);
figure; imagesc(newDepthMap);