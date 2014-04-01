function scene = s3dRenderScene(fullfname, sceneName)
% Use a pinhole camera model to calculate scene radiance
%
% fullfname: full path (file name) of the pbrt file that we render
% path:  
%
% This function renders a PBRT scene using the given a full fname, then
% returns the data as an optical image. The pbrt file should be previously 
% generated and be in datapath/generatedPbrtFiles/.  The output file is placed
% in datapath/generatedPbrtFiles/tempPBRT/.  The
% proper optics is also placed into this optical image that corresponds to
% the focal length and field of view (FOV).

%
% Todo: We are considering copying the pbrt file that was used to generate
% the scene when we save the scene for future use.
%
% Todo: figure out consistency between s3dRenderScen and s3dRenderOi
%%  


if (ieNotDefined('fullfname')) || ~exist(fullfname,'file')
    error('PBRT full file name required.  File not found');
end

% Use pinhole and pbrt to create the scene data
pbrtExe = fullfile(pbrtRootPath, 'src','bin','pbrt');
if ~exist(pbrtExe,'file')
    error('PBRT executable not found');
end

% Make a tempPBRT directory where the output files will go
generatedDir = fullfile(dataPath, 'generatedPbrtFiles', 'tempPBRT');
if exist(generatedDir,'dir')
    unix(['rm ' fullfile(generatedDir, '*')]);
else
    mkdir(generatedDir);
end
outfile  = fullfile(generatedDir, 'temp_out.dat');

% [p,n,ext] = fileparts(fullfname);
cmd = sprintf('%s %s --outfile %s',pbrtExe,fullfname,outfile);

% chdir(p)
unix(cmd)

%% ISET will read the PBRT output and convert to a scene

scene = pbrt2scene(outfile);
%rename the oi, if a name is given
if (~ieNotDefined('sceneName'))
    scene = sceneSet(scene, 'name', sceneName);
    scene = sceneAdjustLuminance(scene,100);
end

end