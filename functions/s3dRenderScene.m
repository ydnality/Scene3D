function scene = s3dRenderScene(fullfname, sceneName)
% Use a pinhole camera model to calculate scene radiance
%
% fullfname: full path (file name) of the pbrt file that we render
% path:  
%
% This function renders a PBRT scene using the given pbrt fname, then
% returns the data as an optical image. A temporary directory is created to
% render the PBRT scene.  The contents of this directory are deleted at the
% function call to the particular scene.  The oi is added to an oiWindow,
% but the oiWindow is not displayed by default.  The directory that this
% function looks for the pbrt file is s3dRooth/scripts/pbrtFiles/.  The
% proper optics is also placed into this optical image that corresponds to
% the focal length and field of view (FOV).
%
% PBRT creates some temporary files.  We put these in a local ./tempPBRT
% and then we delete them at the end.
%
% Todo: We are considering copying the pbrt file that was used to generate
% the scene when we save the oi for future use.

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
    unix(['rm ' generatedDir]);
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