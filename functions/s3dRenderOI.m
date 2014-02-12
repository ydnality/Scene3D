%% oi= s3dRenderOI(fname)
% fname: file name of the pbrt file to render
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
% Todo: We are considering copying the pbrt file that was used to generate
% the scene when we save the oi for future use.
function oi = s3dRenderOI(inputPbrt, focalLength, sceneName)

    if (ieNotDefined('focalLength'))
        focalLength = .050;
    end

    fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
    inputPbrt.writeFile(fullfname);
    
    
    generatedDir = fullfile(dataPath, 'generatedPbrtFiles', 'tempPBRT');
    if exist(generatedDir,'dir')
        unix(['rm ' generatedDir]);
    else
        mkdir(generatedDir);
    end
    outfile  = fullfile(generatedDir, 'temp_out.dat');

    % dMapFile = 'temp_out_DM.dat'; 
    unix([fullfile(pbrtHome, '/src/bin/pbrt ') fullfname ' --outfile ' outfile])

    % ISET will read the PBRT output
    oi = pbrt2oi(outfile);
    %rename the oi, if a name is given
    if (~ieNotDefined('oiName'))
        oi = oiSet(oi, 'name', oiName);
    end
    oi = s3dFixOi(oi, focalLength);
%     chdir('..');
end