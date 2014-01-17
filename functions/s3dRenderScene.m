%% oi= s3dRenderScene(fname)
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
function oi = s3dRenderScene(fname, focalLength, path)

    if (ieNotDefined('focalLength'))
        focalLength = .050;
    end
    if (ieNotDefined('path'))
        % PBRT will run the PBRT script - initializing
        chdir(fullfile(s3dRootPath, 'data', 'pbrtScenes'));
    else
        chdir(path)
    end   
    

    mkdir('tempOutput');
    chdir('tempOutput');
    unix('rm *');

    % scene rendering
    % unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
    outfile = 'temp_out.dat';
    dMapFile = 'temp_out_DM.dat'; 
    unix([fullfile(pbrtHome, '/src/bin/pbrt ') '../' fname ' --outfile ' outfile]);

    % ISET will read the PBRT output
    oi = pbrt2oi(outfile);
    oi = s3dFixOi(oi, focalLength);
    chdir('..');
end