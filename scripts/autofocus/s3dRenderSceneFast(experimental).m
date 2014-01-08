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
function oi = s3dRenderSceneFast(fname)

    unix('rm temp_out.dat');
    % scene rendering
    % unix([fullfile(pbrtHome, '/src/bin/pbrt') fname '--outfile output.dat']);
    outfile = 'temp_out.dat';
    unix([fullfile(pbrtHome, '/src/bin/pbrt ') '../' fname ' --outfile ' outfile]);

    % ISET will read the PBRT output
    oi = pbrt2oi(outfile);
end