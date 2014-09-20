%% output = s3dRenderDepthMap(fname)
%
% Returns the ground truth depth map from given fname 
% (relative to s3droot/scripts/pbrtFiles/)
% fname must correspond to the pbrt file whose scene's depth map will be
% rendered.  **Note that the number of pixel samples here must be set to 1,
% and the number of reflections must be set to 0.

% This function works by rendering lots of 1 sample pbrt scenes, and taking
% the median value of those rendered depth maps, to create one that is of
% high quality.  There is some discussion of whether or not this produces a
% valid depth map, but for all intensive purposes, it is close enough.
function output = s3dRenderDepthMap(fname, numRenders)

    if (ieNotDefined('numRenders'))
        numRenders = 31;
    end
%     chdir(fullfile(s3dRootPath, 'scripts', 'pbrtFiles'));
    mkdir('tempOutput');
    chdir('tempOutput');
    unix('rm -rf*');
    
    outfile = 'depthRender_out.dat';
    dMapFile = 'depthRender_out_DM.dat';
    imageWidth = 300;
    imageHeight = 200;
    depthMap = zeros(imageHeight, imageWidth, numRenders);

    for i = 1:numRenders
        unix([fullfile(pbrtHome, '/src/bin/pbrt ') '../' fname ' --outfile ' outfile]);
        if (i ==1)
            oi = pbrt2oi(outfile);
            imageHeight = oiGet(oi, 'rows');
            imageWidth = oiGet(oi, 'cols');
        end
        depthMap(:,:, i) = s3dReadDepthMapFile(dMapFile, [imageHeight imageWidth]);
        unix('rm *');
    end

    depthMapProcessedMedian = median(depthMap, 3);
    output = depthMapProcessedMedian;
end