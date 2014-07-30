%% oi= s3dRenderOI(fname)
% sceneName: file name of the pbrt file to render (without the .pbrt
% extension).  sceneName will also be used to name the optical image.
%
% This function renders a PBRT oi using the given pbrt object, then
% returns the data as an optical image. The pbrt file should be previously 
% generated and be in datapath/generatedPbrtFiles/.  The output file is placed
% in datapath/generatedPbrtFiles/tempPBRT/.  The
% proper optics is also placed into this optical image that corresponds to
% the focal length and field of view (FOV).
%
% Todo: We are considering copying the pbrt file that was used to generate
% the scene when we save the oi for future use.
% 
% Todo: backwards compatibility to render pbrt files

function oi = s3dRenderOI(inputPbrt, focalLength, sceneName)

    if (ieNotDefined('focalLength'))
        focalLength = .050;
    end
    if (ieNotDefined('sceneName'))
        sceneName = 'deleteMe';
    end
    
    if(isa(inputPbrt, 'pbrtObject'))
        fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
        inputPbrt.writeFile(fullfname);
    elseif (ischar(inputPbrt))
        %if inputPbrt is a char, then it becomes the input file
        fullfname = inputPbrt;  
    else
        error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject');
    end
    
    generatedDir = fullfile(dataPath, 'generatedPbrtFiles', 'tempPBRT');
    if exist(generatedDir,'dir')
        unix(['rm ' generatedDir]);
    else
        mkdir(generatedDir);
    end
    outfile  = fullfile(generatedDir, 'temp_out.dat');

    % dMapFile = 'temp_out_DM.dat'; 
    unix(['pbrt ' fullfname ' --outfile ' outfile])

    % ISET will read the PBRT output
    oi = pbrt2oi(outfile);
    %rename the oi, if a name is given
    if (~ieNotDefined('oiName'))
        oi = oiSet(oi, 'name', oiName);
    end
    oi = oiSet(oi, 'name', sceneName);
    oi = s3dFixOi(oi, focalLength);
%     chdir('..');
end