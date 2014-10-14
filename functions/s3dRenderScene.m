function scene = s3dRenderScene(inputPbrt, sceneName)
% scene = s3dRenderScene(inputPbrt, sceneName)
% Uses a pinhole camera model to calculate scene radiance
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
%
% Note: pbrt looks for imported image and data files at
% s3dRootPath/data/generatedPbrtFiles
%%  


    if (ieNotDefined('inputPbrt'))
        error('PBRT full file name required.  File not found');
    end

    if (ieNotDefined('sceneName'))
        sceneName = 'deleteMe';
    end
    
    if(isa(inputPbrt, 'pbrtObject'))
        fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
        inputPbrt.writeFile(fullfname);
    elseif (ischar(inputPbrt))
        if (~exist(inputPbrt,'file'))
            error('PBRT full file name required.  File not found');
        end
        %if inputPbrt is a char, then it becomes the input file
        fullfname = inputPbrt;
    else
        error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject');
    end
    
    
    % Use pinhole and pbrt to create the scene data
    %pbrtExe = fullfile(pbrtRootPath, 'src','bin','pbrt');
    pbrtExe = 'pbrt';
    
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