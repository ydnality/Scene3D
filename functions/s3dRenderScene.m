function scene = s3dRenderScene(inputPbrt, sceneName, noScale, dockerFlag)
% scene = s3dRenderScene(inputPbrt, sceneName, noScale, dockerFlag)
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
%% Set up parameters

if ieNotDefined('dockerFlag'), dockerFlag = 0; end
if (ieNotDefined('inputPbrt')), error('PBRT full file name required.');end
if (ieNotDefined('sceneName')), sceneName = 'deleteMe'; end
if (ieNotDefined('noScale')), noScale = false; end

%% Make a tempPBRT directory where the output files will go
generatedDir = tempname;
mkdir(generatedDir);

% Make the pbrt files from the pbrtObject
if(isa(inputPbrt, 'pbrtObject'))
    % Strip the path off of each materialArray entry before writing the
    % pbrt file to disk.
    for ii = 1:numel(inputPbrt.materialArray)
        if ischar(inputPbrt.materialArray{ii}) && exist(inputPbrt.materialArray{ii},'file')
            copyfile(inputPbrt.materialArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.materialArray{ii});
            inputPbrt.materialArray{ii} = [fName,extension];
        end
    end
    fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
    inputPbrt.writeFile(fullfname);
elseif (ischar(inputPbrt))
    if (~exist(inputPbrt,'file'))
        error('PBRT full file name required.  File not found');
    end
    %if inputPbrt is a char, then it becomes the input file
    fullfname = inputPbrt;
    if dockerFlag
        %error('inputPbrt as char is not supported. You must pass a pbrt object');
        
        %copy all relavent files into the temp directory
        [directory, ~, ~] = fileparts(fullfname);
        
        [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir); 
        [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
        [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
        [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)
    end
else
    error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject');
end


%copy the base .pbrt file into the temp directory for processing
copyfile(fullfname,generatedDir);
outfile  = fullfile(generatedDir, 'temp_out.dat');

%% Execute either via docker or via local instance of pbrt
if dockerFlag
    % We assume docker is installed on this system and we execute the
    % function in a docker container
    s = system('which docker');
    if s, error('Docker not found'); end
    
    dHub = 'vistalab/pbrt';  % Docker container at dockerhub
    dCommand = 'pbrt';       % Command run in the docker
    
    [~,n,e] = fileparts(fullfname); % Get name of pbrt input file
    [~,outstem,outext] = fileparts(outfile); % Get name of output file
    
    cmd = sprintf('docker run -t -i -v %s:/data %s %s /data/%s --outfile /data/%s',generatedDir,dHub,dCommand,[n,e],[outstem,outext]);
    
    % Execute the docker call
    [s,r] = system(cmd);
    if s, error('Docker execution failure %s\n',r); end
    disp(r);
    
    % Tell the user where the result iss
    fprintf('Wrote: %s',outfile);
    
else
    % We assume pbrt exists in the local file system and we run it
    % here.
    
    % Use pinhole and pbrt to create the scene data
    %pbrtExe = fullfile(pbrtRootPath, 'src','bin','pbrt');
    
    [s,pbrtExe] = system('which pbrt');
    if s, error('PBRT not found'); end
    
    % [p,n,ext] = fileparts(fullfname);
    cmd = sprintf('%s %s --outfile %s\n',pbrtExe,fullfname,outfile);
    % chdir(p)
    unix(cmd)
end


%% ISET will read the PBRT output and convert to a scene

scene = pbrt2scene(outfile);

%rename the scene, if a name is given
if (~ieNotDefined('sceneName'))
    scene = sceneSet(scene, 'name', sceneName);
    if(~noScale)
        % By default we scale to a mean luminance of 100 cd/m2.
        scene = sceneAdjustLuminance(scene,100);
    end
end

end