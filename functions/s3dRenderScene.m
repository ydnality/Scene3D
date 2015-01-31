function output = s3dRenderScene(inputPbrt, name, noScale, dockerFlag, oiFlag)
% scene = s3dRenderScene(inputPbrt, sceneName, noScale, dockerFlag) Uses a
% pinhole camera model to calculate scene radiance
%
% inputPbrt: either (1) a pbrtObject, or (2) full path (file name) of the
% pbrt file that we render path.  **All includes and referenced files in
% the pbrt file should be of files in the same directory as that pbrt file.
%
% This function renders a PBRT scene using the given a full fname, then
% returns the data as an ISET scene. 
% The pbrt file is generated in a temporary directory. The output file is
% placed in that same directory. The proper optics is also placed into this
% optical image that corresponds to the focal length and field of view
% (FOV).

% Todo: figure out consistency between s3dRenderScen and s3dRenderOi

%% Set up parameters

if ieNotDefined('dockerFlag'), dockerFlag = 0; end
if (ieNotDefined('inputPbrt')), error('PBRT full file name or pbrtObject required.');end
if (ieNotDefined('name')), name = 'default'; end
if (ieNotDefined('noScale')), noScale = false; end
if (ieNotDefined('oiFlag')), oiFlag = false; end

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
            
            [directory, ~, ~] = fileparts(inputPbrt.materialArray{ii});
            inputPbrt.materialArray{ii} = [fName,extension];
            
             %TODO: make into helper function
            [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir); 
            [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.exr'), generatedDir);  %image textures
            [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)
        end
    end
    %do the same thing for the lights
    for ii = 1:numel(inputPbrt.lightSourceArray)
        if ischar(inputPbrt.lightSourceArray{ii}) && exist(inputPbrt.lightSourceArray{ii},'file')
            copyfile(inputPbrt.lightSourceArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.lightSourceArray{ii});
            
             [directory, ~, ~] = fileparts(inputPbrt.lightSourceArray{ii});
            inputPbrt.lightSourceArray{ii} = [fName,extension];
            
             %TODO: make into helper function
            [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir); 
            [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.exr'), generatedDir);  %image textures
            [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)

        end
    end
   %do the same thing for geometry
    for ii = 1:numel(inputPbrt.geometryArray)
        if ischar(inputPbrt.geometryArray{ii}) && exist(inputPbrt.geometryArray{ii},'file')
            copyfile(inputPbrt.geometryArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.geometryArray{ii});
            
            [directory, ~, ~] = fileparts(inputPbrt.geometryArray{ii});
            inputPbrt.geometryArray{ii} = [fName,extension];

            %TODO: make into helper function
            [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir); 
            [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.exr'), generatedDir);  %image textures
            [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
            [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)

        end
    end
    fullfname = fullfile(dataPath, 'generatedPbrtFiles', [name '.pbrt']);
    inputPbrt.writeFile(fullfname);
elseif (ischar(inputPbrt))
    if (~exist(inputPbrt,'file'))
        error('PBRT full file name required.  File not found');
    end
    %if inputPbrt is a char, then it becomes the input file
    fullfname = inputPbrt;

    %copy all relavent files into the temp directory
    [directory, ~, ~] = fileparts(fullfname);
    
    [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.exr'), generatedDir);  %image textures
    [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)
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
    
   % [s,pbrtExe] = system('which pbrt');  %this doesn't work.  not sure
   % why.
    pbrtExe = 'pbrt';
   %if s, error('PBRT not found'); end
    
    % [p,n,ext] = fileparts(fullfname);
    %cmd = sprintf('%s %s --outfile %s\n',pbrtExe,fullfname,outfile);
    [~,n,e] = fileparts(fullfname); % Get name of pbrt input file
    tempInputFile = fullfile(generatedDir, [n e]);
    cmd = sprintf('%s %s --outfile %s\n',pbrtExe,tempInputFile,outfile);
    % chdir(p)
    unix(cmd)
end


%% ISET will read the PBRT output and convert to a scene

if (oiFlag)
    % ISET will read the PBRT output
    output = pbrt2oi(outfile);
    %rename the oi
    output = oiSet(output, 'name', name);
else
    output = pbrt2scene(outfile);
    %rename the scene, if a name is given
    if (~ieNotDefined('name'))
        output = sceneSet(output, 'name', name);
        if(~noScale)
            % By default we scale to a mean luminance of 100 cd/m2.
            output = sceneAdjustLuminance(output,100);
        end
    end
end


end