function output = s3dRenderScene(inputPbrt, name, noScale, dockerFlag, oiFlag)
% Render a scene using pbrt from a docker container
%
% sscene = s3dRenderScene(inputPbrt, sceneName, noScale, dockerFlag, oiFlag) 
%
% Renders a PBRT scene using the given full fname, then returns the data as
% an ISET scene.
%
% When called with oiFlag true, the output format is an OI.  Otherwise
% the output format is a scene.
% 
% To create scene radiance, the inputPBRT should use a pinhole camera
% model. We should probably enforce this at some point. **** TODO ***
%
%  inputPbrt: either (1) a pbrtObject, or (2) full path (file name) of the
%    pbrt file that we render path.  **All includes and referenced files in
%    the pbrt file should be of files in the same directory as that pbrt file.
%  name:      Name of the scene
%  noScale:   Scales the pbrt values or not.  (true)s
%             Sometimes we do not want to
%             scale, when the relative intensities matter (e.g., two flash
%             estimation).
%  dockerFlag: Use the docker container at vistalab/pbrt (false)
%  oiFlag:     Return on oi or a scene (false)
%
% The pbrt file is generated in a temporary directory. The output file is
% placed in that same directory. The proper optics is also placed into this
% optical image that corresponds to the focal length and field of view
% (FOV).
%
% Example
%   inputPbrt = pbrtObject;
%   scene = s3dRenderScene(inputPbrt,'deleteMe.mat',false, true, false);
%
%  Bad idea.  
%   oi = s3dRenderScene(inputPbrt,'deleteMe.mat',false, true, true);
%
% See also: s3dRenderOi
%
% Todo: figure out consistency between s3dRenderScene and s3dRenderOi

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
            
            copyRelFiles(directory, generatedDir);        
        end
    end
    %do the same thing for the lights
    for ii = 1:numel(inputPbrt.lightSourceArray)
        if ischar(inputPbrt.lightSourceArray{ii}) && exist(inputPbrt.lightSourceArray{ii},'file')
            copyfile(inputPbrt.lightSourceArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.lightSourceArray{ii});
            
            [directory, ~, ~] = fileparts(inputPbrt.lightSourceArray{ii});
            inputPbrt.lightSourceArray{ii} = [fName,extension];
            
            copyRelFiles(directory, generatedDir);
        end
    end
   %do the same thing for geometry
    for ii = 1:numel(inputPbrt.geometryArray)
        if ischar(inputPbrt.geometryArray{ii}) && exist(inputPbrt.geometryArray{ii},'file')
            copyfile(inputPbrt.geometryArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.geometryArray{ii});
            
            [directory, ~, ~] = fileparts(inputPbrt.geometryArray{ii});
            inputPbrt.geometryArray{ii} = [fName,extension];

            copyRelFiles(directory, generatedDir);
        end
    end
    
    %do the same thing for include files
    for ii = 1:numel(inputPbrt.includeArray)
        if ischar(inputPbrt.includeArray{ii}) && exist(inputPbrt.includeArray{ii},'file')
            copyfile(inputPbrt.includeArray{ii},generatedDir)
            [~,fName,extension] = fileparts(inputPbrt.includeArray{ii});
            
            [directory, ~, ~] = fileparts(inputPbrt.includeArray{ii});
            inputPbrt.includeArray{ii} = [fName,extension];

            copyRelFiles(directory, generatedDir);
        end
    end
    
    % Sometimes we need to copy the camera lens file.  This is a draft of
    % how we might approach it.
    if isa(inputPbrt.camera.lens,'pbrtLensRealisticObject')
        lensFileName = inputPbrt.camera.lens.specFile;
        lensFilePath = which(lensFileName);
        lensFileG = fullfile(generatedDir,lensFileName);
        status = copyfile(lensFilePath,lensFileG);
        if ~status, disp('Copying lens spec file seems problematic'); end
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
    
    copyRelFiles(directory, generatedDir);
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
    if s
        %, error('Docker not found'); 
        warning('Docker not found! Trying on system pbrt instead!');
        s3dPbrtLocalCall(fullfname, generatedDir, outfile);
    else
        % Initialize the docker container
        dHub = 'vistalab/pbrt';  % Docker container at dockerhub
        fprintf('Checking for most recent docker container\n');
        system(sprintf('docker pull %s',dHub));
        
        % Start the docker container that runs pbrt
        dCommand = 'pbrt';       % Command run in the dockers
        [~,n,e] = fileparts(fullfname); % Get name of pbrt input file
        [~,outstem,outext] = fileparts(outfile); % Get name of output file

        % rm = clears the container when it is finished running
        % -t = terminal to grab tty output
        % -i = interactive (not sure it's needed)
        cmd = sprintf('docker run -t -i --rm -v %s:/data %s %s /data/%s --outfile /data/%s',generatedDir,dHub,dCommand,[n,e],[outstem,outext]);

        % Execute the docker call
        [s,r] = system(cmd);
        if s, error('Docker execution failure %s\n',r); end
        disp(r);

        % Tell the user where the result iss
        fprintf('Wrote: %s',outfile);
    end
else
    s3dPbrtLocalCall(fullfname, generatedDir, outfile);
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

%remove the temporary directory and all the files
rmdir(generatedDir, 's');

end

