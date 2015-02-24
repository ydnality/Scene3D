%% output = s3dRenderDepthMap(fname)
%
%s3dRenderDepthMap(inputPbrt, numRenders, sceneName, dockerFlag)
%
% Accepts 2 types of input for inputPbrt: 
% (1) pbrtObject: This is the preferred type, as the function will change
% the sampler to the correct one, and should run as is.
% (2) pbrt text file. 
% Returns the ground truth depth map from given fname.  
% This mode requires some manipulation of the .pbrt file before execution!
% ****Note that the number of pixel samples here must be set to 1,
% and the number of reflections must be set to 0 for a correct output to be
% produced!


function output = s3dRenderDepthMap(inputPbrtIn, numRenders, sceneName, dockerFlag)

if ieNotDefined('dockerFlag'), dockerFlag = 0; end

if (ieNotDefined('numRenders'))
    numRenders = 1;
end

if (ieNotDefined('sceneName'))
    sceneName = 'deleteMe';
end

%make a temporary directory
generatedDir = tempname();
mkdir(generatedDir);

% Make the pbrt file from the pbrtObject
if(isa(inputPbrtIn, 'pbrtObject'))
    % first, make sure the sampler is of type stratified, so that we get a
    % clean depth map.  Also, set samples to 1, and jitter off to eliminate
    % noise.
    
    %clone pbrt object first
    inputPbrt = pbrtObject();
    inputPbrt.makeDeepCopy(inputPbrtIn);
    
    inputPbrt.sampler.setType('stratified');
    for i = 1:length(inputPbrt.sampler.propertyArray)
        inputPbrt.sampler.removeProperty();
    end
    inputPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '1'));
    inputPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '1'));
    inputPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

    % Strip the path off of each materialArray entry before writing the
    % pbrt file to disk.
%     for ii = 1:numel(inputPbrt.materialArray)
%         if ischar(inputPbrt.materialArray{ii}) && exist(inputPbrt.materialArray{ii},'file')
%             copyfile(inputPbrt.materialArray{ii},generatedDir)
%             [~,fName,extension] = fileparts(inputPbrt.materialArray{ii});
%             inputPbrt.materialArray{ii} = [fName,extension];
%         end
%         
%     end
%     
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
            [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)
            [status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), generatedDir);   %copies all .dat files (lens files)

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
            [status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), generatedDir);   %copies all .dat files (lens files)

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
            [status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), generatedDir);   %copies all .dat files (lens files)

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
    
    fullfname = fullfile(dataPath, 'generatedPbrtFiles', [sceneName '.pbrt']);
    inputPbrt.writeFile(fullfname);
elseif (ischar(inputPbrtIn))
    pbrtExe = 'pbrt';
    
    if ~exist(pbrtExe,'file')
        error('PBRT executable not found');
    end
    if (~exist(inputPbrtIn,'file'))
        error('PBRT full file name required.  File not found');
    end
    %if inputPbrt is a char, then it becomes the input file
    fullfname = inputPbrtIn;

    %copy all relavent files into the temp directory
    [directory, ~, ~] = fileparts(fullfname);
    
    [status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.tga'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), generatedDir);
    [status,message,messageid] = copyfile(fullfile(directory, '*.dat'), generatedDir);   %copies all .dat files (lens files)
    [status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), generatedDir);   %copies all .dat files (lens files)

else
    error('invalid inputPbrt type.  Must be either a character array of the pbrt file, or a pbrtObject');
end


copyfile(fullfname,generatedDir);
outfile  = fullfile(generatedDir, 'temp_out.dat');
dMapFile  = fullfile(generatedDir, 'temp_out_DM.dat');


%renders a few times then takes the median.  This median operation is
%only necessary if not using stratified pbrt sampling.  The default of
%numRenders is 1, and this should be what we use 99% of the time.
for i = 1:numRenders
    if (dockerFlag)
        % We assume docker is installed on this system and we execute the
        % function in a docker container
        s = system('which docker');
        if s
            warning('Docker not found! Trying on system pbrt instead!');
            s3dPbrtLocalCall(fullfname, generatedDir, outfile);
        else
            dHub = 'vistalab/pbrt:spherical';  % Docker container at dockerhub
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
        end
    else
        warning('Docker not found! Trying on system pbrt instead!');
        s3dPbrtLocalCall(fullfname, generatedDir, outfile);
    end
    if (i ==1)
        oi = pbrt2oi(outfile);
        imageHeight = oiGet(oi, 'rows');
        imageWidth = oiGet(oi, 'cols');
        depthMap = zeros(imageHeight, imageWidth, numRenders);
    end
    depthMap(:,:, i) = s3dReadDepthMapFile(dMapFile, [imageHeight imageWidth]);
end

depthMapProcessedMedian = median(depthMap, 3);
output = depthMapProcessedMedian;
end