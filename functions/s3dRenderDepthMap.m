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
function output = s3dRenderDepthMap(inputPbrt, numRenders, sceneName, dockerFlag)

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
if(isa(inputPbrt, 'pbrtObject'))
    % first, make sure the sampler is of type stratified, so that we get a
    % clean depth map.  Also, set samples to 1, and jitter off to eliminate
    % noise.
    inputPbrt.sampler.setType('stratified');
    for i = 1:length(inputPbrt.sampler.propertyArray)
        inputPbrt.sampler.removeProperty();
    end
    inputPbrt.sampler.addProperty(pbrtPropertyObject('integer xsamples', '1'));
    inputPbrt.sampler.addProperty(pbrtPropertyObject('integer ysamples', '1'));
    inputPbrt.sampler.addProperty(pbrtPropertyObject('bool jitter', '"false"'));

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
    pbrtExe = 'pbrt';
    
    if ~exist(pbrtExe,'file')
        error('PBRT executable not found');
    end
    if (~exist(inputPbrt,'file'))
        error('PBRT full file name required.  File not found');
    end
    %if inputPbrt is a char, then it becomes the input file
    fullfname = inputPbrt;
    if dockerFlag
        error('inputPbrt as char is not supported. You must pass a pbrt object');
    end
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
        %cmd = sprintf('%s %s --outfile %s',pbrtExe,fullfname,outfile);
        [s,pbrtExe] = system('which pbrt');
        if s, error('PBRT not found'); end
        
        % [p,n,ext] = fileparts(fullfname);
        cmd = sprintf('%s %s --outfile %s',pbrtExe,fullfname,outfile);
        unix(cmd);
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