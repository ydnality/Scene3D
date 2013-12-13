% renders batch files in PBRT.
% templateFileArray: an array of type batchFileClass that contains the
% template files for this batch process
% conditionsFile: text file containing the conditions (replacement
% conditions) for the batch process
% returns the renderFileArray (files that were rendered in pbrt)
function renderFileArray = s3dRenderBatch(templateFileArray, conditionsFile)

%initialize this variable
clear('renderFileArray');
renderFileArray = {};

% read first line and determine file format
[fid  msg]= fopen(conditionsFile,'r'); % open .pbrt file

if ~feof(fid)
     tline = fgets(fid);
end
delimiterFind = regexp(tline,'\t')
numVars = length(delimiterFind) + 1;
fclose(fid);

%make format string
formatString = '';
for i = 1:numVars
   formatString = strcat(formatString, '%s'); 
end


%import data
fid = fopen(conditionsFile);
import = textscan(fid, formatString, 'delimiter' , '\t');
fclose(fid);

%assign variable names by looking at imported data
%TODO: add error handling for import{1}
numConditons = length(import{1}) - 1;
varArray = cell(1, numVars);
for i = 1:numVars
    temp = import{i};
    varArray{i}.name = temp(1);
end

%replace variable templates in file
%loop through conditions (rows)
for i = 2:numConditons + 1
    i-1
    
    %loop through variables (columns)
    for j = 1:numVars
        %assign data structure for this condition
        tempVarArray = import{j};
        tempVar = tempVarArray(i);
        varArray{j}.value = tempVar{1};  %remove cell array structure
        
        %if the variable is fileName, then set the stem for the template
        %array to the current filename
        if (strcmp(varArray{j}.name, 'fileName'))
            for l = 1:length(templateFileArray)
                templateFileArray{l}.setOutputStem(varArray{j}.value); 
            end
        end
    end
    
    %for this particular condition, replace variables with respective
    %values
    for k = 1:length(templateFileArray)
        s3dReplaceVar(templateFileArray{k}.getInputFile(), templateFileArray{k}.getOutputFile() , varArray); 
        
        if (strcmp(templateFileArray{k}.getPostFix(), '.pbrt'))
           renderFileArray = cat(1, renderFileArray, templateFileArray{k}.getOutputFile());
        end
    end
end

% 
% %run the new pbrt files
for i = 1:length(renderFileArray)
    oi = s3dRenderScene(renderFileArray{i}, .050, pwd);
    oi = oiSet(oi, 'name', renderFileArray{i});

    vcAddAndSelectObject(oi);
    oiWindow;
    opticalimage = oi;
    vcSaveObject(opticalimage, [renderFileArray{i} '.mat']);
    
end

end