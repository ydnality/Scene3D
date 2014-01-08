%% returnVal = s3dReplaceVar(pbrtFile, outFile, varArray)
% Replaces variables in template pbrt files with proper values, given by conditions file
% pbrtFile: the template pbrt file, with <variableName> as the general form for variables 
% outFile: the output file name
% varArray: a cell array containing the variable values 
% TODO: write more documentation on this function 
% TODO: make varArray an array of a new type pbrtVars as a better programming practice
% **adapted from code originally written by my cousin Yeong-Jeh(James) Huang

function returnVal = s3dReplaceVar(pbrtFile, outFile, varArray)
    
    %we need to figure out where exactly to place the new temp file
%     copyfile(pbrtFile, outFile);
       
    [fid msg] = fopen(pbrtFile,'r'); % open original .pbrt file
    if ~isempty(msg)
        returnVal = 1;
        disp(['Template file, ' pbrtFile ' not found!']);
        return;
    end

%     mkdir ('batchOutput');
    nfid = fopen(outFile,'w');
    while ~feof(fid)
         tline = fgets(fid); % fgetl() will read row by row, and will not ignore enpty row
         if~ischar(tline),break,end % check if success
         %disp(tline);
         
         for i = 1:length(varArray)
            tline = regexprep(tline, strcat('<', varArray{i}.name, '>'), varArray{i}.value);
         end
         
         fprintf(nfid,'%s',tline);
         
         %xresfind = strfind(tline,'xresolution'); % find keyword "xresolution"
         %yresfind = strfind(tline,'yresolution'); % find keyword "yresolution"
         
    end
    fclose(nfid);
    fclose(fid); % close .pbrt file

    returnVal = 0;
    return;
end