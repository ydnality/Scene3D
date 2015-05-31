function [ output_args ] = copyRelFiles( directory, newDir )
% Copies all relevant pbrt files from one directory, to another directory
%
% We use this function in s3dRenderScene and s3dRenderOI for placing
% essential data files into the directory for the docker container.
%
% When we make the pbrt directory, sometimes we just specify the core file
% name.  But the other files need to stay with it.  This function moves all
% the files together, including the core (pbrt) file, and any potentially
% associated files that are in the pbrt directory.
%
% AL Vistasoft

[status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), newDir); 
[status,message,messageid] = copyfile(fullfile(directory, '*.tga'), newDir);
[status,message,messageid] = copyfile(fullfile(directory, '*.exr'), newDir);  %image textures
[status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), newDir);
[status,message,messageid] = copyfile(fullfile(directory, '*.dat'), newDir);   %copies all .dat files (lens files)
[status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), newDir);   %copies all .dat files (lens files)

output_args = status;

end

