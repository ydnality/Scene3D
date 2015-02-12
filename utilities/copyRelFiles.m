function [ output_args ] = copyRelFiles( directory, newDir )
%copies all relavent pbrt files from one directory, to another directory
%   Detailed explanation goes here

[status,message,messageid] = copyfile(fullfile(directory, '*.pbrt'), newDir); 
[status,message,messageid] = copyfile(fullfile(directory, '*.tga'), newDir);
[status,message,messageid] = copyfile(fullfile(directory, '*.exr'), newDir);  %image textures
[status,message,messageid] = copyfile(fullfile(directory, '*.jpg'), newDir);
[status,message,messageid] = copyfile(fullfile(directory, '*.dat'), newDir);   %copies all .dat files (lens files)
[status,message,messageid] = copyfile(fullfile(directory, '*.brdf'), newDir);   %copies all .dat files (lens files)

output_args = status;

end

