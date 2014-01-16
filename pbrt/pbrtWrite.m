function pbrtWrite(pbrt,fname)
% Write a text file from a pbrt structure
%
% 
% Example:
%    pbrt = pbrtCreate;
%    pbrtWrite(pbrt,'deleteMe.pbrt');
%
% AL VISTASOFT Team, Copyright 2013

%%

fid = fopen(fname,'w');

fprintf(fid,'# PBRT v2.0 Blender Scene file (written from Scene3D)\n#\n\n');

%% Scale
% s = pbrtGet(pbrt,'Scale');
fprintf(fid,'Scale ');
fprintf(fid,'%f  ',pbrt.scale);
fprintf(fid,'# account for fixed lookat bug\n');


%% Camera position

fprintf(fid,'\n\nLookAt\n');
fprintf(fid,'\t%f %f %f\n',pbrt.cameraPosition');

fclose(fid);

%% 
end
