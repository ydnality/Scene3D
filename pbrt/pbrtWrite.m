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

%% Lens file
fprintf(fid,'\n\nInclude "%s"\n', pbrt.lensFile);

%% Image resolution

fprintf(fid,'\n\nFilm "image"\n');
fprintf(fid,'\t"integer xresolution" [%i]\n',pbrt.film.xresolution);
fprintf(fid,'\t"integer yresolution" [%i]\n',pbrt.film.yresolution);

%% Sampler 

fprintf(fid,'\n\nSampler "%s"\n', pbrt.sampler.type);
fprintf(fid,'\t"integer pixelsamples" [%i]\n',pbrt.sampler.pixelsamples);

%% SurfaceIntegrator

fprintf(fid,'\n\nSurfaceIntegrator "%s"\n', pbrt.surfaceIntegrator.type);
fprintf(fid,'\t"integer maxdepth" [%i]\n',pbrt.surfaceIntegrator.maxdepth);

%% Renderer

fprintf(fid,'\n\nRenderer "%s"\n', pbrt.renderer);
%% World Begin

fprintf(fid,'\n\nWorldBegin\n');

%% Lightsource

%loop through each light source
for i = 1:length(pbrt.lightSource)
    fprintf(fid,'\n\nAttributeBegin\n');
    fprintf(fid,'\n\tLightSource\n');
    fprintf(fid,'\t\t"%s"\n', pbrt.lightSource{i}.type);
    
    fprintf(fid,'\t\t"%s" [', pbrt.lightSource{i}.spectrum.type);
    fprintf(fid,'%f ', pbrt.lightSource{i}.spectrum.value);
    fprintf(fid,']\n');
    
    fprintf(fid,'\t\t"float coneangle" %f\n', pbrt.lightSource{i}.coneangle);
    fprintf(fid,'\t\t"float conedeltaangle" %f\n', pbrt.lightSource{i}.conedeltaangle);
    
    fprintf(fid,'\t\t"point from" [');
    fprintf(fid,'%f ', pbrt.lightSource{i}.from);
    fprintf(fid,']\n');
    
    fprintf(fid,'\t\t"point to" [');
    fprintf(fid,'%f ', pbrt.lightSource{i}.to);
    fprintf(fid,']\n');    
    
    fprintf(fid,'\nAttributeEnd\n');
end

%% Materials File
fprintf(fid,'\n\nInclude "%s"\n', pbrt.materialsFile);

%% Geometry File
fprintf(fid,'\n\nInclude "%s"\n', pbrt.geometryFile);

%% World End
fprintf(fid,'\n\nWorldEnd\n');

%%
fclose(fid);
end
