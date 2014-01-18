% pbrtObject contains all the basic properties of a Scene3D pbrt scene.
% All this information will be used to create a pbrt File.  However, note
% that key properties such as materials and geometry are still handled by
% importing text files, since these are very complex.
classdef pbrtObject <  handle
    
    %matlab.mixin.Copyable
    properties (SetAccess = private)
        name;
        scale;
        camera; %this will contain position, file, film
        sampler;
        surfaceIntegrator;
        renderer;
        lightSourceArray;
        materials;  %we will use .fileName for now to avoid declaring all the materials stuff
        geometry;   %we will use .fileName for now to avoid declaring all the geometry stuff
    end
    methods
        % default constructor
        function obj = pbrtObject()
            obj.name = 'default';
            
            % This is some blender to pbrt coordinate frames (AL)
            obj.scale = [-1 1 1];
            
            % Placement of the camera.  Defined by a vector that looks in a certain
            % direction, and then tells you which way is up.
            obj.camera.position = [
                4.5 -80 7 % Starting up
                4.5 -79 7 % Ending up
                0 0 1];   % Which way is up
            
            % Example lens
            obj.camera.lens.file = 'idealLens-50mm.pbrt';
            
            % Typical "sensor"
            obj.camera.film.name = 'image';
            obj.camera.film.xresolution = 200;
            obj.camera.film.yresolution = 200;
            
            % Sampler
            obj.sampler.type = 'lowdiscrepancy';
            obj.sampler.pixelsamples = 512;
            
            % SurfaceIntegrator
            obj.surfaceIntegrator.type  = 'directlighting';
            obj.surfaceIntegrator.maxdepth = 0;
            
            % Renderer
            obj.renderer = 'sample';
            
            % World Begin
            
            %  Attribute Begin
            %   Light source
            obj.lightSourceArray = cell(1,1);
            
            %default white light at camera position
            whiteLight.type = 'spot';
            whiteLight.spectrum.type = 'rgb I';
            whiteLight.spectrum.value = [1000 1000 1000];
            whiteLight.coneangle = 180;
            whiteLight.conedeltaangle = 180;
            whiteLight.from = [4.5 -90 8.5];
            whiteLight.to = [4.5 -89 8.5];
            
            obj.lightSourceArray{1} = whiteLight;
            %  Attribute End
            
            %  Material file
            obj.materials.file = 'depthTargetDepths-mat.pbrt';
            %  Geometry file
            obj.geometry.file = 'depthTargetDepths-geom.pbrt';
            % WorldEnd
        end
        
        % Write a text file from a pbrt structure
        function returnVal = writeFile(obj, fname)
            fid = fopen(fname,'w');
            fprintf(fid,'# PBRT v2.0 Blender Scene file (written from Scene3D)\n#\n\n');
            
            %% Scale
            % s = pbrtGet(pbrt,'Scale');
            fprintf(fid,'Scale ');
            fprintf(fid,'%f  ',obj.scale);
            fprintf(fid,'# account for fixed lookat bug\n');
            
            %% Camera position
            
            fprintf(fid,'\n\nLookAt\n');
            
            if (isfield(obj.camera.position, 'fileName'))
                fprintf(fid,'\n\nInclude "%s"\n', obj.camera.position.fileName);
            else
                fprintf(fid,'\t%f %f %f\n',obj.camera.position');
            end
            
            %% Lens file 
            % TODO check for lens file or not a lens file
            fprintf(fid,'\n\nInclude "%s"\n', obj.camera.lens.file);
            
            %% Image resolution
            
            fprintf(fid,'\n\nFilm "image"\n');
            fprintf(fid,'\t"integer xresolution" [%i]\n',obj.camera.film.xresolution);
            fprintf(fid,'\t"integer yresolution" [%i]\n',obj.camera.film.yresolution);
            
            %% Sampler
            
            fprintf(fid,'\n\nSampler "%s"\n', obj.sampler.type);
            fprintf(fid,'\t"integer pixelsamples" [%i]\n',obj.sampler.pixelsamples);
            
            %% SurfaceIntegrator
            
            fprintf(fid,'\n\nSurfaceIntegrator "%s"\n', obj.surfaceIntegrator.type);
            fprintf(fid,'\t"integer maxdepth" [%i]\n',obj.surfaceIntegrator.maxdepth);
            
            %% Renderer
            
            fprintf(fid,'\n\nRenderer "%s"\n', obj.renderer);
            %% World Begin
            
            fprintf(fid,'\n\nWorldBegin\n');
            
            %% Lightsource
            
            %loop through each light source
            for i = 1:length(obj.lightSourceArray)
                
                %check to see if a file or directly defined
                %if it is a filename, use an include, if not
                if (isfield(obj.lightSourceArray{i}, 'fileName'))
                    if ~strcmp(obj.lightSourceArray{i}.fileName , '')
                        fprintf(fid,'\n\nInclude "%s"\n', obj.lightSourceArray{i}.fileName);
                    end
                else
                    %this is the case where the lightsource is explicitly declared
                    fprintf(fid,'\n\nAttributeBegin\n');
                    fprintf(fid,'\n\tLightSource\n');
                    fprintf(fid,'\t\t"%s"\n', obj.lightSourceArray{i}.type);
                    
                    fprintf(fid,'\t\t"%s" [', obj.lightSourceArray{i}.spectrum.type);
                    fprintf(fid,'%f ', obj.lightSourceArray{i}.spectrum.value);
                    fprintf(fid,']\n');
                    
                    fprintf(fid,'\t\t"float coneangle" %f\n', obj.lightSourceArray{i}.coneangle);
                    fprintf(fid,'\t\t"float conedeltaangle" %f\n', obj.lightSourceArray{i}.conedeltaangle);
                    
                    fprintf(fid,'\t\t"point from" [');
                    fprintf(fid,'%f ', obj.lightSourceArray{i}.from);
                    fprintf(fid,']\n');
                    
                    fprintf(fid,'\t\t"point to" [');
                    fprintf(fid,'%f ', obj.lightSourceArray{i}.to);
                    fprintf(fid,']\n');
                    
                    fprintf(fid,'\nAttributeEnd\n');
                end
            end
            
            %% Materials File
            %TODO check for file or not
            fprintf(fid,'\n\nInclude "%s"\n', obj.materials.file);
            
            %% Geometry File
            %TODO check for file or not
            fprintf(fid,'\n\nInclude "%s"\n', obj.geometry.file);
            
            %% World End
            fprintf(fid,'\n\nWorldEnd\n');
            
            %%
            fclose(fid);
            
            returnVal = 1; %1 for success, 0 for failure
        end
        
        %example code
        %     function obj = batchFileClass(inStem, inPostfix)
        %         obj.inputStem = inStem;
        %         obj.postfix = inPostfix;
        %         obj.outputStem = '';
        %     end
        %
        %     function obj=setInputStem(obj, inStem)
        %         obj.inputStem = inStem;
        %
        %     end
        %
        %     function obj=setOutputStem(obj, outStem)
        %         obj.outputStem = outStem;
        %     end
        %
        %     function inputFile = getInputFile(obj)
        %         inputFile = strcat(obj.inputStem, obj.postfix);
        %     end
        %
        %     function outString = getOutputFile(obj)
        %         outString = strcat(obj.outputStem, obj.postfix);
        %     end
        %
        %     function postFix = getPostFix(obj)
        %         postFix = obj.postfix;
        %     end
    end
end

