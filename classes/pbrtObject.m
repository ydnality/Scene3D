% pbrtObject contains all the basic properties of a Scene3D pbrt scene.
% All this information will be used to create a pbrt File.  However, note
% that key properties such as materials and geometry are still handled by
% importing text files, since these are very complex.
%TODO: check types
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
        materialArray;  %we will use .fileName for now to avoid declaring all the materials stuff
        geometryArray;   %we will use .fileName for now to avoid declaring all the geometry stuff
    end
    methods
        % default constructor
        function obj = pbrtObject()
            obj.name = 'default';
            
            % This is some blender to pbrt coordinate frames (AL)
            obj.scale = [-1 1 1];
            obj.camera = cameraObject(); %assign the default camera
            
            % Sampler
            obj.sampler = samplerObject(); 
            
            % SurfaceIntegrator
            %             obj.surfaceIntegrator.type  = 'surfaceIntegrator';
            %             obj.surfaceIntegrator.surfIntType = 'directlighting';
            %             obj.surfaceIntegrator.maxdepth = 0;
            obj.surfaceIntegrator = surfaceIntegratorObject(); 

            % Renderer
            obj.renderer = 'sample';
            
            % World Begin
            
            %  Attribute Begin
            %   Light source
            obj.lightSourceArray = cell(1,1);
            whiteLight = lightObject();
            
            obj.lightSourceArray{1} = whiteLight;
            %  Attribute End
            
            %  Material file
            obj.materialArray = cell(1,1);
            obj.materialArray{1} = 'depthTargetSpheres-mat.pbrt';
            
            %example materials object
            tempPropertyArray = cell(1,1);
            tempKd.type = 'color Kd';   %TODO: the user shouldn't need to know pbrt syntax...
            tempKd.value = [0 0.374624 0];
            tempPropertyArray{1} = tempKd;
            greenLambertian = materialObject('greenLambertian', 'matte', tempPropertyArray);
            
            obj.materialArray{2} = greenLambertian;
            
            % Geometry file
            obj.geometryArray = cell(1,1);
            obj.geometryArray{1} = 'depthTargetSpheres-geom.pbrt';
           
            %example geometry object
            %             examplePlane.type = 'geometry';
            %             examplePlane.name = 'Example Plane';
            %             examplePlane.material = 'greenLambertian';
            %             examplePlane.triangleMesh = [0 1 2;
            %                 0 2 3];
            %             examplePlane.points = [1.000000 1.000000 0.000000;
            %                 -1.000000 1.000000 0.000000;
            %                 -1.000000 -1.000000 0.000000;
            %                 1.000000 -1.000000 0.000000];
            %             examplePlane.transform = [-4.37113882867e-008 -0.0 -1.0 0.0;
            %                 1.0 -4.37113882867e-008 -4.37113882867e-008 0.0 ;
            %                 -4.37113882867e-008 -1.0 1.91068567692e-015 0.0;
            %                 0.0 0.0 8.87306690216 1.0];
            
            examplePlane = geometryObject();
            obj.geometryArray{2} = examplePlane;
            
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
            fprintf(fid,'\n\nInclude "%s"\n', obj.camera.lens);
            
            %% Image resolution
            
            fprintf(fid,'\n\nFilm "image"\n');
            fprintf(fid,'\t"integer xresolution" [%i]\n',obj.camera.film.xresolution);
            fprintf(fid,'\t"integer yresolution" [%i]\n',obj.camera.film.yresolution);
            
            %% Sampler
            
            fprintf(fid,'\n\nSampler "%s"\n', obj.sampler.type);    %TODO: consider putting this inside each of the objects? i'm not sure yet
            fprintf(fid,'\t"integer pixelsamples" [%i]\n',obj.sampler.pixelSamples);
            
            %% SurfaceIntegrator
            
            fprintf(fid,'\n\nSurfaceIntegrator "%s"\n', obj.surfaceIntegrator.type);
            fprintf(fid,'\t"integer maxdepth" [%i]\n',obj.surfaceIntegrator.maxDepth);
            
            %% Renderer
            
            fprintf(fid,'\n\nRenderer "%s"\n', obj.renderer);
            %% World Begin
            
            fprintf(fid,'\n\nWorldBegin\n');
            
            %% Lightsource
            
            %loop through each light source
            for i = 1:length(obj.lightSourceArray)
                
                %check to see if a file or directly defined
                %if it is a filename, use an include, if not
                if (isa(obj.lightSourceArray{i}, 'fileObject'))
                    if ~strcmp(obj.lightSourceArray{i}.file , '')
                        fprintf(fid,'\n\nInclude "%s"\n', obj.lightSourceArray{i}.fileName);
                    end
                else
                    %this is the case where the lightsource is explicitly declared
                    if (isa(obj.lightSourceArray{i}, 'lightObject'))
                        fprintf(fid,'\n\nAttributeBegin\n');
                        fprintf(fid,'\n\tLightSource\n');
                        fprintf(fid,'\t\t"%s"\n', obj.lightSourceArray{i}.type);

                        fprintf(fid,'\t\t"%s" [', obj.lightSourceArray{i}.spectrum.type);
                        fprintf(fid,'%f ', obj.lightSourceArray{i}.spectrum.value);
                        fprintf(fid,']\n');

                        fprintf(fid,'\t\t"float coneangle" %f\n', obj.lightSourceArray{i}.coneAngle);
                        fprintf(fid,'\t\t"float conedeltaangle" %f\n', obj.lightSourceArray{i}.coneDeltaAngle);

                        fprintf(fid,'\t\t"point from" [');
                        fprintf(fid,'%f ', obj.lightSourceArray{i}.from);
                        fprintf(fid,']\n');

                        fprintf(fid,'\t\t"point to" [');
                        fprintf(fid,'%f ', obj.lightSourceArray{i}.to);
                        fprintf(fid,']\n');

                        fprintf(fid,'\nAttributeEnd\n');
                    else
                        disp('Error! Non-light type placed in light array!');
                        return;
                    end
                end
            end
            
            %% Materials File
            %TODO check for file or not
            
            for i = 1:length(obj.materialArray)
                if (isa(obj.materialArray{i},'materialObject')); 
                    fprintf(fid,'\n\nMakeNamedMaterial "%s"\n', obj.materialArray{i}.name);
                    fprintf(fid,'\t"string type" ["%s"]\n', obj.materialArray{i}.type);
                    
                    for j = 1:length(obj.materialArray{i}.propertyArray)
                        fprintf(fid,'\t"%s" [', obj.materialArray{i}.propertyArray{j}.type);
                        fprintf(fid,'%f ', obj.materialArray{i}.propertyArray{j}.value);
                        fprintf(fid,']\n');
                    end
                    
                    %                     if (strcmp(obj.materialArray{i}.matType, 'matte'))
                    %                         fprintf(fid,'\t"%s" [', obj.materialArray{i}.kd.kdType);
                    %                         fprintf(fid,'%f ', obj.materialArray{i}.kd.value);
                    %                         fprintf(fid,']\n');
                    %                     end
                else
                    fprintf(fid,'\n\nInclude "%s"\n', obj.materialArray{i});
                end
            end
            %% Geometry File
            %TODO check for file or not
            
            
            for i = 1:length(obj.geometryArray)
                
                curGeometry =obj.geometryArray{i}; 
                if (isa(curGeometry,'geometryObject')); 
                    fprintf(fid,'\n\nAttributeBegin #%s\n', curGeometry.name);
                    
                    fprintf(fid,'\n\tTransform \n\t[\n');
                    fprintf(fid,'\t%f %f %f %f \n', curGeometry.transform' );
                    fprintf(fid,'\t]\n');
                    
                    fprintf(fid,'\tNamedMaterial "%s"\n', curGeometry.material);
                    
                    fprintf(fid,'\tShape "trianglemesh" "integer indices" \n\t[\n');
                    fprintf(fid,'\t%i %i %i\n', curGeometry.triangleMesh');
                    fprintf(fid,'\t]\n');
                    
                    fprintf(fid,'\t"point P" \n\t[\n');
                    fprintf(fid,'\t%f %f %f\n', curGeometry.points');
                    fprintf(fid,'\t]\n');                   
                    
                    fprintf(fid,'\n\nAttributeEnd\n');
                else
                    fprintf(fid,'\n\nInclude "%s"\n', curGeometry);
                end
            end
            %% World End
            fprintf(fid,'\n\nWorldEnd\n');
            
            %%
            fclose(fid);
            
            returnVal = 1; %1 for success, 0 for failure
        end
        
        function addLightSource(obj,newLightSource)
            arrayLength = length(obj.lightSourceArray);
            obj.lightSourceArray{arrayLength+1} = newLightSource; 
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

