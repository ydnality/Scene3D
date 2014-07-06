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
%         shapeArray; %built in shapes
    end
    methods
        
        function obj = pbrtObject()
        % Default constructor. Creates default camera, sampler, surface
        % integrator, light sources, materials, and geometries for the pbrt
        % scene.
            obj.name = 'default';
            
            % This is some blender to pbrt coordinate frames (AL)
            obj.scale = [-1 1 1];
            obj.camera = pbrtCameraObject(); %assign the default camera
            
            % Sampler
            obj.sampler = pbrtSamplerObject('lowdiscrepancy', pbrtPropertyObject('integer pixelsamples', 128)); 
            
            % SurfaceIntegrator
            %             obj.surfaceIntegrator.type  = 'surfaceIntegrator';
            %             obj.surfaceIntegrator.surfIntType = 'directlighting';
            %             obj.surfaceIntegrator.maxdepth = 0;
            obj.surfaceIntegrator = pbrtSurfaceIntegratorObject(); 

            % Renderer
            obj.renderer = 'sample';
            
            % World Begin
            
            %  Attribute Begin
            %   Light source
            obj.lightSourceArray = cell(1,1);
            whiteLight = pbrtLightSpotObject();
            obj.lightSourceArray{1} = whiteLight;
            %  Attribute End
            
            %  Material file
            obj.materialArray = cell(1,1);
            datapath = fullfile(s3dRootPath,'data');
            obj.materialArray{1} = fullfile(datapath,'validate', 'pbrtObject', 'depthTargetSpheres-mat.pbrt');
            
            %example materials object
            %             tempProperty = pbrtPropertyObject('color Kd', [0 0.374624 0]);  %TODO: the user shouldn't need to know pbrt syntax...
            %             greenLambertian = pbrtMaterialObject('greenLambertian', 'matte', tempProperty);
            %             obj.addMaterial(greenLambertian);
             
            % Geometry file
            obj.geometryArray = cell(1,1);
            obj.geometryArray{1} = fullfile(datapath,'validate', 'pbrtObject', 'depthTargetSpheres-geom.pbrt');
            %             examplePlane = pbrtGeometryObject();
            %             obj.addGeometry = examplePlane;
            
            % Shapes
%             obj.shapeArray = cell(1,1);
%             obj.shapeArray{1}.shape = pbrtShapeObject();
%             obj.shapeArray{1}.transform = pbrtTransformObject('translate', [4 0 3]);
            % WorldEnd
        end
        
       
        function returnVal = writeFile(obj, fname)
        %writeFile(obj, fname)
        %Write a text file from a pbrt structure
            fid = fopen(fname,'w');
            fprintf(fid,'# PBRT v2.0 Blender Scene file (written from Scene3D)\n#\n\n');
            
            %% Scale
            % s = pbrtGet(pbrt,'Scale');
            fprintf(fid,'Scale ');
            fprintf(fid,'%f  ',obj.scale);
            fprintf(fid,'# account for fixed lookat bug\n');
            
            %% Camera position & Lens & Image resolution
            obj.camera.writeFile(fid);

            %% Sampler
            obj.sampler.writeFile(fid);
            
            %% SurfaceIntegrator
            obj.surfaceIntegrator.writeFile(fid);
            
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
                    if (isa(obj.lightSourceArray{i}, 'pbrtLightObject'))
                        fprintf(fid,'\n\nAttributeBegin\n');
                        obj.lightSourceArray{i}.writeFile(fid);
                        fprintf(fid,'\nAttributeEnd\n');
                    else
                        error('Error! Non-light type placed in light array!');
                    end
                end
            end
            
            %% Materials File
            for i = 1:length(obj.materialArray)
                if (isa(obj.materialArray{i},'pbrtMaterialObject')); 
                    obj.materialArray{i}.writeFile(fid);
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
            for i = 1:length(obj.geometryArray)
                curGeometry =obj.geometryArray{i}; 
                if (isa(curGeometry,'pbrtGeometryObject')); 
                    curGeometry.writeFile(fid);
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
        %addLightSource(obj,newLightSource)
        %
        %adds a light source to the light source array
        %
        %newLightSource: must be of type pbrtLightObject 
        %TODO: error checking
            obj.lightSourceArray{end+1} = newLightSource; 
        end

        function addMaterial(obj, newMaterial)
        %addMaterial(obj, newMaterial)
        %
        %adds a material to the material array
        %
        %newMateria: must be of type pbrtMaterialObject
        %TODO: error checking
            obj.materialArray{end+1} = newMaterial;
        end
            
        function addGeometry(obj, newGeometry)
        %addGeometry(obj, newGeometry)
        %
        %%adds geoemtery to the geometry array
        %
        %newGeometry: must be of type pbrtGeometryObject
        %TODO: error checking     
            obj.geometryArray{end+1} = newGeometry;
        end        
        
%         function addShape(obj, newShape, newTransform)
%         %addShape(obj, newShape, newTransform)
%         %
%         %%adds shape to the shape array
%         %
%         %newShape: must be of type pbrtGeometryObject
%         %newTransform: must be of type pbrtTransformObject
%         %TODO: error checking     
%             obj.shapeArray{end+1} = newShape;
%             obj.shapeArray{end+1} = newTransform;
%         end    
        
        function removeGeometry(obj, deleteIndex)
        %removeGeometry(obj, deleteIndex)
        %
        %removes the geometry corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %TODO: error checking
            if (ieNotDefined('deleteIndex'))
                obj.geometryArray(end) = [];
            else
                obj.geometryArray(deleteIndex) = [];
            end
        end

        function returnVal = removeLight(obj, deleteIndex)
        %removeLight(obj, deleteIndex)
        %
        %removes the light source corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %returns the deleted value
            if (ieNotDefined('deleteIndex'))
                returnVal = obj.lightSourceArray(end);
                obj.lightSourceArray(end) = [];
            else
                returnVal = obj.lightSourceArray(deleteIndex);
                obj.lightSourceArray(deleteIndex) = [];
            end
        end
        
        function returnVal = removeMaterial(obj, deleteIndex)
        %removeMaterial(obj, deleteIndex)
        %
        %removes the light source corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %returns the deleted value
            if (ieNotDefined('deleteIndex'))
                returnVal = obj.materialArray(end);
                obj.materialArray(end) = [];
            else
                returnVal = obj.materialArray(deleteIndex);
                obj.materialArray(deleteIndex) = [];
            end
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

