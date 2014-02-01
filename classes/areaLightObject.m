% arealightObject contains the subclass that makes pbrt arealightObject
classdef areaLightObject <  lightObject
    
    properties (SetAccess = private)
        transformArray;   
        shapeArray;   %TODO: make a shape object
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = areaLightObject(inName, inSpectrum)
            %superclass properties
            obj@lightObject();
            obj.setType('arealight');
            
            if (~ieNotDefined('inName'))
                obj.setName(inName);
            end
            if (~ieNotDefined('inSpectrum'))
                obj.setSpectrum(inSpectrum);
            end
            
            obj.transformArray = cell(1,1);
            obj.addTransform(transformObject());
            obj.shapeArray = cell(1,1);
        
        end

        function addTransform(obj,inTransform)
            validateattributes(inTransform, {'transformObject'}, {'nonempty'});
            obj.transformArray{length(obj.transformArray) + 1} = inTransform;
        end
        
        function addShape(obj, inShape)
            validateattributes(inShape, {'shapeObject'}, {'nonempty'});
            obj.shapeArray{length(obj.shapeArray) + 1} = inShape;
        end       

        %sets the direction of the camera
        %inFrom: how much the light should be moved. Note that the
        %direction stays the same
        function move(obj, offset)
            %TODO: error check
            obj.addTransform(obj, transformObject('translate', offset));
           return;
        end        
        
        %writes the pbrt file corresponding to this object
        function writeFile(obj, fid)
           %call superclass function 
%            writeFile@lightObject(obj, fid);
           
           fprintf(fid,'\n\tAreaLightSource "%s" "%s" [', obj.type,obj.spectrum.type);
           fprintf(fid,'%f]\n', obj.spectrum.value);
           
           %write transforms
           for i = 1:length(obj.transformArray)
               fprintf(fid,'\t\t');
               obj.transformArray{i}.writeFile();
           end
           
           %write shapes
           for i = 1:length(obj.shapeArray)
                obj.shapeArray{i}.writeFile();
           end           
        end
    end

end