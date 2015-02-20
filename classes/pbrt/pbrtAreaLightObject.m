% pbrtAreaLightObject contains the subclass that makes pbrt pbrtAreaLightObject
classdef pbrtAreaLightObject <  pbrtLightObject
    
    properties (SetAccess = private)
        transformArray;   
        shapeArray;   %TODO: make a shape object
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = pbrtAreaLightObject(inName, inSpectrum)
            %superclass properties
            obj@pbrtLightObject();
            obj.setType('area');
            
            if (~ieNotDefined('inName'))
                obj.setName(inName);
            end
            if (~ieNotDefined('inSpectrum'))
                obj.setSpectrum(inSpectrum);
            end
            
            obj.transformArray = cell(0,1);
            %obj.addTransform(pbrtTransformObject());
            obj.shapeArray = cell(0,1);
        
        end

        function addTransform(obj,inTransform)
            validateattributes(inTransform, {'pbrtTransformObject'}, {'nonempty'});
            obj.transformArray{length(obj.transformArray) + 1} = inTransform;
        end
        
        function addShape(obj, inShape)
            validateattributes(inShape, {'pbrtShapeObject'}, {'nonempty'});
            obj.shapeArray{length(obj.shapeArray) + 1} = inShape;
        end       

        %sets the direction of the camera
        %inFrom: how much the light should be moved. Note that the
        %direction stays the same
        function move(obj, offset)
            %TODO: error check
            obj.addTransform(obj, pbrtTransformObject('translate', offset));
           return;
        end        
        
        %writes the pbrt file corresponding to this object
        function writeFile(obj, fid)
           %call superclass function 
%            writeFile@lightObject(obj, fid);
           
           fprintf(fid,'\n\tAreaLightSource "%s" "%s" [', obj.type,obj.spectrum.type);
           fprintf(fid,'%f ', obj.spectrum.value);
           fprintf(fid, ']\n');
           
           %write transforms
           for i = 1:length(obj.transformArray)
               fprintf(fid,'\t\t');
               obj.transformArray{i}.writeFile(fid);
           end
           
           %write shapes
           for i = 1:length(obj.shapeArray)
                obj.shapeArray{i}.writeFile(fid);
           end
           
           %write properties
           for i = 1:length(obj.propertyArray)
               obj.propertyArray{i}.writeFile(fid);
           end
        end
    end

end