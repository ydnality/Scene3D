% arealightObject contains the subclass that makes pbrt arealightObject
classdef lightAreaObject <  lightObject
    
    properties (SetAccess = private)
        transformArray;   
        shapeArray;   %TODO: make a shape object
    end
    methods
        
        %default constructor.  The input variables may be omitted or left
        %with empty arguments if the user does not wish to specify them.  A
        %default value will be assumed.  
        function obj = lightAreaObject(inName, inSpectrum)
            %superclass properties
            obj@lightObject();
            obj.setType('diffuse');  %this is the only option for area lights right now in pbrt
            
            if (~ieNotDefined('inName'))
                obj.setName(inName);
            end
            if (~ieNotDefined('inSpectrum'))
                obj.setSpectrum(inSpectrum);
            end
            
            %obj.transformArray = cell(1,1);
            obj.addTransform(transformObject());
            %obj.shapeArray = cell(1,1);
            obj.addShape(shapeObject('sphere', 'radius', 1));  %default shape is a tiny sphere
        end

        %adds a transform to the transformArray
        function addTransform(obj,inTransform)
            validateattributes(inTransform, {'transformObject'}, {'nonempty'});
            obj.transformArray{length(obj.transformArray) + 1} = inTransform;
        end
        
        %adds a shape to the shapeArray
        function addShape(obj, inShape)
            validateattributes(inShape, {'shapeObject'}, {'nonempty'});
            obj.shapeArray{length(obj.shapeArray) + 1} = inShape;
        end       

        %removes the shape corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %returns the deleted value
        function returnVal = removeShape(obj, removeIndex)
            if (ieNotDefined('removeIndex'))
                returnVal = obj.shapeArray{end};
                obj.shapeArray(end)= [];
            else
                validateattributes(inShape, {'shapeObject'}, {'nonempty'});
                returnVal = obj.shapeArray{removeIndex};
                obj.shapeArray(removeIndex)= [];
            end
        end
        
        %removes the transform source corresponding to the specified index
        %if deleteIndex is undefined, remove from the end
        %returns the deleted value
        function returnVal = removeTransform(obj, removeIndex)
            if (ieNotDefined('removeIndex'))
                returnVal = obj.transformArray{end};
                obj.transformArray(end)= [];
            else
                validateattributes(inShape, {'shapeObject'}, {'nonempty'});
                returnVal = obj.transformArray{removeIndex};
                obj.transformArray(removeIndex)= [];
            end
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
        %TODO: check fid error checking
        function writeFile(obj, fid)
           %call superclass function 
%            writeFile@lightObject(obj, fid);
           
           fprintf(fid,'\n\tAreaLightSource "%s" "%s" [', obj.type,obj.spectrum.type);
           fprintf(fid,'%f ', obj.spectrum.value);
           fprintf(fid,']\n');
           
           %write transforms
           for i = 1:length(obj.transformArray)
               fprintf(fid,'\t');
               obj.transformArray{i}.writeFile(fid);
           end
           
           %write shapes
           for i = 1:length(obj.shapeArray)
                obj.shapeArray{i}.writeFile(fid);
           end           
        end
    end

end