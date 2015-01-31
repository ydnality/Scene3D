% spotlightObject contains the subclass that makes pbrt spotlights
% under construction.  Does not work yet!
classdef pbrtLightInfiniteObject <  pbrtLightObject
    
    properties (SetAccess = private)
        nSamples; 
        mapName; 
        transformArray;
        %additional properties may be specified within propertyArray
    end
    methods
        
        function obj = pbrtLightInfiniteObject(inName, inNSamples, inSpectrum, inMapName, inTransform)
            %obj = pbrtLightInfiniteObject(inName, inNSamples, scale, inSpectrum, inMapName, inTransform)
            %default constructor.  The input variables may be omitted or left
            %with empty arguments if the user does not wish to specify them.  A
            %default value will be assumed.
        
            %superclass properties
            obj@pbrtLightObject();
            obj.setType('infinite');
            
            if (~ieNotDefined('inName'))
                obj.setName(inName);
            end
            
            %infinite light specific properties
            if(ieNotDefined('inSpectrum'))
                obj.setSpectrum(pbrtSpectrumObject('color L', [1 1 1]));
            else
                obj.setSpectrum(inSpectrum);
            end
            
            if(ieNotDefined('inNSamples'))
                obj.setNSamples(8);
            else
                obj.setDeltaAngle(inNSamples);
            end
            if(~ieNotDefined('inMapName'))
                obj.setMapName(inMapName);
            end
            if(~ieNotDefined('inTransform'))
                obj.addTransform(inTransform);
            end            
        end

        function setNSamples(obj, inNSamples)
           obj.nSamples = inNSamples; 
           return;
        end    
        
%         function setSpectrum(obj, inSpectrum)
%            obj.spectrum = inSpectrum; 
%            return;
%         end           
        
        function setMapName(obj, inMapName)
           obj.mapName = inMapName; 
           return;
        end      
 
        function addTransform(obj, inTransform)
            %addTransform(obj, inTransform)
            %
            %%adds transform to the end of transform array
            %
            %inTransform: must be of type pbrtTransformObject
            assert(isa(inTransform, 'pbrtTransformObject'), 'pbrtTransformObject required as input!');
            obj.materialArray{end+1} = inTransform;
            return;
        end
        
        function returnVal = removeTransform(obj, deleteIndex)
            %removeTransform(obj, deleteIndex)
            %
            %removes the transformcorresponding to the specified index
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
                
        
        %writes the pbrt file corresponding to this object
        function writeFile(obj, fid)           
            
           %fprintf(fid, '\t\tAttribute Begin\n');
           for i = 1:length(obj.transformArray)
               obj.transformArray(i).writeFile(fid);
           end

           fprintf(fid,'\n\tLightSource\n');
           fprintf(fid,'\t\t"infinite"\n');
           fprintf(fid,'\t\t"integer nsamples" [%i]', obj.nSamples);
           
           if(~isempty(obj.mapName))
                fprintf(fid,'\t\t"string mapname" ["%s"]', obj.mapName);
           end
           fprintf(fid, '\n');
          % fprintf(fid, '\t\tAttribute End\n');
        end
    end

end