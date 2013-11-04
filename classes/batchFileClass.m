% This class helps with the management of file names for the batch render
% step.  We supply a file stem, and postfixes, which apend to the stem to
% create files with the same file base, but with different endings.  This
% class deals with "input" and "output' names.  The "input" name is the
% name of the template, while the "output" name are the generated files.
% These file names will share the same postfix, but have different stems.
% This class handles this relationship and makes it easier to manage.
%
% To use the object, first supply the stem, the input stem, then the post fix
% Next, call obj.getInputFile to obtain the input file names, and 
% obj.getOutputFile to obtain the appropriate output file name.

classdef batchFileClass <  handle
    
    %matlab.mixin.Copyable
    properties
        inputStem;
        outputStem;
        postfix; 
    end
    methods
        function obj = batchFileClass(inStem, inPostfix)
            obj.inputStem = inStem;
            obj.postfix = inPostfix;
            obj.outputStem = '';
        end
        
        function obj=setInputStem(obj, inStem)
            obj.inputStem = inStem;
            
        end
        
        function obj=setOutputStem(obj, outStem)
            obj.outputStem = outStem;
        end

        function inputFile = getInputFile(obj)
            inputFile = strcat(obj.inputStem, obj.postfix);
        end
        
        function outString = getOutputFile(obj)
            outString = strcat(obj.outputStem, obj.postfix);
        end
        
        function postFix = getPostFix(obj)
            postFix = obj.postfix;
        end
    end
    
    
end