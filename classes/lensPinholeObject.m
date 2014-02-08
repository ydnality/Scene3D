% this is a pinhole "lens object.  The only 2 parameters are film distance
% and film diagonal. 
classdef lensPinholeObject <  lensObject
    properties 
        %additional properties go here
    end
    methods
        
        %default constructor
        function obj = lensPinholeObject(varargin)
%         function obj = lensPinholeObject(inFilmDistance, inFilmDiag)
            %call the superclass lensObject constructor
%             obj = obj@lensObject(inFilmDistance, inFilmDiag);
            obj = obj@lensObject(varargin);
        end
        
        %writes the section of text corresponding to this object
        function returnVal = writeFile(obj, fid)
            fprintf(fid,'\nCamera "pinhole"\n');
            fprintf(fid,'\t"float filmdistance" %f\n', obj.filmDistance);
            fprintf(fid,'\t"float filmdiag" %f\n', obj.filmDiag);
            returnVal = 1;
        end        
    end
end