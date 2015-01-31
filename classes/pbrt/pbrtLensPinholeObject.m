% this is a pinhole "lens object.  The only 2 parameters are film distance
% and film diagonal. 
classdef pbrtLensPinholeObject <  pbrtLensObject
    properties 
        %additional properties go here
    end
    methods
        
        %default constructor
        function obj = pbrtLensPinholeObject(inFilmDistance, inFilmDiag)
%         function obj = lensPinholeObject(inFilmDistance, inFilmDiag)
            %call the superclass lensObject constructor
            
            
            if (ieNotDefined('inFilmDistance'))
                % Example lens
                obj.filmDistance = 140;
            else
                obj.filmDistance = inFilmDistance;
            end
            if (ieNotDefined('inFilmDiag'))
                % Example lens
                obj.filmDiag = 43.267;
            else
                obj.filmDiag = inFilmDiag;
            end
%             obj = obj@pbrtLensObject(inFilmDistance, inFilmDiag);


%             obj = obj@pbrtLensObject(varargin);
        end
        
        
        function returnVal = writeFile(obj, fid)
            %writes the section of text corresponding to this object
            fprintf(fid,'\nCamera "pinhole"\n');
            fprintf(fid,'\t"float filmdistance" %f\n', obj.filmDistance);
            fprintf(fid,'\t"float filmdiag" %f\n', obj.filmDiag);
            returnVal = 1;
        end        
    end
end