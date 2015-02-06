% this is the simplest lens - no aperture is required because pinhole cameras
% have no aperture and the pinhole lens will inherit this superclas.
% This will be a superclass that will be inherited by other classes in the future
classdef pbrtLensRealisticObject < pbrtLensObject
    properties
        %filmDistance;
        %filmDiag;
        apertureDiameter;  %in mm
        diffractionEnabled = false;
        chromaticAberrationEnabled = false;
        specFile = '';
        curveRadius = 0;
    end
    methods
        
        %default constructor
        function obj = pbrtLensRealisticObject(inFilmDistance, inFilmDiag, ...
                specFile, apertureDiameter, diffractionEnabled, ...
                chromaticAberrationEnabled, curveRadius)
            

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
            
            %add additional properties
            if (ieNotDefined('specFile'))
                obj.specFile = 'default.pbrt';
            else
                obj.specFile = specFile;
            end
            
            if (ieNotDefined('apertureDiameter'))
                obj.apertureDiameter = 1;
            else
                obj.apertureDiameter = apertureDiameter;
            end
            
            if (ieNotDefined('diffractionEnabled'))
                obj.diffractionEnabled = false;
            else
                obj.diffractionEnabled = diffractionEnabled;
            end
            
            if (ieNotDefined('chromaticAberrationEnabled'))
                obj.chromaticAberrationEnabled = false;
            else
                obj.chromaticAberrationEnabled = chromaticAberrationEnabled;
            end         
            
            if (ieNotDefined('curveRadius'))
                obj.curveRadius = false;
            else
                obj.curveRadius = curveRadius;
            end               
        end
        
        
        %writes the section of text corresponding to this object
        function returnVal = writeFile(obj, fid)
            
             %writes the section of text corresponding to this object
            fprintf(fid,'\nCamera "realisticDiffraction"\n');
            fprintf(fid,'\t"float filmdistance" %f\n', obj.filmDistance);
            fprintf(fid,'\t"float filmdiag" %f\n', obj.filmDiag);
            fprintf(fid,'\t"float aperture_diameter" %f\n', obj.apertureDiameter);    
            fprintf(fid,'\t"float diffractionEnabled" %f\n', obj.diffractionEnabled);    
            fprintf(fid,'\t"float chromaticAberrationEnabled" %f\n', obj.chromaticAberrationEnabled); 
            fprintf(fid,'\t"float curveRadius" %f\n', obj.curveRadius); 
            fprintf(fid,'\t"string specfile" "%s"', obj.specFile); 
            
            returnVal = 1;
        end
    end
end