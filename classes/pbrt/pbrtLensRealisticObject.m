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
        pinholeExitApLoc;
        filmCenter = [0 0]; %because of a pbrt oddity, we had to put this here although it doesn't make intuitive sense
        numPinholesW = -1;   %similar oddity here.  These 2 properties are for a pinhole array
        numPinholesH = -1;  
        microlensMode = false;  %flag for microlenses
    end
    methods
        
        %default constructor
        function obj = pbrtLensRealisticObject(inFilmDistance, inFilmDiag, ...
                specFile, apertureDiameter, diffractionEnabled, ...
                chromaticAberrationEnabled, curveRadius, pinholeExitApLoc, filmCenter, ...
                numPinholesW, numPinholesH, microlensMode)
            

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
            if(~ieNotDefined('pinholeExitApLoc'))
               obj.pinholeExitApLoc  = pinholeExitApLoc;  
            end
            if(~ieNotDefined('filmCenter'))
               obj.filmCenter  = filmCenter;  
            end
            if(~ieNotDefined('numPinholesW'))
                obj.numPinholesW = numPinholesW;
            end
            if(~ieNotDefined('numPinholesH'))
                obj.numPinholesH = numPinholesH;
            end
            
            if(~ieNotDefined('microlensMode'))
                obj.microlensMode = microlensMode;
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
            fprintf(fid,'\t"string specfile" "%s"\n', obj.specFile); 
            fprintf(fid,'\t"float film_center_x" %f\n', obj.filmCenter(1));
            fprintf(fid,'\t"float film_center_y" %f\n', obj.filmCenter(2));
            
            if (~isempty(obj.pinholeExitApLoc))
                fprintf(fid,'\t"float pinhole_exit_x" %f\n', obj.pinholeExitApLoc(1));
                fprintf(fid,'\t"float pinhole_exit_y" %f\n', obj.pinholeExitApLoc(2));
                fprintf(fid,'\t"float pinhole_exit_z" %f\n', obj.pinholeExitApLoc(3));
            end
            
            if( obj.numPinholesH > 0 && obj.numPinholesW > 0)
                 fprintf(fid,'\t"float num_pinholes_w" %i\n', obj.numPinholesW);
                 fprintf(fid,'\t"float num_pinholes_h" %i\n', obj.numPinholesH);
            end
            
            if (obj.microlensMode)
                fprintf(fid,'\t"float microlens_enabled" 1\n');
            end
            
            returnVal = 1;
        end
    end
end