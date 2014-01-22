% lightObject contains the object to create a light in PBRT
classdef lightObject <  handle
    
    properties %(SetAccess = private)
        name;
        spectrum; 
        coneAngle;
        coneDeltaAngle;
        from;
        to;
        
        % declare an additional spotlight
%         light.type = 'light';
%         light.lightType = 'spot';
%         light.spectrum.type = 'rgb I';
%         light.spectrum.value = [1000 1000 1000];
%         light.coneangle = 180;
%         light.conedeltaangle = 180;
%         light.from = [4.5 -90 8.5];
%         light.to = [4.5 -89 8.5];

    end
    methods
        
        %default constructor
        function obj = lightObject()
            obj.name = 'defaultLight';
            obj.spectrum.type = 'rgb I';
            obj.spectrum.value = [1000 1000 1000];
            obj.coneAngle = 180;
            obj.coneDeltaAngle = 180;
            obj.from =  [4.5 -90 8.5];
            obj.to = [4.5 -89 8.5];
        end
        
        
    end
    
end