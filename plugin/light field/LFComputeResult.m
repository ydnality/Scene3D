% Compute several interesting features for a given LF 

function [result]=LFComputeResult(LF)


%INPUT
%LFobj: it an object describing a light light field 


%OUTPUT
%result: object with several field describing the relevant features



%% 
for li=1:size(LF.wave,1)
    % Parameter of eccentricity
    maxY=max(max(LF.Y(:,:,li)));
    minY=min(min(LF.Y(:,:,li)));
    DeltaY(li,1)=maxY-minY;
    
    % Parameter of incident angles
    maxU=max(max(LF.U(:,:,li)));
    minU=min(min(LF.U(:,:,li)));
    DeltaU(li,1)=maxU-minU;
end


%% OTHERs COULD BE ADDED ..........



%% Append field to the output

result.spotSize=DeltaY;
result.spotAperture=DeltaU;
